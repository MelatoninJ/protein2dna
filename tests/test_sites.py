"""Choosing two individual restriction sites, and the numbers that follow from the choice.

The size and digest arithmetic is checked against Biopython's own cut positions rather
than against the code under test, so a wrong cut offset cannot hide.
"""

from __future__ import annotations

import copy
import os
from pathlib import Path

import pytest
from Bio.Restriction import RestrictionBatch
from Bio.Seq import Seq

import p2d
from p2d import vectorio
from p2d.analysis import analyze_vector, confirm
from p2d.enzymes import (
    backbone_conflicts,
    find_cut_site,
    list_cut_sites,
    plan_cut,
)
from p2d.report import construct_report
from p2d.vector import StartSource, Vector

SAMPLES = Path(os.environ["P2D_VECTOR_SAMPLES"]) if os.environ.get("P2D_VECTOR_SAMPLES") else None
PAYLOAD = "SKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLKFICTTGKLPVPWPTL"


def sample(name: str) -> Path:
    if SAMPLES is None or not (SAMPLES / name).is_file():
        pytest.skip(f"{name} not available (set P2D_VECTOR_SAMPLES)")
    return SAMPLES / name


@pytest.fixture(scope="module")
def vector():
    return p2d.load_bundled()


@pytest.fixture(scope="module")
def host():
    return p2d.load_host("ecoli_bl21")


def design(vector, host, up, down, payload=PAYLOAD, **kw):
    site = plan_cut(vector, find_cut_site(vector, up), find_cut_site(vector, down))
    plan = p2d.ConstructPlan.single_chain(payload)
    return p2d.assemble(p2d.design_insert(plan, vector, site, host, seed=1, **kw))


def bio_fragments(seq: str, names: list[str]) -> list[int]:
    """Circular digest lengths from Biopython's cut positions (1-based, first base after cut)."""
    cuts = sorted(
        {
            p - 1
            for hits in RestrictionBatch(names).search(Seq(seq), linear=False).values()
            for p in hits
        }
    )
    if not cuts:
        return [len(seq)]
    n = len(seq)
    return sorted(((cuts[(i + 1) % len(cuts)] - cuts[i]) % n) or n for i in range(len(cuts)))


def with_second_site(vector: Vector, enzyme_site: str, at: int) -> Vector:
    """A copy of the vector with a same-length second recognition site written at ``at``.

    Same length, so the annotations stay valid and are kept.
    """
    record = copy.deepcopy(vector.record)
    record.seq = Seq(vector.seq[:at] + enzyme_site + vector.seq[at + len(enzyme_site) :])
    v = Vector(record=record, name="multi")
    v.expression_start = vector.expression_start
    return v


def test_every_site_is_listed_with_its_position_and_count(vector):
    sites = list_cut_sites(vector, (80, 200))
    assert [(c.enzyme, c.start) for c in sites] == [
        ("NcoI", 85),
        ("NdeI", 144),
        ("BamHI", 150),
        ("Acc65I", 156),
        ("KpnI", 156),
        ("EcoRI", 168),
        ("XhoI", 174),
    ]
    xho = sites[-1]
    assert xho.id == "XhoI:174" and xho.label == "XhoI at 175" and xho.count == 1
    assert xho.cut_offset == 1 and xho.overhang == "TCGA" and "SlaI" in xho.aliases


def test_blunt_type_iis_and_short_sites_are_left_out(vector):
    names = {c.enzyme for c in list_cut_sites(vector)}
    assert not names & {"SmaI", "BsaI"}
    assert all(len(c.site) >= 6 for c in list_cut_sites(vector))


def test_a_site_in_several_places_is_numbered(vector):
    twice = with_second_site(vector, "CTCGAG", 400)
    xho = [c for c in list_cut_sites(twice) if c.enzyme == "XhoI"]
    assert [(c.start, c.number, c.count) for c in xho] == [(174, 1, 2), (400, 2, 2)]
    assert xho[1].label == "XhoI at 401 (2 of 2)"


def test_window_keeps_only_sites_wholly_inside(vector):
    assert {c.enzyme for c in list_cut_sites(vector, (150, 156))} == {"BamHI"}
    assert list_cut_sites(vector, (151, 156)) == []


def test_site_ids_resolve_by_alias_and_reject_nonsense(vector):
    assert find_cut_site(vector, "SlaI:174").enzyme == "XhoI"
    for bad in ["XhoI", "XhoI:abc", ":174", "XhoI:175", "Nope:174"]:
        with pytest.raises(ValueError):
            find_cut_site(vector, bad)


@pytest.mark.parametrize(
    "up, down, message",
    [
        ("XhoI:174", "NdeI:144", "must lie upstream"),
        ("BamHI:150", "BamHI:150", "must lie upstream"),
    ],
)
def test_impossible_pairs_say_why(vector, up, down, message):
    with pytest.raises(ValueError, match=message):
        plan_cut(vector, find_cut_site(vector, up), find_cut_site(vector, down))


def test_a_start_codon_is_required():
    seq = p2d.load_bundled().seq
    bare = Vector.from_sequence(seq, name="bare")
    with pytest.raises(ValueError, match="no confirmed start codon"):
        plan_cut(bare, find_cut_site(bare, "NdeI:144"), find_cut_site(bare, "XhoI:174"))


def test_a_site_before_the_start_codon_is_refused(vector):
    moved = Vector.from_sequence(vector.seq, name="moved")
    moved.expression_start = 160  # not an ATG context, but the check order is what matters
    with pytest.raises(ValueError, match="lies before the start codon"):
        plan_cut(moved, find_cut_site(moved, "BamHI:150"), find_cut_site(moved, "XhoI:174"))


def test_the_site_holding_the_start_codon_may_be_the_5prime_site(vector):
    site = plan_cut(vector, find_cut_site(vector, "NcoI:85"), find_cut_site(vector, "XhoI:174"))
    assert site.start_source is StartSource.UPSTREAM_SITE and site.start_offset_in_site == 2
    assert (site.upstream_pos, site.downstream_pos) == (85, 174)


def test_a_stop_codon_between_start_and_site_is_refused(vector):
    seq = vector.seq[:120] + "TAA" + vector.seq[123:]  # a stop in frame inside the His tag
    v = Vector.from_sequence(seq, name="stopped")
    v.expression_start = 87
    with pytest.raises(ValueError, match="stop codon at 121"):
        plan_cut(v, find_cut_site(v, "BamHI:150"), find_cut_site(v, "XhoI:174"))


def test_a_second_site_inside_the_replaced_piece_is_harmless(vector, host):
    """XhoI also cuts at 162, between BamHI and the chosen XhoI at 174: it is replaced."""
    v = with_second_site(vector, "CTCGAG", 162)
    assert [c.start for c in list_cut_sites(v) if c.enzyme == "XhoI"] == [162, 174]
    a = design(v, host, "BamHI:150", "XhoI:174")
    assert a.ok, a.errors
    final = Vector.from_sequence(a.seq, name="final")
    assert len(final.find_all("XhoI")) == 1, "the extra site left with the replaced piece"
    real = str(Seq(a.seq[87:][: 3 * len(a.protein)]).translate())
    assert real == a.protein


def test_a_second_site_in_the_kept_backbone_is_refused(vector):
    """Choosing the XhoI at 162 leaves the one at 174 in the backbone: a digest cuts there too."""
    v = with_second_site(vector, "CTCGAG", 162)
    with pytest.raises(ValueError, match=r"also cut the backbone you keep \(XhoI at 175\)"):
        plan_cut(v, find_cut_site(v, "BamHI:150"), find_cut_site(v, "XhoI:162"))


def test_a_distant_second_site_is_refused_for_either_enzyme(vector):
    far = with_second_site(vector, "GGATCC", 700)  # a second BamHI, far away
    with pytest.raises(ValueError, match=r"BamHI at 701"):
        plan_cut(far, find_cut_site(far, "BamHI:150"), find_cut_site(far, "XhoI:174"))
    far = with_second_site(vector, "CTCGAG", 700)
    with pytest.raises(ValueError, match=r"XhoI at 701"):
        plan_cut(far, find_cut_site(far, "BamHI:150"), find_cut_site(far, "XhoI:174"))
    # an enzyme the pair does not use does not matter
    far = with_second_site(vector, "GAATTC", 700)
    assert plan_cut(far, find_cut_site(far, "BamHI:150"), find_cut_site(far, "XhoI:174"))


def test_conflicts_name_every_offender(vector):
    v = with_second_site(with_second_site(vector, "CTCGAG", 700), "GGATCC", 800)
    up, down = find_cut_site(v, "BamHI:150"), find_cut_site(v, "XhoI:174")
    assert backbone_conflicts(v, up, down) == [("XhoI", 700), ("BamHI", 800)]


PAIRS = [
    ("NdeI:144", "XhoI:174"),
    ("BamHI:150", "XhoI:174"),
    ("NdeI:144", "EcoRI:168"),
    ("NcoI:85", "XhoI:174"),
]


@pytest.mark.parametrize("up, down", PAIRS)
def test_sizes_add_up_and_match_biopython(vector, host, up, down):
    a = design(vector, host, up, down)
    assert a.ok, a.errors
    r = construct_report(a)
    s = r["sizes"]
    assert s["backbone_bp"] + s["replaced_bp"] == s["vector_bp"] == 1149
    assert s["backbone_bp"] + s["insert_after_digest_bp"] == s["plasmid_bp"] == len(a.seq)
    assert s["plasmid_bp"] - s["vector_bp"] == s["net_change_bp"]
    names = r["digest_check"]["enzymes"]
    assert sorted(f["length"] for f in r["digest_check"]["vector"]) == bio_fragments(
        vector.seq, names
    )
    assert sorted(f["length"] for f in r["digest_check"]["final"]) == bio_fragments(a.seq, names)


def test_digest_roles_name_the_backbone_and_the_insert(vector, host):
    r = construct_report(design(vector, host, "BamHI:150", "XhoI:174"))
    roles = {f["role"]: f["length"] for f in r["digest_check"]["vector"]}
    assert roles == {"backbone (kept)": 1125, "replaced piece": 24}
    final = {f["role"]: f["length"] for f in r["digest_check"]["final"]}
    assert (
        final["backbone (kept)"] == 1125 and final["insert"] == r["sizes"]["insert_after_digest_bp"]
    )


def test_the_insert_parts_rebuild_the_insert(vector, host):
    a = design(vector, host, "NdeI:144", "EcoRI:168")
    parts = construct_report(a)["insert_parts"]
    assert sum(len(p["dna"]) if "dna" in p else p["bp"] for p in parts) == len(a.insert_dna)
    assert parts[0]["dna"] == "CATATG" and parts[-1]["dna"] == "GAATTC"


def test_frame_pads_are_reported_when_the_sites_need_them(vector, host):
    shifted = Vector.from_sequence(vector.seq[:162] + "A" + vector.seq[162:], name="shift1")
    shifted.expression_start = 87
    a = design(
        shifted,
        host,
        "BamHI:150",
        next(c.id for c in list_cut_sites(shifted) if c.enzyme == "XhoI"),
    )
    r = construct_report(a)
    assert r["frame"]["pad_3prime"] == a.design.pad_down != ""
    assert any(p["name"] == "frame pad" for p in r["insert_parts"])


def test_protein_numbers_are_sane(vector, host):
    r = construct_report(design(vector, host, "NdeI:144", "XhoI:174"))["protein"]
    assert r["payload_residues"] == len(PAYLOAD) and r["residues"] > r["payload_residues"]
    assert 6 < r["mass_kda"] < 12 and 4 < r["pi"] < 10


def test_replaced_features_warn_when_something_important_goes(vector, host):
    """Choose a far second XhoI and the replaced stretch swallows the origin of replication."""
    v = with_second_site(vector, "CTCGAG", 800)
    a = design(v, host, "BamHI:150", "XhoI:800")
    gone = {
        f["label"]
        for f in construct_report(a)["replaced_features"]
        if f["fully_removed"] and f["essential"]
    }
    assert any("ori" in label for label in gone)
    near = construct_report(design(vector, host, "BamHI:150", "XhoI:174"))["replaced_features"]
    assert not [f for f in near if f["essential"]], "a normal MCS swap removes nothing essential"


def real(name):
    return confirm(analyze_vector(vectorio.read_vector_file(sample(name))))


def test_real_pet28a_every_valid_pair_matches_biopython(host):
    v = real("pET-28a_plus.dna")
    mcs = list_cut_sites(v, (v.expression_start, 5260))
    tried = 0
    for up in mcs:
        for down in mcs:
            try:
                site = plan_cut(v, up, down)
            except ValueError:
                continue
            a = p2d.assemble(
                p2d.design_insert(p2d.ConstructPlan.single_chain(PAYLOAD), v, site, host, seed=3)
            )
            assert a.ok, (up.label, down.label, a.errors)
            real_protein = str(Seq(a.seq[v.expression_start :][: 3 * len(a.protein)]).translate())
            assert real_protein == a.protein, (up.label, down.label)
            r = construct_report(a)
            names = r["digest_check"]["enzymes"]
            assert sorted(f["length"] for f in r["digest_check"]["final"]) == bio_fragments(
                a.seq, names
            )
            assert r["sizes"]["plasmid_bp"] == len(a.seq)
            tried += 1
    assert tried >= 20


def test_real_pet28a_known_products(host):
    v = real("pET-28a_plus.dna")
    tagged = design(v, host, "NdeI:5127", "XhoI:5206")
    assert tagged.protein.startswith("MGSSHHHHHHSSGLVPRGSHM" + PAYLOAD)
    bam = design(v, host, "BamHI:5166", "XhoI:5206")
    assert bam.protein.startswith("MGSSHHHHHHSSGLVPRGSHMASMTGGQQMGRGS" + PAYLOAD)


def test_real_pet28a_c_terminus_is_a_choice(host):
    """XhoI and the His6 tag are one base off the start codon's frame in pET-28a(+)."""
    v = real("pET-28a_plus.dna")
    keep = design(v, host, "BamHI:5166", "XhoI:5206")  # default: the vector's own frame
    # the pad is the vector's own base, so the residue is the one the vector itself reads here
    assert keep.protein.endswith(PAYLOAD + "TRAPPPPPLRSGC"), "reads vector junk, not the tag"
    tag = design(v, host, "BamHI:5166", "XhoI:5206", cterm=p2d.CTerm.TAG)
    assert tag.ok and tag.protein.endswith(PAYLOAD + "LEHHHHHH")
    # the tag follows XhoI in its own frame, so reaching it needs no extra base; staying
    # in the start codon's frame is what costs one (and shifts the tag out of register)
    assert tag.design.pad_down == "" and len(keep.design.pad_down) == 1
    stop = design(v, host, "BamHI:5166", "XhoI:5206", cterm=p2d.CTerm.STOP)
    assert stop.ok and stop.protein.endswith(PAYLOAD) and stop.design.stop_in_insert
    for a in (tag, stop):
        real_protein = str(Seq(a.seq[v.expression_start :][: 3 * len(a.protein)]).translate())
        assert real_protein == a.protein
        assert construct_report(a)["sizes"]["plasmid_bp"] == len(a.seq)


def test_tag_mode_reads_the_tag_in_frame_even_when_the_frame_is_shifted(vector, host):
    shifted = Vector.from_sequence(vector.seq[:162] + "A" + vector.seq[162:], name="shift1")
    shifted.expression_start = 87
    xho = next(c.id for c in list_cut_sites(shifted) if c.enzyme == "XhoI")
    default = design(shifted, host, "BamHI:150", xho)
    tag = design(shifted, host, "BamHI:150", xho, cterm=p2d.CTerm.TAG)
    assert tag.ok and tag.protein.endswith("HHHHHH")
    assert not default.protein.endswith("HHHHHH"), "without the choice the tag is out of frame"
    assert tag.design.pad_down != default.design.pad_down


def test_stop_mode_ends_the_protein_and_needs_no_downstream_pad(vector, host):
    a = design(vector, host, "NdeI:144", "XhoI:174", cterm=p2d.CTerm.STOP)
    assert a.ok and a.protein.endswith(PAYLOAD) and a.design.pad_down == ""
    assert a.design.stop_in_insert and a.design.stop_dna == "TAA"
    assert a.insert_dna.endswith("TAACTCGAG"), "one stop codon, then the 3' site"
    assert not design(vector, host, "NdeI:144", "XhoI:174").design.stop_in_insert


def test_two_stop_codons_are_available_on_request(vector, host):
    up, down = (find_cut_site(vector, "NdeI:144"), find_cut_site(vector, "XhoI:174"))
    site = plan_cut(vector, up, down)
    plan = p2d.ConstructPlan.single_chain(PAYLOAD)

    def build(n):
        return p2d.assemble(
            p2d.design_insert(plan, vector, site, host, cterm=p2d.CTerm.STOP, stop_codons=n)
        )

    one, two = build(1), build(2)
    assert one.ok and two.ok and two.protein == one.protein
    assert (one.design.stop_dna, two.design.stop_dna) == ("TAA", "TAATGA")
    assert len(two.insert_dna) == len(one.insert_dna) + 3
    with pytest.raises(ValueError, match="one or two"):
        build(3)


def test_the_stop_codons_are_the_hosts_most_used_in_order(host):
    assert host.stop_codons(1) == "TAA" and host.stop_codons(2) == "TAATGA"
    other = p2d.load_host("hsapiens")
    assert other.stop_codons(1) == max(other.stop_usage, key=other.stop_usage.get)


def test_tag_mode_without_a_tag_says_so(vector, host):
    seq = vector.seq[:180] + "GCCGCCGCCGCCGCCGCC" + vector.seq[198:]  # the His6 codons become Ala
    v = Vector.from_sequence(seq, name="notag")
    v.expression_start = 87
    with pytest.raises(ValueError, match="No C-terminal tag"):
        design(v, host, "NdeI:144", "XhoI:174", cterm=p2d.CTerm.TAG)


def test_find_cterm_tag_reports_where_the_tag_starts(vector):
    from p2d.frame import find_cterm_tag

    tag = find_cterm_tag(vector, 180)
    assert tag and tag.position == 180 and tag.residues == "HHHHHH"
    assert find_cterm_tag(vector, 300) is None


CODON = {"A": "GCT", "D": "GAT", "E": "GAA", "G": "GGT", "H": "CAT", "I": "ATT", "K": "AAA",
         "L": "CTG", "N": "AAT", "P": "CCG", "Q": "CAG", "S": "TCT", "T": "ACC", "V": "GTT",
         "W": "TGG", "Y": "TAT", "F": "TTT"}  # fmt: skip


def dna(aa):
    return "".join(CODON[c] for c in aa)


@pytest.mark.parametrize("frame", [0, 1, 2])
@pytest.mark.parametrize(
    "name,aa",
    [("His6", "HHHHHH"), ("FLAG", "DYKDDDDK"), ("Myc", "EQKLISEEDL"), ("Strep-II", "WSHPQFEK"),
     ("HA", "YPYDVPDYA"), ("V5", "GKPIPNPLLGLDST")],
)  # fmt: skip
def test_every_known_tag_is_found_in_any_frame(name, aa, frame):
    from p2d.frame import find_cterm_tag

    lead = "C" * frame + dna("LEAAA")
    seq = "ACGT" * 5 + lead + dna(aa) + "TAA" + "ACGT" * 20
    found = find_cterm_tag(Vector.from_sequence(seq, name="t"), 20)
    assert found is not None and found.name == name and found.residues == aa
    assert found.position == 20 + frame + len(dna("LEAAA"))


def test_the_earliest_tag_wins_and_a_stop_hides_what_follows():
    from p2d.frame import find_cterm_tag

    both = "A" * 10 + dna("GSDYKDDDDKGSHHHHHH") + "TAA" + "C" * 12
    assert find_cterm_tag(Vector.from_sequence(both, name="t"), 10).name == "FLAG"
    hidden = "A" * 10 + dna("GS") + "TAA" + dna("HHHHHH") + "TAA" + "C" * 12
    assert find_cterm_tag(Vector.from_sequence(hidden, name="t"), 10) is None


def test_a_his_run_is_reported_whole():
    from p2d.frame import find_cterm_tag

    seq = "A" * 10 + dna("GSHHHHHHHH") + "TAA" + "C" * 12
    assert find_cterm_tag(Vector.from_sequence(seq, name="t"), 10).residues == "HHHHHHHH"


@pytest.mark.parametrize("up, down", [("NdeI:144", "BamHI:150"), ("EcoRI:168", "XhoI:174")])
def test_touching_sites_are_allowed_and_warned_about(vector, host, up, down):
    a = design(vector, host, up, down)
    codes = {i.code: i for i in a.issues}
    assert a.ok and "sites_adjacent" in codes
    assert codes["sites_adjacent"].severity.value == "warning"
    assert "gel" in codes["sites_adjacent"].message
    real = str(Seq(a.seq[vector.expression_start :][: 3 * len(a.protein)]).translate())
    assert real == a.protein, "the protein from the finished plasmid is the one reported"


def test_sites_with_a_spacer_carry_no_adjacency_warning(vector, host):
    codes = {i.code for i in design(vector, host, "NdeI:144", "XhoI:174").issues}
    assert "sites_adjacent" not in codes


def test_sites_that_share_bases_are_still_refused():
    import random

    rng = random.Random(4)
    rand = lambda n: "".join(rng.choice("ACGT") for _ in range(n))  # noqa: E731
    seq = rand(200).replace("GCATGC", "GACTGC") + "GCATGCAT" + rand(200)
    v = Vector.from_sequence(seq, name="overlap")
    v.expression_start = 10
    sites = {c.enzyme: c for c in list_cut_sites(v)}
    with pytest.raises(ValueError, match="overlap.*share bases"):
        plan_cut(v, sites["SphI"], sites["NsiI"])


def test_a_touching_pair_ranks_below_an_equal_one_with_a_spacer(vector, host):
    from p2d.recommend import rank_pairs

    options, _ = rank_pairs(vector, PAYLOAD, host)
    touching = [o for o in options if o.adjacent]
    assert touching and all(o.up.end == o.down.start for o in touching)
    keys = [o.sort_key for o in options]
    assert keys == sorted(keys)

    import dataclasses

    base = next(o for o in options if not o.adjacent)
    worse = dataclasses.replace(base, adjacent=True)
    assert sorted([worse, base], key=lambda o: o.sort_key) == [base, worse]
    assert base.sort_key < worse.sort_key


def _stop_issue(vector, host):
    a = design(vector, host, "NdeI:144", "XhoI:174", cterm=p2d.CTerm.STOP)
    return next((i for i in a.issues if i.code == "stop_context"), None)


def test_the_base_after_the_stop_is_reported_with_its_source(vector, host):
    issue = _stop_issue(vector, host)
    assert issue and issue.severity.value == "info"
    assert "TAAC" in issue.message and "Poole" in issue.message and "TAAT" in issue.message


def test_a_strong_context_is_not_flagged(vector, tmp_path):
    text = (Path(p2d.__file__).parent / "hosts" / "ecoli_bl21.yaml").read_text()
    assert "strong: [TAAT]" in text
    other = tmp_path / "other.yaml"
    other.write_text(text.replace("strong: [TAAT]", "strong: [TAAC]"))
    assert _stop_issue(vector, p2d.load_host(other)) is None


def test_hosts_without_termination_data_say_nothing(vector):
    assert _stop_issue(vector, p2d.load_host("hsapiens")) is None
    plain = design(vector, p2d.load_host("ecoli_bl21"), "NdeI:144", "XhoI:174")
    assert "stop_context" not in {i.code for i in plain.issues}, "no stop codon in the insert"
