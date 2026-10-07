"""Enzyme analysis and arbitrary enzyme pairs, checked against independent translation."""

import itertools

import pytest
from Bio.Seq import Seq

import p2d
from p2d import enzymes
from p2d.enzymes import EnzymeCandidate, candidate_enzymes, ends_compatible, plan_site
from p2d.vector import StartSource, Vector

PAYLOAD = "SKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLKFICTTGKLPVPWPTL"
STRESS = "ACDEFGHIKLMNPQRSTVWY" * 3


@pytest.fixture(scope="module")
def vector():
    return p2d.load_bundled()


@pytest.fixture(scope="module")
def host():
    return p2d.load_host("ecoli_bl21")


def translate_from(seq: str, start: int) -> str:
    tail = seq[start:]
    return str(Seq(tail[: len(tail) - len(tail) % 3]).translate(to_stop=True))


# ----------------------------------------------------------------- analysis


def test_expression_window_is_the_in_frame_orf(vector):
    w = enzymes.expression_window(vector)
    assert (w.start, w.stop_end) == (87, 201)
    assert vector.seq[w.stop_start : w.stop_end] in {"TAA", "TAG", "TGA"}


def test_candidates_are_the_single_cutters_in_the_window(vector):
    found = {c.name: c for c in candidate_enzymes(vector)}
    assert {"NcoI", "NdeI", "BamHI", "KpnI", "Acc65I", "EcoRI", "XhoI"} == set(found)
    assert [c.start for c in candidate_enzymes(vector)] == sorted(c.start for c in found.values())
    assert found["XhoI"].start == 174 and found["XhoI"].site == "CTCGAG"
    for c in found.values():
        assert len(vector.find_all(c.name)) == 1, "every candidate cuts the vector once"
        assert vector.seq[c.start : c.end] == c.site


def test_isoschizomers_are_merged_under_one_name(vector):
    found = {c.name: c for c in candidate_enzymes(vector)}
    assert "FauNDI" in found["NdeI"].aliases
    assert {"PaeR7I", "SlaI"} <= set(found["XhoI"].aliases)
    assert "FauNDI" not in found, "an alias must not also appear as its own entry"


def test_kpni_and_acc65i_share_a_site_but_not_ends(vector):
    found = {c.name: c for c in candidate_enzymes(vector)}
    assert found["KpnI"].site == found["Acc65I"].site
    assert (found["KpnI"].overhang_kind, found["Acc65I"].overhang_kind) == ("3'", "5'")
    assert not ends_compatible(found["KpnI"], found["Acc65I"])


def test_multi_cutters_and_blunt_and_outside_window_are_excluded(vector):
    names = {c.name for c in candidate_enzymes(vector)}
    assert "SmaI" not in names  # blunt
    assert "BsaI" not in names  # Type IIS, not palindromic
    assert "HindIII" not in names and "SacI" not in names  # do not cut, or cut outside


def test_vector_without_a_curated_start_is_refused_clearly(vector):
    bare = Vector.from_sequence(vector.seq, name="bare")
    with pytest.raises(ValueError, match="expression start"):
        candidate_enzymes(bare)


# ------------------------------------------------------------ end compatibility


def cand(name, kind, ovhg):
    return EnzymeCandidate(name, "X", 0, 1, kind, ovhg)


def test_ends_compatible_matches_known_cases():
    assert ends_compatible(cand("BamHI", "5'", "GATC"), cand("BglII", "5'", "GATC"))
    assert ends_compatible(cand("NcoI", "5'", "CATG"), cand("BspHI", "5'", "CATG"))
    assert not ends_compatible(cand("EcoRI", "5'", "AATT"), cand("XhoI", "5'", "TCGA"))
    assert not ends_compatible(cand("KpnI", "3'", "GTAC"), cand("NcoI", "5'", "GTAC"))


# ------------------------------------------------------------------- plan_site


def test_start_source_follows_where_the_start_codon_is(vector):
    ndei = plan_site(vector, "NdeI", "XhoI")
    assert ndei.start_source is StartSource.VECTOR and ndei.vector_orf_start == 87
    ncoi = plan_site(vector, "NcoI", "XhoI")
    assert ncoi.start_source is StartSource.UPSTREAM_SITE and ncoi.start_offset_in_site == 2


def test_aliases_are_accepted_by_name(vector):
    assert plan_site(vector, "FauNDI", "SlaI").upstream_enzyme == "NdeI"


@pytest.mark.parametrize(
    "up, down, message",
    [
        ("XhoI", "NdeI", "must lie upstream"),
        ("BamHI", "BamHI", "must lie upstream"),
        ("NdeI", "BamHI", "touch or overlap"),
        ("NcoI", "BspHI", "not a single-cutter"),
        ("EcoRI", "SmaI", "not a single-cutter"),
        ("Nonexistent", "XhoI", "not a single-cutter"),
    ],
)
def test_impossible_pairs_say_why(vector, up, down, message):
    with pytest.raises(ValueError, match=message):
        plan_site(vector, up, down)


def test_compatible_ends_are_refused_even_when_in_order(vector, monkeypatch):
    base = candidate_enzymes(vector)
    twin = [
        c if c.name != "XhoI" else EnzymeCandidate(**{**c.__dict__, "overhang": "GATC"})
        for c in base
    ]
    monkeypatch.setattr(enzymes, "candidate_enzymes", lambda v: twin)
    with pytest.raises(ValueError, match="compatible"):
        plan_site(vector, "BamHI", "XhoI")


# ------------------------------------------- every valid pair, checked end to end


def valid_pairs(vector):
    cs = candidate_enzymes(vector)
    for up, down in itertools.permutations(cs, 2):
        try:
            plan_site(vector, up.name, down.name)
        except ValueError:
            continue
        yield up.name, down.name


def test_there_are_several_valid_pairs(vector):
    pairs = set(valid_pairs(vector))
    assert {("NdeI", "XhoI"), ("BamHI", "EcoRI"), ("NcoI", "XhoI"), ("KpnI", "XhoI")} <= pairs
    assert ("XhoI", "NdeI") not in pairs and ("BamHI", "BamHI") not in pairs
    assert ("EcoRI", "XhoI") not in pairs, "abutting sites cannot be double-digested"


@pytest.mark.parametrize("payload", [PAYLOAD, STRESS])
def test_every_valid_pair_designs_and_translates_from_the_real_start(vector, host, payload):
    """The protein p2d reports must be what the plasmid really makes from its own ATG."""
    plan = p2d.ConstructPlan.single_chain(payload)
    for up, down in valid_pairs(vector):
        site = plan_site(vector, up, down)
        a = p2d.assemble(p2d.design_insert(plan, vector, site, host, seed=1))
        assert a.ok, (up, down, a.errors)
        real = translate_from(a.seq, vector.expression_start)
        assert payload in real, (up, down)
        if site.start_source is StartSource.VECTOR:
            assert a.protein == real, (up, down)
        else:
            # the upstream enzyme removes the vector's tag, so p2d's own start is the real one
            assert a.design.orf_start == vector.expression_start
            assert a.protein == real


def test_n_terminal_tag_follows_the_enzyme_choice(vector, host):
    plan = p2d.ConstructPlan.single_chain(PAYLOAD)

    def protein(up):
        site = plan_site(vector, up, "XhoI")
        return p2d.assemble(p2d.design_insert(plan, vector, site, host)).protein

    assert protein("NdeI").startswith("MGSSHHHHHHSSGLVPRGSHM" + PAYLOAD)  # tag kept
    # tag removed, ATG donated by the site; CCATGG reads M-G, so NcoI adds a Gly after the Met
    assert protein("NcoI").startswith("MG" + PAYLOAD)


def test_downstream_enzyme_decides_whether_the_c_terminal_tag_survives(vector, host):
    plan = p2d.ConstructPlan.single_chain(PAYLOAD)
    tagged = p2d.assemble(p2d.design_insert(plan, vector, plan_site(vector, "NdeI", "XhoI"), host))
    assert tagged.protein.endswith("LEHHHHHH")
    # EcoRI lies upstream of XhoI: the XhoI scar and the tag are both kept
    eco = p2d.assemble(p2d.design_insert(plan, vector, plan_site(vector, "NdeI", "EcoRI"), host))
    assert eco.protein.endswith("EFLEHHHHHH")


# -------------------------------------- vectors whose downstream frame is shifted


def shifted_vector(vector, extra: int) -> Vector:
    """Pad the stuffer between BamHI and EcoRI by ``extra`` bases: sites downstream of it
    now start out of frame with the start codon, as in the pET-28 a/b/c series."""
    seq = vector.seq[:162] + "A" * extra + vector.seq[162:]
    v = Vector.from_sequence(seq, name=f"shift{extra}")
    v.expression_start = 87
    return v


@pytest.mark.parametrize("extra", [1, 2])
def test_downstream_frame_is_preserved_when_the_site_is_out_of_frame(vector, host, extra):
    """The vector's own C-terminal reading must continue exactly as before the splice."""
    v = shifted_vector(vector, extra)
    plan = p2d.ConstructPlan.single_chain(PAYLOAD)
    checked = 0
    for up, down in valid_pairs(v):
        d = p2d.design_insert(plan, v, plan_site(v, up, down), host, seed=2)
        a = p2d.assemble(d)
        assert a.ok, (up, down, a.errors)
        boundary = d.down.end + (87 - d.down.end) % 3
        tail = translate_from(v.seq, boundary)
        assert tail and a.protein.endswith(tail), (extra, up, down)
        assert (d.coding_start - d.orf_start) % 3 == 0
        checked += 1
    assert checked, "the shifted vector must still offer pairs to test"


def test_out_of_frame_downstream_site_gets_a_three_prime_pad(vector, host):
    v = shifted_vector(vector, 1)
    plan = p2d.ConstructPlan.single_chain(PAYLOAD)
    pads = {
        len(p2d.design_insert(plan, v, plan_site(v, "NdeI", dn), host).pad_down)
        for dn in ("EcoRI", "XhoI")
        if dn in {c.name for c in candidate_enzymes(v)}
    }
    assert pads and pads != {0}
