"""The DNA to order: clamps, and the in-silico digest and ligation that proves it clones."""

import re

import pytest
from Bio.Restriction import RestrictionBatch
from Bio.Seq import Seq
from Bio.SeqUtils import gc_fraction
from hypothesis import given, settings
from hypothesis import strategies as st

import p2d
from oracle import hits
from p2d import make_order
from p2d.enzymes import find_cut_site, list_cut_sites, plan_cut
from p2d.order import MAX_CLAMP, simulate_cloning
from p2d.vector import Vector

PAYLOAD = "KVFLDWINEAYQRGTRVLAEMAKRGDEFVKRLIAEGHDPFEVLKELGYSE"


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


PAIRS = [
    ("NdeI:144", "XhoI:174"),
    ("BamHI:150", "XhoI:174"),
    ("NdeI:144", "EcoRI:168"),
    ("NcoI:85", "XhoI:174"),
]


def test_the_order_is_clamp_insert_clamp(vector, host):
    a = design(vector, host, "NdeI:144", "XhoI:174")
    o = make_order(a, 6, 4)
    assert o.dna == o.clamp_5 + a.insert_dna + o.clamp_3
    assert (len(o.clamp_5), len(o.clamp_3)) == (6, 4) and o.length == len(a.insert_dna) + 10
    assert o.dna[o.insert_start : o.insert_end] == a.insert_dna
    assert o.dna[slice(*o.coding_span)] == a.design.coding, "coding span follows the clamp"


def test_no_clamp_means_the_insert_itself(vector, host):
    a = design(vector, host, "NdeI:144", "XhoI:174")
    o = make_order(a, 0, 0)
    assert o.dna == a.insert_dna and o.ligation.ok


def test_the_default_is_six_bases_each_side(vector, host):
    a = design(vector, host, "NdeI:144", "XhoI:174")
    o = make_order(a)
    assert (len(o.clamp_5), len(o.clamp_3)) == (6, 6)


def test_same_design_same_clamps_and_different_designs_differ(vector, host):
    a = design(vector, host, "NdeI:144", "XhoI:174")
    assert make_order(a).dna == make_order(a).dna
    b = design(vector, host, "BamHI:150", "XhoI:174")
    assert (make_order(a).clamp_5, make_order(a).clamp_3) != (
        make_order(b).clamp_5,
        make_order(b).clamp_3,
    )


@pytest.mark.parametrize("n", [2, 3, 5, 6, 8, 12, 20, 30])
def test_clamps_are_not_palindromic_or_repetitive(vector, host, n):
    a = design(vector, host, "NdeI:144", "XhoI:174")
    o = make_order(a, n, n)
    for clamp in (o.clamp_5, o.clamp_3):
        assert str(Seq(clamp).reverse_complement()) != clamp
        assert not re.search(r"(.)\1\1", clamp)
        if n >= 6:
            assert 0.34 <= gc_fraction(clamp) <= 0.67


@pytest.mark.parametrize("n", [0, 1, 6, 12, 30])
@pytest.mark.parametrize("up, down", PAIRS)
def test_clamps_never_create_a_cloning_site_or_host_motif(vector, host, up, down, n):
    a = design(vector, host, up, down)
    o = make_order(a, n, n)
    patterns = [
        a.design.vector.seq[a.design.up.start : a.design.up.end],
        a.design.vector.seq[a.design.down.start : a.design.down.end],
        *host.soft_motifs(),
    ]
    for start, end in hits(o.dna, patterns):
        assert (
            start >= o.insert_start - 0
            and end <= o.insert_end
            or (start >= o.insert_start and end <= o.insert_end)
        ), f"a site at {start}..{end} touches a clamp"


@pytest.mark.parametrize("n", [0, 1, 3, 6, 12, 30])
@pytest.mark.parametrize("up, down", PAIRS)
def test_the_order_digests_and_ligates_into_the_planned_plasmid(vector, host, up, down, n):
    a = design(vector, host, up, down)
    o = make_order(a, n, n)
    assert o.ok and o.ligation.ok, o.ligation.message
    assert "exactly" in o.ligation.message


def test_an_independent_digest_agrees(vector, host):
    """Biopython's catalyse, run by the test, gives the same middle fragment."""
    a = design(vector, host, "BamHI:150", "XhoI:174")
    o = make_order(a, 6, 6)
    up, down = (RestrictionBatch([n]).get(n) for n in ("BamHI", "XhoI"))
    after_up = up.catalyse(Seq(o.dna))  # [clamp piece, rest]
    middle = down.catalyse(after_up[1])[0]
    v = a.design.vector.seq
    cut_a, cut_b = 150 + up.fst5, 174 + down.fst5
    assert v[:cut_a] + str(middle) + v[cut_b:] == a.seq


def test_enzymes_that_cut_the_fragment_elsewhere_are_caught(vector, host):
    a = design(vector, host, "NdeI:144", "XhoI:174")
    o = make_order(a, 6, 6)
    bad = "CTCGAG" + o.dna  # a second XhoI in front, outside the intended sites
    check = simulate_cloning(a, bad)
    assert not check.ok and "XhoI cuts the ordered fragment 2 times" in check.message


def test_a_dropped_base_is_caught(vector, host):
    a = design(vector, host, "NdeI:144", "XhoI:174")
    o = make_order(a, 6, 6)
    bad = o.dna[:60] + o.dna[61:]  # one base missing from the coding region
    check = simulate_cloning(a, bad)
    assert not check.ok and "different plasmid" in check.message


def test_a_fragment_with_no_site_is_caught(vector, host):
    a = design(vector, host, "NdeI:144", "XhoI:174")
    check = simulate_cloning(a, "ACGT" * 40)
    assert not check.ok and "cuts the ordered fragment 0 times" in check.message


def test_a_failed_design_cannot_produce_an_orderable_sequence(vector, host, monkeypatch):
    a = design(vector, host, "NdeI:144", "XhoI:174")
    monkeypatch.setattr(type(a), "ok", property(lambda self: False))
    o = make_order(a)
    assert not o.ok and any(i.code == "design_has_errors" for i in o.issues)


def test_order_metrics(vector, host):
    a = design(vector, host, "NdeI:144", "XhoI:174")
    o = make_order(a, 6, 6)
    assert o.gc == pytest.approx(gc_fraction(o.dna))
    assert o.longest_homopolymer == max(len(m.group()) for m in re.finditer(r"(.)\1*", o.dna))
    assert 1 <= o.longest_homopolymer <= 8


@pytest.mark.parametrize("bad", [-1, MAX_CLAMP + 1, 2.5, "6", None, True])
def test_bad_clamp_lengths_are_refused(vector, host, bad):
    a = design(vector, host, "NdeI:144", "XhoI:174")
    with pytest.raises(ValueError, match="whole number"):
        make_order(a, bad, 6)
    with pytest.raises(ValueError, match="whole number"):
        make_order(a, 6, bad)


def test_pads_and_clamps_together_still_clone(vector, host):
    """A shifted downstream frame needs a 3' pad; the clamp must sit outside it."""
    shifted = Vector.from_sequence(vector.seq[:162] + "A" + vector.seq[162:], name="shift1")
    shifted.expression_start = 87
    xho = next(c.id for c in list_cut_sites(shifted) if c.enzyme == "XhoI")
    a = design(shifted, host, "BamHI:150", xho)  # keeping the vector's frame needs a 3' pad
    assert a.design.pad_down != "", "premise: this design needs a pad"
    for n in (0, 6, 12):
        o = make_order(a, n, n)
        assert o.ok, o.ligation.message
        assert o.dna[o.insert_start : o.insert_end] == a.insert_dna


@settings(max_examples=60, deadline=None)
@given(
    aa=st.text(alphabet="ACDEFGHIKLMNPQRSTVWY", min_size=1, max_size=40),
    n5=st.integers(0, 12),
    n3=st.integers(0, 12),
    pair=st.sampled_from(PAIRS),
    cterm=st.sampled_from([p2d.CTerm.STOP, p2d.CTerm.TAG]),
)
def test_property_every_design_orders_and_clones(aa, n5, n3, pair, cterm):
    v, h = p2d.load_bundled(), p2d.load_host("ecoli_bl21")
    a = design(v, h, *pair, payload=aa, cterm=cterm)
    o = make_order(a, n5, n3)
    assert o.ok, (pair, o.ligation.message)
    assert o.dna[o.insert_start : o.insert_end] == a.insert_dna
    up = RestrictionBatch([a.design.site.upstream_enzyme]).get(a.design.site.upstream_enzyme)
    assert len(up.search(Seq(o.dna), linear=True)) == 1, "the 5' enzyme cuts the fragment once"


def test_real_pet28a_every_valid_pair_with_clamps(host):
    import os
    from pathlib import Path

    from p2d.analysis import analyze_vector, confirm
    from p2d.vectorio import read_vector_file

    samples = os.environ.get("P2D_VECTOR_SAMPLES")
    path = Path(samples) / "pET-28a_plus.dna" if samples else None
    if not path or not path.is_file():
        pytest.skip("pET-28a_plus.dna not available (set P2D_VECTOR_SAMPLES)")
    v = confirm(analyze_vector(read_vector_file(path)))
    sites = list_cut_sites(v, (v.expression_start, 5260))
    tried = 0
    for up in sites:
        for down in sites:
            try:
                plan_cut(v, up, down)
            except ValueError:
                continue
            for cterm in (p2d.CTerm.STOP, p2d.CTerm.TAG):
                try:
                    a = design(v, host, up.id, down.id, cterm=cterm)
                except ValueError:
                    continue  # e.g. TAG with no tag after this 3' site
                o = make_order(a, 6, 6)
                assert o.ok, (up.label, down.label, cterm, o.ligation.message)
                tried += 1
    assert tried >= 30
