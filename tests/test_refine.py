"""GC and repeat refinement: measured against brute force, and never at the cost of a site."""

import random

import pytest
from hypothesis import given, settings
from hypothesis import strategies as st

import p2d
from oracle import region_hits
from p2d.optimize import CodonStrategy, Limits, OptimizationRequest, optimize_coding
from p2d.optimize.refine import measure, refine

SITES = ["GAATTC", "CTCGAG", "CCATGG"]


@pytest.fixture(scope="module")
def host():
    return p2d.load_host("ecoli_bl21")


def rc(s):
    return s[::-1].translate(str.maketrans("ACGT", "TGCA"))


def translate(dna):
    from Bio.Seq import Seq

    return str(Seq(dna).translate())


def brute_window_excess(dna, lim):
    w, total = lim.window_size, 0.0
    for s in range(len(dna) - w + 1):
        share = sum(c in "GC" for c in dna[s : s + w]) / w
        total += 100 * max(0.0, lim.window_min - share, share - lim.window_max)
    return total


def brute_has_repeat(dna, limit, inverted):
    """Is there a repeat longer than ``limit``: two non-overlapping copies, direct or inverted."""
    k = limit + 1
    for i in range(len(dna) - k + 1):
        kmer = dna[i : i + k]
        target = rc(kmer) if inverted else kmer
        j = dna.find(target)
        while j != -1:
            if abs(j - i) >= k:
                return True
            j = dna.find(target, j + 1)
    return False


def request(host, aa, limits, **kw):
    return OptimizationRequest(
        aa=aa,
        host=host,
        prefix="CCATGG",
        suffix="CTCGAG",
        forbidden=SITES,
        soft_forbidden=host.soft_motifs(),
        strategy=CodonStrategy.MAX_CAI,
        limits=limits,
        **kw,
    )


def protein(seed, letters, n):
    rng = random.Random(seed)
    return "".join(rng.choice(letters) for _ in range(n))


# ---------------------------------------------------------------------- the limits


def test_the_host_profile_supplies_the_defaults(host):
    lim = host.limits()
    assert (lim.gc_min, lim.gc_max) == (0.32, 0.6)
    assert (lim.window_size, lim.window_min, lim.window_max) == (50, 0.28, 0.7)
    assert (lim.direct_max, lim.inverted_max) == (22, 200)


def test_single_limits_can_be_replaced(host):
    lim = host.limits(gc_max=0.55, direct_max=30)
    assert lim.gc_max == 0.55 and lim.direct_max == 30 and lim.gc_min == 0.32


@pytest.mark.parametrize(
    "bad",
    [{"gc_min": 0.7}, {"window_min": 0.9}, {"window_size": 5}, {"direct_max": 3}, {"gc_max": 1.5}],
)
def test_impossible_limits_are_refused(host, bad):
    with pytest.raises(ValueError):
        host.limits(**bad)


# --------------------------------------------------------------------- measuring


@given(st.text("ACGT", min_size=0, max_size=300), st.integers(10, 60))
@settings(max_examples=80, deadline=None)
def test_window_gc_matches_brute_force(dna, size):
    lim = Limits(window_size=size)
    m = measure(dna, lim)
    assert m.window_excess == pytest.approx(brute_window_excess(dna, lim))
    if dna:
        assert m.gc == pytest.approx(sum(c in "GC" for c in dna) / len(dna))


@given(st.text("ACGT", min_size=0, max_size=160))
@settings(max_examples=120, deadline=None)
def test_repeat_detection_matches_brute_force(dna):
    lim = Limits(direct_max=8, inverted_max=8)
    m = measure(dna, lim)
    expected = brute_has_repeat(dna, 8, False) or brute_has_repeat(dna, 8, True)
    assert (m.repeats > 0) == expected


def test_overlapping_copies_do_not_count_as_a_repeat():
    lim = Limits(direct_max=20)
    assert measure("A" * 40, lim).repeats == 0, "too short to hold two separate copies"
    assert measure("A" * 100, lim).repeats > 0, "long enough to hold two"


def test_copies_that_only_overlap_are_not_a_repeat():
    flank = protein(11, "ACGT", 30)
    dna = flank + "A" * 24 + flank[::-1].replace("A", "C")
    lim = Limits(direct_max=20, inverted_max=20)
    assert not brute_has_repeat(dna, 20, False)
    assert measure(dna, lim).repeats == 0


def test_an_inverted_repeat_is_found_and_counted_separately():
    unit = protein(12, "ACGT", 26)
    dna = unit + "TTTTTTTTTTTT" + rc(unit)
    assert brute_has_repeat(dna, 20, True) and not brute_has_repeat(dna, 20, False)
    assert measure(dna, Limits(direct_max=20, inverted_max=20)).repeats > 0
    assert measure(dna, Limits(direct_max=20, inverted_max=40)).repeats == 0


# --------------------------------------------------------------------- refinement


@pytest.mark.parametrize("letters", ["KNIFYMKNIL", "AGPRAGPRWH"], ids=["AT-rich", "GC-rich"])
def test_extreme_proteins_are_brought_inside_the_limits(host, letters):
    aa, lim = protein(1, letters, 300), host.limits()
    plain = optimize_coding(request(host, aa, None))
    refined = optimize_coding(request(host, aa, lim))
    assert measure(plain.dna, lim).penalty > 0
    after = measure(refined.dna, lim)
    assert after.gc_excess == 0 and after.window_excess == 0
    assert translate(refined.dna) == aa
    assert refined.metrics["refined_swaps"] > 0
    assert any(i.code == "refined" for i in refined.issues)


def test_nothing_changes_when_the_sequence_is_already_inside(host):
    aa = protein(2, "ACDEFGHIKLMNPQRSTVWY", 200)
    wide = host.limits(gc_min=0.0, gc_max=1.0, window_min=0.0, window_max=1.0, direct_max=200)
    assert (
        optimize_coding(request(host, aa, wide)).dna == optimize_coding(request(host, aa, None)).dna
    )


def test_repeats_made_of_identical_codons_are_broken(host):
    unit = protein(3, "ACDEFGHIKLMNPQRSTVWY", 30)
    aa = unit + protein(4, "ACDEFGHIKLMNPQRSTVWY", 20) + unit
    lim = host.limits()
    assert measure(optimize_coding(request(host, aa, None)).dna, lim).repeats > 0
    refined = optimize_coding(request(host, aa, lim))
    assert measure(refined.dna, lim).repeats == 0 and translate(refined.dna) == aa


def test_refinement_never_adds_a_site(host):
    for letters in ("KNIFYMKNIL", "AGPRAGPRWH", "ACDEFGHIKLMNPQRSTVWY"):
        aa = protein(5, letters, 250)
        for lim in (None, host.limits()):
            r = optimize_coding(request(host, aa, lim))
            found = region_hits("CCATGG", r.dna, "CTCGAG", SITES)
            assert not found, (letters, found)


def test_refinement_does_not_add_host_motifs_it_did_not_have(host):
    motifs = list(host.soft_motifs())
    for seed in range(6):
        aa = protein(seed, "KNIFYMKNILAGPR", 250)
        plain = optimize_coding(request(host, aa, None))
        refined = optimize_coding(request(host, aa, host.limits()))
        before = region_hits("CCATGG", plain.dna, "CTCGAG", motifs)
        after = region_hits("CCATGG", refined.dna, "CTCGAG", motifs)
        assert len(after) <= len(before)


def test_refinement_is_deterministic(host):
    aa = protein(6, "KNIFYMKNIL", 200)
    a = optimize_coding(request(host, aa, host.limits())).dna
    assert a == optimize_coding(request(host, aa, host.limits())).dna


def test_what_the_protein_forces_is_reported_not_hidden(host):
    all_w = optimize_coding(request(host, "W" * 60, host.limits()))  # one codon, 67% GC
    codes = {i.code for i in all_w.issues}
    assert {"gc_out_of_range"} <= codes and translate(all_w.dna) == "W" * 60
    strict = host.limits(gc_min=0.4, window_min=0.4)
    all_k = optimize_coding(request(host, "K" * 60, strict))
    codes = {i.code for i in all_k.issues}
    assert {"gc_out_of_range", "gc_window_out_of_range"} <= codes
    assert measure(all_k.dna, strict).gc == pytest.approx(1 / 3), "the most GC it can reach"


def test_repeats_the_protein_forces_are_reported(host):
    aa = "M" * 40  # a single codon: nothing to swap
    r = optimize_coding(request(host, aa, host.limits()))
    assert "repeat_kept" in {i.code for i in r.issues}


def test_stricter_limits_are_honoured(host):
    aa = protein(7, "ACDEFGHIKLMNPQRSTVWY", 300)
    tight = host.limits(gc_min=0.48, gc_max=0.52)
    r = optimize_coding(request(host, aa, tight))
    assert 0.48 <= measure(r.dna, tight).gc <= 0.52


def test_refine_swaps_only_synonymous_codons(host):
    aa = protein(8, "KNIFYMKNIL", 120)
    codons = [host.best_codon(c) for c in aa]
    usage = {c: host.codons_for(c) for c in set(aa)}
    out = refine(aa, codons, usage, "CCATGG", "CTCGAG", SITES, host.limits())
    assert all(new in usage[a] for new, a in zip(out.codons, aa, strict=True))
    assert out.after.penalty < out.before.penalty


@given(st.text("ACDEFGHIKLMNPQRSTVWY", min_size=1, max_size=120), st.integers(0, 3))
@settings(max_examples=40, deadline=None)
def test_any_protein_keeps_its_translation_and_stays_site_free(aa, tightness):
    host = p2d.load_host("ecoli_bl21")
    lim = host.limits(gc_min=0.40 + 0.03 * tightness, gc_max=0.60 - 0.03 * tightness)
    r = optimize_coding(request(host, aa, lim))
    assert translate(r.dna) == aa
    if r.ok:
        assert not region_hits("CCATGG", r.dna, "CTCGAG", SITES)


# ------------------------------------------------------------ through the design


def test_a_full_design_uses_the_host_limits_by_default(host):
    from p2d.enzymes import list_cut_sites, plan_cut

    aa = protein(9, "KNIFYMKNIL", 150)
    plan = p2d.ConstructPlan.single_chain(aa)
    v = p2d.load_bundled()
    sites = list_cut_sites(v, None)
    site = plan_cut(
        v,
        next(s for s in sites if s.enzyme == "NcoI"),
        next(s for s in sites if s.enzyme == "XhoI"),
    )

    def build(**kw):
        return p2d.design_insert(
            plan, v, site, host, strategy=CodonStrategy.MAX_CAI, cterm=p2d.CTerm.STOP, **kw
        )

    assert measure(build().coding, host.limits()).window_excess == 0
    loose = host.limits(gc_min=0, gc_max=1, window_min=0, window_max=1, direct_max=200)
    assert measure(build(limits=loose).coding, host.limits()).window_excess > 0


def test_fixing_repeats_never_costs_gc(host):
    """A protein made of one repeated unit: repeats and GC pull in different directions."""
    aa = "KNIFYMKNIL" * 6 + "KNIF"
    lim = host.limits()
    r = optimize_coding(request(host, aa, lim))
    best = "".join(max(host.codons_for(a), key=lambda c: sum(x in "GC" for x in c)) for a in aa)
    top = measure(best, lim).gc
    got = measure(r.dna, lim)
    assert got.gc >= min(top, lim.gc_min) - 0.02, "GC was sacrificed for repeats"
    gc_only = optimize_coding(request(host, aa, host.limits(direct_max=200, inverted_max=200)))
    assert got.gc_excess <= measure(gc_only.dna, lim).gc_excess + 1e-9
