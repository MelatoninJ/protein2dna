"""The DFA + DP codon optimiser, tested against an independent regex oracle."""

import math

import pytest
from Bio.Seq import Seq
from hypothesis import given, settings
from hypothesis import strategies as st

import p2d
from oracle import region_hits
from p2d import frame
from p2d.frame import reverse_translate
from p2d.optimize import CodonStrategy, OptimizationRequest, optimize_coding

ALL_20 = "ACDEFGHIKLMNPQRSTVWY"
PAYLOAD = "KVFLDWINEAYQRGTRVLAEMAKRG"
SITES = ["GAATTC", "CTCGAG", "GGTCTC", "CATATG", "GRATCY"]


@pytest.fixture(scope="module")
def host():
    return p2d.load_host("ecoli_bl21")


def request(host, aa, **kw):
    kw.setdefault("prefix", "")
    kw.setdefault("suffix", "")
    kw.setdefault("forbidden", [])
    return OptimizationRequest(aa=aa, host=host, **kw)


def cai(host, aa, dna):
    logs = []
    for i, residue in enumerate(aa):
        usage = host.codons_for(residue)
        logs.append(math.log(usage[dna[3 * i : 3 * i + 3]] / max(usage.values())))
    return math.exp(sum(logs) / len(logs))


def translate(dna):
    return str(Seq(dna).translate())


STRATEGIES = [CodonStrategy.MAX_CAI, CodonStrategy.WEIGHTED_SAMPLE]


# ---------------------------------------------------------- back-translation


@pytest.mark.parametrize("strategy", STRATEGIES)
@pytest.mark.parametrize(
    "aa", [PAYLOAD, "M", "W", "MW" * 30, "L" * 40, "SRLGP" * 12, "A" * 7 + "K" * 7]
)
def test_back_translation(host, aa, strategy):
    r = optimize_coding(request(host, aa, strategy=strategy, seed=3, forbidden=SITES))
    assert r.ok, r.issues
    assert translate(r.dna) == aa


@pytest.mark.parametrize("strategy", STRATEGIES)
def test_all_twenty_residues_round_trip(host, strategy):
    for aa in (ALL_20, ALL_20[::-1], ALL_20 * 3):
        r = optimize_coding(request(host, aa, strategy=strategy, seed=1, forbidden=SITES))
        assert r.ok, r.issues
        assert translate(r.dna) == aa


def test_every_single_residue_is_reachable(host):
    for residue in ALL_20:
        r = optimize_coding(request(host, residue, strategy=CodonStrategy.MAX_CAI))
        assert translate(r.dna) == residue


def test_empty_payload(host):
    r = optimize_coding(request(host, ""))
    assert r.ok and r.dna == "" and r.metrics["n_codons"] == 0


def test_unknown_residue_is_a_clear_error(host):
    with pytest.raises(ValueError, match="no codons for residue 'Z'"):
        optimize_coding(request(host, "AZA"))


# ------------------------------------------------------------ forbidden sites


@pytest.mark.parametrize("strategy", STRATEGIES)
def test_no_forbidden_site_in_assembled_context(host, strategy):
    prefix, suffix = "TTTCATATG", "CTCGAGAAA"
    r = optimize_coding(
        request(host, PAYLOAD * 3, prefix=prefix, suffix=suffix, forbidden=SITES, strategy=strategy)
    )
    assert r.ok
    assert region_hits(prefix, r.dna, suffix, SITES) == set()


def test_site_present_only_in_fixed_context_is_not_a_failure(host):
    """The enzyme's own site in prefix/suffix is intentional; only overlaps count."""
    r = optimize_coding(request(host, "KK", prefix="GAATTC", suffix="CTCGAG", forbidden=SITES))
    assert r.ok


def test_reverse_complement_strand_is_avoided(host):
    # GGTTTC is not palindromic; its reverse complement GAAACC is E+T at their best codons
    aa = "ET"
    naive = reverse_translate(aa, host)
    assert region_hits("", naive, "", ["GGTTTC"]), "premise: naive path creates the RC site"
    r = optimize_coding(request(host, aa, forbidden=["GGTTTC"], strategy=CodonStrategy.MAX_CAI))
    assert r.ok and translate(r.dna) == aa
    assert region_hits("", r.dna, "", ["GGTTTC"]) == set()


def test_ambiguity_codes_are_expanded(host):
    aa = "EI"  # best codons GAA ATT -> GAAATT, which GAAATN matches
    assert region_hits("", reverse_translate(aa, host), "", ["GAAATN"])
    r = optimize_coding(request(host, aa, forbidden=["GAAATN"], strategy=CodonStrategy.MAX_CAI))
    assert r.ok and region_hits("", r.dna, "", ["GAAATN"]) == set()


def test_pattern_spanning_prefix_and_coding_is_avoided(host):
    """Optimising the region in isolation would miss this one."""
    prefix, aa, forbidden = "AAGA", "IL", ["GAATTC"]  # GA + ATT + CTG at best codons
    assert region_hits(prefix, reverse_translate(aa, host), "", forbidden)

    r = optimize_coding(
        request(host, aa, prefix=prefix, forbidden=forbidden, strategy=CodonStrategy.MAX_CAI)
    )
    assert r.ok and translate(r.dna) == aa
    assert region_hits(prefix, r.dna, "", forbidden) == set()
    assert any(i.code == "constraint_swaps" for i in r.issues)


def test_pattern_spanning_coding_and_suffix_is_avoided(host):
    aa, suffix, forbidden = "KI", "C", ["AAATTC"]  # AAA + ATT at best codons, then C
    assert region_hits("", reverse_translate(aa, host), suffix, forbidden)

    r = optimize_coding(
        request(host, aa, suffix=suffix, forbidden=forbidden, strategy=CodonStrategy.MAX_CAI)
    )
    assert r.ok
    assert region_hits("", r.dna, suffix, forbidden) == set()


def test_prefix_length_beyond_the_longest_pattern_changes_nothing(host):
    """Only the last ``max_len`` bases of context can matter; MAX_CAI has no seed to vary."""
    kw = dict(forbidden=["GAATTC"], strategy=CodonStrategy.MAX_CAI)
    short = optimize_coding(request(host, "IL", prefix="AAGA", **kw))
    long = optimize_coding(request(host, "IL", prefix="TTTTGGGGCCCC" * 9 + "AAGA", **kw))
    assert short.dna == long.dna
    assert region_hits("AAGA", short.dna, "", ["GAATTC"]) == set()


# ----------------------------------------------------------------- no path


def test_unavoidable_site_is_reported_not_raised(host):
    # methionine has exactly one codon, so ATG cannot be avoided
    r = optimize_coding(request(host, "AMA", forbidden=["ATG"]))
    assert not r.ok
    (issue,) = [i for i in r.issues if i.code == "unavoidable_site"]
    assert issue.severity is p2d.Severity.ERROR
    assert "ATG" in issue.message
    assert "residue 2" in issue.message and "(M)" in issue.message
    assert "different enzyme pair" in issue.message


def test_unavoidable_site_across_the_fixed_junction(host):
    # W is TGG only; with a fixed 3' context of AT that always completes TGGAT
    r = optimize_coding(request(host, "W", suffix="AT", forbidden=["TGGAT"]))
    assert not r.ok
    (issue,) = [i for i in r.issues if i.code == "unavoidable_site"]
    assert "3' junction" in issue.message and "TGGAT" in issue.message


def test_unavoidable_site_across_the_5prime_junction(host):
    r = optimize_coding(request(host, "WA", prefix="GGA", forbidden=["GGATGG"]))
    assert not r.ok
    (issue,) = [i for i in r.issues if i.code == "unavoidable_site"]
    assert "residue 1" in issue.message


def test_failure_still_returns_translatable_fallback(host):
    r = optimize_coding(request(host, "AMA", forbidden=["ATG"]))
    assert translate(r.dna) == "AMA"
    assert r.metrics["feasible"] is False


def test_oversized_ambiguous_pattern_is_a_clear_error(host):
    with pytest.raises(ValueError, match="expands to more than"):
        optimize_coding(request(host, "AA", forbidden=["N" * 8]))


# --------------------------------------------------------------- strategies


def test_same_seed_is_byte_identical(host):
    kw = dict(forbidden=SITES, prefix="CATATG", suffix="CTCGAG")
    a = optimize_coding(request(host, PAYLOAD * 2, seed=42, **kw))
    b = optimize_coding(request(host, PAYLOAD * 2, seed=42, **kw))
    assert a.dna == b.dna and a.metrics == b.metrics


def test_different_seed_gives_different_dna(host):
    aa = PAYLOAD * 3
    dnas = {optimize_coding(request(host, aa, seed=s)).dna for s in range(5)}
    assert len(dnas) == 5
    assert {translate(d) for d in dnas} == {aa}


def test_unseeded_call_is_still_reproducible(host):
    a = optimize_coding(request(host, PAYLOAD * 2, forbidden=SITES))
    b = optimize_coding(request(host, PAYLOAD * 2, forbidden=SITES))
    assert a.dna == b.dna


def test_max_cai_is_never_worse_than_the_naive_baseline(host):
    """Without sites in the way, MAX_CAI is exactly the naive best-codon choice."""
    naive = reverse_translate(PAYLOAD, host)
    r = optimize_coding(request(host, PAYLOAD, strategy=CodonStrategy.MAX_CAI))
    assert r.dna == naive
    assert r.metrics["cai"] == pytest.approx(cai(host, PAYLOAD, naive)) == pytest.approx(1.0)


def test_max_cai_beats_sampling_under_the_same_constraints(host):
    aa = PAYLOAD * 4
    best = optimize_coding(request(host, aa, forbidden=SITES, strategy=CodonStrategy.MAX_CAI))
    for seed in range(5):
        sampled = optimize_coding(
            request(host, aa, forbidden=SITES, strategy=CodonStrategy.WEIGHTED_SAMPLE, seed=seed)
        )
        assert best.metrics["cai"] >= sampled.metrics["cai"]


def test_constrained_max_cai_is_optimal(host):
    """Brute force over every synonymous choice for a short peptide."""
    from itertools import product

    aa, forbidden = "ILE", ["GAATTC", "AATTCT", "TTCTCG"]
    best = optimize_coding(request(host, aa, forbidden=forbidden, strategy=CodonStrategy.MAX_CAI))
    scores = []
    for combo in product(*(host.codons_for(r) for r in aa)):
        dna = "".join(combo)
        if not region_hits("", dna, "", forbidden):
            scores.append(cai(host, aa, dna))
    assert best.metrics["cai"] == pytest.approx(max(scores))


def test_sampling_respects_the_usage_floor(host):
    aa = ALL_20 * 6
    r = optimize_coding(request(host, aa, seed=9, usage_floor=0.10))
    for i, residue in enumerate(aa):
        assert host.codons_for(residue)[r.dna[3 * i : 3 * i + 3]] >= 0.10


def test_floor_is_relaxed_before_declaring_a_site_unavoidable(host):
    # leucine's only codon under 0.10 is CTA; forbidding the rest leaves just CTA
    forbidden = ["CTG", "TTA", "TTG", "CTC", "CTT"]
    r = optimize_coding(request(host, "L", forbidden=forbidden, seed=1))
    assert r.ok and r.dna == "CTA"
    assert any(i.code == "usage_floor_relaxed" for i in r.issues)


def test_harmonized_is_explicitly_unimplemented(host):
    with pytest.raises(NotImplementedError):
        optimize_coding(request(host, "AAA", strategy=CodonStrategy.HARMONIZED))


def test_metrics_are_sane(host):
    r = optimize_coding(request(host, PAYLOAD, strategy=CodonStrategy.MAX_CAI))
    assert {"cai", "gc", "n_rare"} <= set(r.metrics)
    assert 0 < r.metrics["cai"] <= 1
    assert r.metrics["gc"] == pytest.approx((r.dna.count("G") + r.dna.count("C")) / len(r.dna))
    assert r.metrics["n_rare"] == sum(
        r.dna[i : i + 3] in host.rare_codons for i in range(0, len(r.dna), 3)
    )


# ----------------------------------------------------------- property-based


@settings(max_examples=150, deadline=None)
@given(
    aa=st.text(alphabet=ALL_20, min_size=1, max_size=80),
    strategy=st.sampled_from(STRATEGIES),
    seed=st.integers(0, 2**32),
    prefix=st.text(alphabet="ACGT", max_size=14),
    suffix=st.text(alphabet="ACGT", max_size=14),
)
def test_property_back_translation_and_pattern_freedom(host, aa, strategy, seed, prefix, suffix):
    r = optimize_coding(
        request(
            host,
            aa,
            strategy=strategy,
            seed=seed,
            prefix=prefix,
            suffix=suffix,
            forbidden=SITES,
        )
    )
    assert translate(r.dna) == aa
    if r.ok:
        assert region_hits(prefix, r.dna, suffix, SITES) == set()
    else:
        # a refusal must be a genuine one: it names a site and a position
        assert any(i.code == "unavoidable_site" for i in r.issues)


# ---------------------------------------------------- through design_insert


def test_design_insert_uses_the_optimiser_and_exposes_metrics():
    a = p2d.run(PAYLOAD)
    assert a.ok
    assert a.design.metrics["strategy"] == "weighted_sample"
    assert 0 < a.design.metrics["cai"] <= 1


def test_design_insert_all_residues_and_strategies():
    vector = p2d.load_bundled()
    host = p2d.load_host("ecoli_bl21")
    plan = p2d.ConstructPlan.single_chain(ALL_20 * 2)
    for name in vector.sites:
        for strategy in STRATEGIES:
            d = p2d.design_insert(plan, vector, vector.sites[name], host, strategy=strategy, seed=2)
            assert p2d.assemble(d).ok, (name, strategy)


def test_design_insert_seed_changes_dna_not_protein():
    vector = p2d.load_bundled()
    host = p2d.load_host("ecoli_bl21")
    plan = p2d.ConstructPlan.single_chain(PAYLOAD * 2)
    site = vector.sites["NdeI-XhoI"]
    a = p2d.assemble(p2d.design_insert(plan, vector, site, host, seed=1))
    b = p2d.assemble(p2d.design_insert(plan, vector, site, host, seed=2))
    assert a.insert_dna != b.insert_dna
    assert a.protein == b.protein


def test_design_insert_reports_unavoidable_site_as_error(monkeypatch):
    """An enzyme site that the protein itself forces must fail loudly through the API."""
    vector = p2d.load_bundled()
    host = p2d.load_host("ecoli_bl21")
    site = vector.sites["NdeI-XhoI"]
    real = frame.enzyme_site
    monkeypatch.setattr(
        frame, "enzyme_site", lambda n: "ATG" if n == site.downstream_enzyme else real(n)
    )
    plan = p2d.ConstructPlan.single_chain("AAMAA")  # Met is ATG only
    d = p2d.design_insert(plan, vector, site, host)
    assert not d.ok
    assert any(i.code == "unavoidable_site" for i in d.errors)
