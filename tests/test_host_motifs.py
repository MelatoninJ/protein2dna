"""Host preferences (internal Shine-Dalgarno, Chi, homopolymer runs) as soft constraints.

They are avoided whenever a path exists and never make a design fail; whatever stays is
named in a warning.  The check that nothing is violated *silently* runs against an
independent regex oracle.
"""

import re
import time

import pytest
from Bio.Seq import Seq
from hypothesis import given, settings
from hypothesis import strategies as st

import p2d
from oracle import region_hits
from p2d.host import _parse
from p2d.optimize import CodonStrategy, OptimizationRequest, optimize_coding

ALL_20 = "ACDEFGHIKLMNPQRSTVWY"
HARD = ["GAATTC", "CTCGAG", "GGTCTC", "CATATG"]


@pytest.fixture(scope="module")
def host():
    return p2d.load_host("ecoli_bl21")


def req(host, aa, **kw):
    kw.setdefault("prefix", "")
    kw.setdefault("suffix", "")
    kw.setdefault("forbidden", [])
    kw.setdefault("strategy", CodonStrategy.MAX_CAI)
    kw.setdefault("seed", 1)
    return OptimizationRequest(aa=aa, host=host, **kw)


def kept(result):
    return [i for i in result.issues if i.code == "soft_site_kept"]


def translate(dna):
    return str(Seq(dna).translate())


# ---------------------------------------------------------------------- the host


def test_soft_motifs_come_from_the_profile(host):
    soft = host.soft_motifs()
    assert soft["AGGAGG"].startswith("internal Shine-Dalgarno")
    assert soft["GCTGGTGG"].startswith("Chi site")
    assert {"A" * 9, "T" * 9, "G" * 6, "C" * 6} <= set(soft), "one longer than the allowed run"
    assert all(reason for reason in soft.values())


def test_a_profile_without_limits_asks_for_nothing(tmp_path):
    bare = tmp_path / "bare.yaml"
    bare.write_text(
        "name: Bare\nslug: bare\ncodon_usage:\n  verified: true\n  table: {M: {ATG: 1.0}}\n"
    )
    assert p2d.load_host(bare).soft_motifs() == {}


def test_a_custom_profile_is_honoured(tmp_path):
    custom = tmp_path / "custom.yaml"
    custom.write_text(
        "name: Custom\nslug: custom\ncodon_usage:\n  verified: true\n  table: {M: {ATG: 1.0}}\n"
        "forbidden_motifs:\n  - {pattern: AATAAA, reason: polyadenylation signal}\n"
        "  - {pattern: ATTTA, reason: AU-rich element}\n"
        "homopolymer_max: {A: 5}\n"
    )
    soft = p2d.load_host(custom).soft_motifs()
    assert soft == {
        "AATAAA": "polyadenylation signal",
        "ATTTA": "AU-rich element",
        "AAAAAA": "run of 6 or more A, longer than this host profile allows",
    }


def test_motif_lookup_is_case_insensitive(tmp_path):
    prof = _parse(
        {
            "name": "x",
            "slug": "x",
            "codon_usage": {"verified": True, "table": {}},
            "forbidden_motifs": [{"pattern": "aggagg"}],
            "homopolymer_max": {"a": 3},
        }
    )
    assert set(prof.soft_motifs()) == {"AGGAGG", "AAAA"}


# ------------------------------------------------------------------ the optimiser


@pytest.mark.parametrize("strategy", [CodonStrategy.MAX_CAI, CodonStrategy.WEIGHTED_SAMPLE])
def test_homopolymer_runs_are_broken_up(host, strategy):
    aa = "K" * 10 + "F" * 10 + "R" * 6  # best codons AAA / TTT / CGC would make 30-base runs
    naive = optimize_coding(req(host, aa, strategy=strategy, forbidden=[]))
    free = re.findall("A{9,}|T{9,}", naive.dna)
    r = optimize_coding(req(host, aa, strategy=strategy, soft_forbidden=host.soft_motifs()))
    assert r.ok and translate(r.dna) == aa and not kept(r)
    assert not re.search("A{9,}|T{9,}|G{6,}|C{6,}", r.dna)
    if strategy is CodonStrategy.MAX_CAI:
        assert free, "premise: without the preference the runs are there"


def test_the_price_of_a_preference_is_visible(host):
    aa = "K" * 10 + "F" * 10
    free = optimize_coding(req(host, aa))
    soft = optimize_coding(req(host, aa, soft_forbidden=host.soft_motifs()))
    assert soft.metrics["cai"] < free.metrics["cai"] == pytest.approx(1.0)
    assert soft.metrics["soft_kept"] == 0


def test_a_site_across_the_junction_is_avoided_when_a_path_exists(host):
    soft = {"AGGAGCGC": "test motif"}  # AGGAG + best Arg codon CGC
    naive = optimize_coding(req(host, "R", prefix="AGGAG"))
    assert naive.dna == "CGC"
    r = optimize_coding(req(host, "R", prefix="AGGAG", soft_forbidden=soft))
    assert r.dna != "CGC" and translate(r.dna) == "R" and not kept(r)
    assert region_hits("AGGAG", r.dna, "", list(soft)) == set()


def test_an_unavoidable_preference_warns_instead_of_failing(host):
    r = optimize_coding(req(host, "AMA", soft_forbidden={"ATG": "test preference"}))
    assert r.ok and r.dna[3:6] == "ATG"
    (w,) = kept(r)
    assert w.severity is p2d.Severity.WARNING and "ATG (test preference)" in w.message
    assert r.metrics["soft_kept"] == 1


def test_only_preferences_that_really_remain_are_reported(host):
    """CTG is avoidable even though ATG is not: it must be honoured, and only ATG named."""
    soft = {"ATG": "forced", "CTG": "avoidable"}
    r = optimize_coding(req(host, "ML", soft_forbidden=soft))
    assert r.dna[:3] == "ATG" and r.dna[3:] != "CTG"
    assert [w.message.split(" ")[0] for w in kept(r)] == ["ATG"]


def test_when_preferences_conflict_the_earlier_ones_win(host):
    leucine = ["TTA", "TTG", "CTC", "CTT", "CTA", "CTG"]  # every Leu codon: all cannot hold
    soft = {p: f"prefer no {p}" for p in leucine}
    r = optimize_coding(req(host, "L", soft_forbidden=soft))
    assert r.ok and r.dna == "CTG", "the last preference is the one given up"
    assert [w.message.split(" ")[0] for w in kept(r)] == ["CTG"]


def test_cloning_sites_always_beat_preferences(host):
    """If a preference would force a cloning site, the preference is the one that yields."""
    r = optimize_coding(
        req(
            host, "L", forbidden=["CTG", "TTA", "TTG", "CTC", "CTA"], soft_forbidden={"CTT": "pref"}
        )
    )
    assert r.ok and r.dna == "CTT" and len(kept(r)) == 1


def test_an_unavoidable_cloning_site_is_still_the_reported_error(host):
    r = optimize_coding(req(host, "AMA", forbidden=["ATG"], soft_forbidden={"GCG": "pref"}))
    assert not r.ok
    (issue,) = [i for i in r.issues if i.code == "unavoidable_site"]
    assert "ATG" in issue.message and "GCG" not in issue.message


def test_same_seed_same_dna_with_preferences(host):
    aa = ALL_20 * 6
    a = optimize_coding(
        req(
            host,
            aa,
            strategy=CodonStrategy.WEIGHTED_SAMPLE,
            seed=9,
            soft_forbidden=host.soft_motifs(),
        )
    )
    b = optimize_coding(
        req(
            host,
            aa,
            strategy=CodonStrategy.WEIGHTED_SAMPLE,
            seed=9,
            soft_forbidden=host.soft_motifs(),
        )
    )
    assert a.dna == b.dna


def test_a_thousand_residues_stay_fast(host):
    aa = (ALL_20 * 50)[:1000]
    t = time.time()
    r = optimize_coding(
        req(
            host,
            aa,
            forbidden=HARD,
            soft_forbidden=host.soft_motifs(),
            strategy=CodonStrategy.WEIGHTED_SAMPLE,
        )
    )
    assert r.ok and time.time() - t < 2.0


# --------------------------------------------- nothing is violated without saying so


@settings(max_examples=120, deadline=None)
@given(
    aa=st.text(alphabet="KFRGAMWLSEC", min_size=1, max_size=60),
    strategy=st.sampled_from([CodonStrategy.MAX_CAI, CodonStrategy.WEIGHTED_SAMPLE]),
    seed=st.integers(0, 2**32),
    prefix=st.text(alphabet="ACGT", max_size=12),
    suffix=st.text(alphabet="ACGT", max_size=12),
)
def test_property_no_silent_violation(host, aa, strategy, seed, prefix, suffix):
    soft = host.soft_motifs()
    r = optimize_coding(
        req(
            host,
            aa,
            strategy=strategy,
            seed=seed,
            prefix=prefix,
            suffix=suffix,
            forbidden=HARD,
            soft_forbidden=soft,
        )
    )
    assert translate(r.dna) == aa
    if not r.ok:  # a cloning site the protein forces, and nothing else
        assert all(i.code == "unavoidable_site" for i in r.issues if i.severity.value == "error")
        return
    assert region_hits(prefix, r.dna, suffix, HARD) == set()
    named = " ".join(w.message for w in kept(r))
    for pattern in soft:
        if region_hits(prefix, r.dna, suffix, [pattern]):
            assert pattern in named, f"{pattern} is in the DNA but was not reported"


# ------------------------------------------------------------------ the whole pipeline


def test_design_insert_applies_the_host_preferences():
    plan = p2d.ConstructPlan.single_chain("K" * 12 + "F" * 12 + "E" * 6)
    v = p2d.load_bundled()
    a = p2d.assemble(p2d.design_insert(plan, v, v.sites["NdeI-XhoI"], p2d.load_host("ecoli_bl21")))
    assert a.ok
    assert not re.search("A{9,}|T{9,}|G{6,}|C{6,}", a.design.coding)
    assert not [i for i in a.issues if i.code == "soft_site_kept"]
    assert a.design.metrics["soft_kept"] == 0


def test_real_protein_is_unchanged_by_the_preferences():
    """Golden proteins are protein-level, so the codon change must not touch them."""
    a = p2d.run("KVFLDWINEAYQRGTRVLAEMAKRGDEFVKRLIAEGHDPFEVLKELGYSE")
    assert a.ok and a.protein.startswith("MGSSHHHHHHSSGLVPRGSHM" + "KVFLDW")
