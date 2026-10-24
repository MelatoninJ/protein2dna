"""Every bundled host profile: a cited, complete codon table, and designs that use it."""

from __future__ import annotations

import dataclasses
import re
import sys
import warnings
from pathlib import Path

import pytest
from Bio.Seq import Seq

import p2d
from oracle import region_hits
from p2d.optimize import CodonStrategy, OptimizationRequest, optimize_coding

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "tools"))
import build_codon_table as tool  # noqa: E402

HOSTS = p2d.available_hosts()
PROTEIN = "KVFLDWINEAYQRGTRVLAEMAKRGDEFVKRLIAEGHDPFEVLKELGYSEHNQPST"


def load(slug):
    with warnings.catch_warnings():
        warnings.simplefilter("error")  # an unverified table would warn
        return p2d.load_host(slug)


def test_the_expected_hosts_ship():
    assert {"ecoli_bl21", "hsapiens", "mmusculus", "scerevisiae", "bsubtilis"} <= set(HOSTS)
    assert set(HOSTS) == set(tool.HOSTS), "tools/build_codon_table.py knows every profile"


@pytest.mark.parametrize("slug", HOSTS)
class TestEveryHost:
    def test_the_table_is_cited_and_verified(self, slug):
        h = load(slug)
        usage = h.raw["codon_usage"]
        assert h.codon_table_verified and h.slug == slug
        table_name = tool.HOSTS[slug][1]
        assert "python_codon_tables" in usage["source"] and table_name in usage["source"]
        assert usage["note"]

    def test_the_table_is_exactly_what_the_cited_package_holds(self, slug):
        """A hand edit of a generated block would show up here."""
        import python_codon_tables as pct

        source = pct.get_codons_table(tool.HOSTS[slug][1])
        h = load(slug)
        for aa, codons in source.items():
            mine = h.stop_usage if aa == "*" else h.codon_usage[aa]
            assert mine == {c: pytest.approx(f, abs=0.005) for c, f in codons.items()}, aa

    def test_the_table_is_complete_and_correct(self, slug):
        h = load(slug)
        seen = [c for codons in h.codon_usage.values() for c in codons]
        assert len(seen) == len(set(seen)) == 61
        for aa, codons in h.codon_usage.items():
            assert all(str(Seq(c).translate()) == aa for c in codons), aa
            assert sum(codons.values()) == pytest.approx(1.0, abs=0.03), aa
        assert sum(h.stop_usage.values()) == pytest.approx(1.0, abs=0.03)
        assert all(str(Seq(c).translate()) == "*" for c in h.stop_usage)

    def test_every_residue_has_a_codon(self, slug):
        h = load(slug)
        assert all(h.best_codon(aa) for aa in "ACDEFGHIKLMNPQRSTVWY")
        assert h.best_codon("M") == "ATG" and h.best_codon("W") == "TGG"

    def test_the_limits_load(self, slug):
        lim = load(slug).limits()
        assert 0 < lim.gc_min < lim.gc_max < 1 and lim.window_size >= 10

    def test_a_design_keeps_its_protein_and_no_site(self, slug):
        h = load(slug)
        sites = ["GAATTC", "CTCGAG", "CATATG"]
        for strategy in (CodonStrategy.MAX_CAI, CodonStrategy.WEIGHTED_SAMPLE):
            r = optimize_coding(
                OptimizationRequest(
                    aa=PROTEIN,
                    host=h,
                    prefix="CATATG",
                    suffix="CTCGAG",
                    forbidden=sites,
                    soft_forbidden=h.soft_motifs(),
                    strategy=strategy,
                    limits=h.limits(),
                )
            )
            assert r.ok and str(Seq(r.dna).translate()) == PROTEIN
            assert not region_hits("CATATG", r.dna, "CTCGAG", sites)

    def test_maximum_cai_follows_that_hosts_preference(self, slug):
        h = load(slug)
        r = optimize_coding(
            OptimizationRequest(
                aa=PROTEIN,
                host=h,
                prefix="",
                suffix="",
                forbidden=[],
                strategy=CodonStrategy.MAX_CAI,
            )
        )
        assert r.metrics["cai"] == pytest.approx(1.0)
        assert [r.dna[i : i + 3] for i in range(0, len(r.dna), 3)] == [
            h.best_codon(a) for a in PROTEIN
        ]


@pytest.mark.parametrize("slug", [s for s in HOSTS if s != "ecoli_bl21"])
def test_rare_codons_are_derived_from_the_table_for_generated_hosts(slug):
    h = load(slug)
    assert sorted(h.rare_codons) == tool.rare_codons(slug)
    for codon in h.rare_codons:
        aa = str(Seq(codon).translate())
        assert h.codon_usage[aa][codon] < tool.RARE_BELOW and len(h.codon_usage[aa]) > 1


def test_hosts_prefer_different_codons():
    best = {s: "".join(load(s).best_codon(a) for a in "KNIRLSG") for s in HOSTS}
    assert len(set(best.values())) >= 4, best


def test_mammalian_hosts_avoid_polyadenylation_and_au_rich_motifs():
    for slug in ("hsapiens", "mmusculus"):
        h = load(slug)
        patterns = set(h.soft_motifs())
        assert {"AATAAA", "ATTAAA", "ATTTA"} <= patterns
    # a protein whose AT-rich codons would spell the motifs; prefer those codons to prove the point
    h = load("hsapiens")
    usage = {aa: dict(c) for aa, c in h.codon_usage.items()}
    for aa, favoured in (("N", "AAT"), ("K", "AAA"), ("I", "ATT")):
        usage[aa] = {c: (0.9 if c == favoured else 0.1 / (len(usage[aa]) - 1)) for c in usage[aa]}
    biased = dataclasses.replace(h, codon_usage=usage)
    plain = dataclasses.replace(biased, forbidden_motifs=[], homopolymer_max={})
    protein = "NKINKINKINKINK" * 3

    def design(host):
        return optimize_coding(
            OptimizationRequest(
                aa=protein,
                host=host,
                prefix="",
                suffix="",
                forbidden=[],
                soft_forbidden=host.soft_motifs(),
                strategy=CodonStrategy.MAX_CAI,
            )
        ).dna

    motifs = ["AATAAA", "ATTAAA", "ATTTA"]
    assert region_hits("", design(plain), "", motifs), "the codons would spell a motif"
    assert not region_hits("", design(biased), "", motifs), "the host profile avoids them"


BACTERIA = ("ecoli_bl21", "bsubtilis")


def test_profiles_say_what_is_and_is_not_cited():
    for slug in set(HOSTS) - set(BACTERIA):
        text = (Path(p2d.__file__).parent / "hosts" / f"{slug}.yaml").read_text()
        assert re.search(r"Only the codon table .*cited source", text.replace("\n# ", " "), re.S)


def test_bacteria_are_released_and_every_other_host_is_in_development():
    released = {s for s in HOSTS if load(s).status == "released"}
    assert released == {"ecoli_bl21", "bsubtilis"}
    assert all(load(s).status == "development" for s in set(HOSTS) - released)
    for slug in HOSTS:
        assert load(slug).status == tool.STATUS[slug], slug


def test_a_profile_without_a_status_counts_as_released(tmp_path):
    path = tmp_path / "mine.yaml"
    path.write_text((Path(p2d.__file__).parent / "hosts" / "ecoli_bl21.yaml").read_text())
    assert p2d.load_host(path).status == "released"


def test_a_design_for_a_host_in_development_says_so_and_still_works():
    vector = p2d.load_bundled()
    for slug, warned in (("ecoli_bl21", False), ("bsubtilis", False), ("hsapiens", True)):
        a = p2d.run("KVFLDWINEAYQRGTRVLAEM", host=slug, vector=vector, site="NdeI-XhoI")
        codes = [i.code for i in a.issues]
        assert ("host_in_development" in codes) is warned and a.ok, slug


def test_the_command_line_marks_hosts_in_development():
    from typer.testing import CliRunner

    from p2d.cli import app

    out = CliRunner().invoke(app, ["hosts"]).output.splitlines()
    marked = {line.split()[0] for line in out if "(in development)" in line}
    assert marked == set(HOSTS) - {"ecoli_bl21", "bsubtilis"}
