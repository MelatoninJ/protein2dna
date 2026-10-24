"""The bacterial limits come from the host's genome, and each cited value can be traced."""

from __future__ import annotations

import os
import random
import sys
from pathlib import Path

import pytest
import yaml
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord

import p2d

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "tools"))
import build_limits as bl  # noqa: E402

SAMPLES = Path(os.environ["P2D_VECTOR_SAMPLES"]) if os.environ.get("P2D_VECTOR_SAMPLES") else None
HOSTS = Path(p2d.__file__).parent / "hosts"
BACTERIA = ("ecoli_bl21", "bsubtilis")


def profile(slug):
    return yaml.safe_load((HOSTS / f"{slug}.yaml").read_text())


def genome(n=600, seed=1, extra=()):
    """A record of ``n`` random genes plus ``extra`` genes; every gene is a clean CDS."""
    rng = random.Random(seed)
    seq, feats = "", []
    for g in [
        "ATG" + "".join(rng.choice("ACGT") for _ in range(297)) + "TAA" for _ in range(n)
    ] + list(extra):
        feats.append(SeqFeature(FeatureLocation(len(seq), len(seq) + len(g), strand=1), type="CDS"))
        seq += g
    return SeqRecord(Seq(seq), id="FAKE", features=feats)


def test_a_run_found_in_too_few_genes_does_not_set_the_limit():
    odd = "ATG" + "A" * 30 + "GC" * 60 + "TAA"
    limits = bl.derive(bl.coding_sequences(genome(extra=[odd])))
    plain = bl.derive(bl.coding_sequences(genome()))
    assert limits["homopolymer"]["A"] == plain["homopolymer"]["A"] < 30, "one gene is not enough"


def test_a_run_found_in_three_genes_does_set_the_limit():
    rich = "ATG" + "A" * 15 + "GC" * 66 + "TAA"
    limits = bl.derive(bl.coding_sequences(genome(extra=[rich] * 3)))
    assert limits["homopolymer"]["A"] == 15
    two = bl.derive(bl.coding_sequences(genome(extra=[rich] * 2)))
    assert two["homopolymer"]["A"] < 15


def test_gc_ranges_are_the_central_99_percent():
    d = bl.derive(bl.coding_sequences(genome()))
    assert d["n_cds"] == 600 and 0.35 < d["gc"][0] < 0.5 < d["gc"][1] < 0.65
    assert d["window"][0] < d["gc"][0] and d["window"][1] > d["gc"][1], (
        "windows vary more than genes"
    )


def test_only_clean_coding_sequences_are_used():
    ok = "ATG" + "GCT" * 60 + "TAA"
    bad = [
        "CTG" + "GCT" * 60 + "TAA",
        "ATG" + "GCN" * 60 + "TAA",
        "ATG" + "GCT" * 10 + "TAA",
        "ATG" + "GCT" * 60 + "TA",
    ]
    got = bl.coding_sequences(genome(n=0, extra=[ok, *bad]))
    assert got == [ok]


def test_too_small_a_record_is_not_taken_for_a_genome():
    with pytest.raises(ValueError, match="not a genome"):
        bl.derive(bl.coding_sequences(genome(n=100)))


def test_the_marker_splice_is_idempotent_and_keeps_what_comes_before():
    text = "name: x\n\n" + bl.MARKER + "\nold: 1\n"
    block = bl.render(
        "ecoli_bl21", "NC_1.1", "Test genome", bl.derive(bl.coding_sequences(genome()))
    )
    once = bl.splice(text, block)
    assert once.startswith("name: x\n\n") and bl.splice(once, block) == once
    loaded = yaml.safe_load(once)
    assert loaded["repeats"] == {"direct_max": 22, "inverted_max": 200}
    with pytest.raises(SystemExit):
        bl.splice("no marker here", block)


@pytest.mark.parametrize("slug", BACTERIA)
def test_every_limit_names_where_it_comes_from(slug):
    data = profile(slug)
    src = data["limits_source"]
    assert src["genome"].startswith(bl.GENOMES[slug]) and "coding sequences" in src["genome"]
    assert "Shen & Huang 1986" in src["repeats"] and "percentile" in src["gc"]
    assert data["repeats"] == {"direct_max": bl.DIRECT_REPEAT_MAX, "inverted_max": bl.INVERTED_OFF}
    for motif in data["forbidden_motifs"]:
        assert motif["source"] and motif["reason"], motif
    text = (HOSTS / f"{slug}.yaml").read_text()
    assert bl.MARKER in text and "placeholder" not in text.lower()


def test_the_e_coli_rare_codons_are_the_ones_rosetta_2_supplies():
    assert set(profile("ecoli_bl21")["rare_codons"]) == {
        "AGA",
        "AGG",
        "ATA",
        "CUA".replace("U", "T"),
        "GGA",
        "CCC",
        "CGG",
    }


def test_the_internal_sd_motif_is_listed_for_both_bacteria_and_chi_only_for_e_coli():
    patterns = {s: {m["pattern"] for m in profile(s)["forbidden_motifs"]} for s in BACTERIA}
    assert "AGGAGG" in patterns["ecoli_bl21"] and "AGGAGG" in patterns["bsubtilis"]
    assert "GCTGGTGG" in patterns["ecoli_bl21"] and "GCTGGTGG" not in patterns["bsubtilis"]


def test_the_limits_are_not_stricter_than_nature():
    """The old defaults forbade runs that the host's own genes contain thousands of times."""
    h = p2d.load_host("ecoli_bl21")
    assert min(h.homopolymer_max.values()) >= 7 and h.homopolymer_max["A"] >= 9
    lim = h.limits()
    assert lim.gc_min <= 0.33 and lim.gc_max >= 0.59 and lim.direct_max == 22


@pytest.mark.parametrize("slug", BACTERIA)
def test_the_shipped_numbers_are_what_the_genome_gives(slug):
    path = SAMPLES / "genomes" / f"{bl.GENOMES[slug]}.gb" if SAMPLES else None
    if not path or not path.is_file():
        pytest.skip(f"{bl.GENOMES[slug]}.gb not available (set P2D_VECTOR_SAMPLES)")
    record = SeqIO.read(str(path), "genbank")
    d = bl.derive(bl.coding_sequences(record))
    data = profile(slug)
    assert data["gc"]["global"] == d["gc"] and data["gc"]["window"]["range"] == d["window"]
    h = d["homopolymer"]
    assert data["homopolymer_max"] == {"A": h["A"], "T": h["T"], "G": h["G"], "C": h["C"]}
    assert str(d["n_cds"]) in data["limits_source"]["genome"]
