"""How annotations help find a vector's start codon and direction, and where they must not."""

from __future__ import annotations

import random

import pytest
from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord

from p2d.analysis import analyze_vector
from p2d.vector import Vector

CODONS = ["GCT", "AAA", "GAA", "CTG", "GGC", "GTT", "ACC", "CAG"]  # no stop codon among them


def record(parts, features, name="v"):
    seq = "".join(parts)
    rec = SeqRecord(Seq(seq), id=name, name=name, annotations={"molecule_type": "DNA"})
    rec.features = features
    return rec


def cds(start, end, label, strand=1):
    return SeqFeature(
        FeatureLocation(start, end, strand=strand), type="CDS", qualifiers={"product": [label]}
    )


def filler(n, seed=0):
    rng = random.Random(seed)
    return "".join(rng.choice("ACGT") for _ in range(n)).replace("ATG", "ACG")


def vector_with(label, *, strand=1, atg=True, extra=(), extra_seq=""):
    """A plasmid with one annotated CDS of about 750 bp and no cloning region."""
    rng = random.Random(7)
    gene = ("ATG" if atg else "CCC") + "".join(rng.choice(CODONS) for _ in range(250))
    left = filler(500, 1)
    body = gene + "TAA"
    seq = (
        left
        + (body if strand == 1 else str(Seq(body).reverse_complement()))
        + filler(500, 2)
        + extra_seq
    )
    feats = [cds(500, 500 + len(body), label, strand), *extra]
    return Vector(record=record([seq], feats), name="synthetic", circular=True)


def test_partner_cds_offered_without_cloning_region():
    a = analyze_vector(vector_with("glutathione S-transferase"))
    best = a.best_start
    assert best and best.position == 500 and best.confidence == "medium"
    assert "no cloning region is annotated" in best.evidence[0]
    assert not a.ready, "a medium start has to be confirmed"


@pytest.mark.parametrize(
    "label", ["beta-lactamase", "LacI repressor", "kanamycin resistance", "URA3"]
)
def test_other_cds_are_never_offered_without_cloning_region(label):
    assert analyze_vector(vector_with(label)).best_start is None


def test_lacz_alpha_or_male_fusion_is_not_mistaken_for_marker():
    a = analyze_vector(vector_with("MalE-lacZalpha fusion protein"))
    assert a.best_start and a.best_start.position == 500


def test_cds_without_atg_not_offered():
    assert analyze_vector(vector_with("glutathione S-transferase", atg=False)).best_start is None


def test_partner_cds_direction_flips_plasmid():
    v = vector_with("glutathione S-transferase", strand=-1)
    a = analyze_vector(v)
    assert a.flipped and any("fusion-partner" in n for n in a.orientation_notes)
    assert a.best_start and a.best_start.position < len(v.seq) - 700, "now read forward"
    assert not analyze_vector(vector_with("glutathione S-transferase")).flipped


def test_stray_t7_does_not_override_partner_cds():
    """pGATA carries a T7 promoter sequence on the other strand than its GST."""
    t7 = "TAATACGACTCACTATAG"
    v = vector_with("glutathione S-transferase", extra_seq=str(Seq(t7).reverse_complement()))
    assert t7 not in v.seq and not analyze_vector(v).flipped


def test_short_or_unlabelled_cds_does_not_vote_on_direction():
    short = Vector(
        record=record(
            [filler(400, 3), "ATG" + "GCT" * 30 + "TAA", filler(400, 4)], [cds(400, 493, "GST", -1)]
        ),
        name="short",
        circular=True,
    )
    assert not analyze_vector(short).flipped
    unlabelled = vector_with("hypothetical protein", strand=-1)
    assert not analyze_vector(unlabelled).flipped


def test_resistance_gene_is_not_the_start():
    rng = random.Random(3)
    gene = "ATG" + "".join(rng.choice(CODONS) for _ in range(433)) + "TAA"  # reads into the site
    seq = filler(100, 5) + gene + filler(600, 6)
    site = SeqFeature(
        FeatureLocation(1300, 1350, strand=1),
        type="misc_feature",
        qualifiers={"note": ["multiple cloning site"]},
    )
    feats = [cds(100, 100 + len(gene), "beta-lactamase"), site]
    v = Vector(record=record([seq], feats), name="m", circular=True)
    a = analyze_vector(v)
    assert a.region is not None, "the site is recognised"
    evidence = a.best_start.evidence if a.best_start else ()
    assert not any("annotated CDS" in e for e in evidence), "a resistance gene is not a start"
    # the same gene labelled as a fusion partner would have been offered
    feats[0] = cds(100, 100 + len(gene), "glutathione S-transferase")
    again = analyze_vector(Vector(record=record([seq], feats), name="m", circular=True))
    assert again.best_start and again.best_start.position == 100
