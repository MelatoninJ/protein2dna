"""Expression vectors with no annotation: the start codon is read from the sequence.

An ATG counts when a ribosome-binding-site-like motif lies 4-10 bases before it and the ATG is
part of an NcoI or NdeI site, or sits behind a T7 promoter or a lac operator.  The real records
(NCBI pTrc99A, pBAD24, pET-3a, pET-11a, pET-30a, pColdI, and four that must give nothing) live in
``unannotated/`` under ``P2D_VECTOR_SAMPLES`` and the tests skip without them.
"""

from __future__ import annotations

import os
import random
from pathlib import Path

import pytest
from Bio.Seq import Seq

from p2d import vectorio
from p2d.analysis import analyze_vector
from p2d.vector import Vector

SAMPLES = Path(os.environ["P2D_VECTOR_SAMPLES"]) if os.environ.get("P2D_VECTOR_SAMPLES") else None
LAC_O = "TTGTGAGCGGATAACAA"
T7 = vectorio.T7_PROMOTER
CODONS = ["GCT", "AAA", "GAA", "CTG", "GGC", "GTT", "ACC", "CAG"]


def rc(s):
    return str(Seq(s).reverse_complement())


def flank(n, seed):
    rng = random.Random(seed)
    return "".join(rng.choice("ACGT") for _ in range(n)).replace("ATG", "ACG")


def cassette(before, *, rbs="AGGAGG", gap="TATATAC", start="ATG", tail="GCTAGC", orf=0, seed=1):
    """Flank + ``before`` + RBS + gap + start codon (+ site tail) + a stretch of codons + stop."""
    body = "".join(random.Random(seed).choice(CODONS) for _ in range(orf))
    return (
        flank(300, seed)
        + before
        + flank(20, seed + 1)
        + rbs
        + gap
        + start
        + tail
        + body
        + "TAA"
        + flank(300, seed + 2)
    )


def analyse(seq):
    return analyze_vector(Vector.from_sequence(seq, name="synthetic"))


# ------------------------------------------------------------------------ synthetic


def test_an_atg_inside_an_ncoi_site_after_an_rbs_is_the_start():
    seq = cassette("", rbs="AGGAGG", gap="AATTCAC", start="CATGG", orf=5).replace(
        "CATGG", "CCATGG", 1
    )
    a = analyse(seq)
    best = a.best_start
    assert best and "NcoI" in best.evidence[0] and best.confidence == "medium"
    assert seq[best.position - 2 : best.position + 4] == "CCATGG"


def test_an_atg_inside_an_ndei_site_is_found_however_short_the_frame():
    seq = cassette(T7, gap="TATAC", start="ATATG", tail="GCTAGCTAA", orf=0).replace(
        "ATATG", "CATATG", 1
    )
    best = analyse(seq).best_start
    assert best and "NdeI" in best.evidence[0]
    assert seq[best.position - 3 : best.position + 3] == "CATATG"
    assert best.confidence == "high"  # site and T7 promoter agree


def test_a_promoter_without_a_site_gives_a_medium_start_needing_a_long_frame():
    long = analyse(cassette(LAC_O, orf=40))
    assert long.best_start.confidence == "medium" and "lac operator" in long.best_start.evidence[0]
    assert not analyse(cassette(LAC_O, orf=3)).starts  # a 12-codon frame is not an expression start


def test_an_rbs_or_atg_alone_is_not_enough():
    assert not analyse(cassette("", orf=40)).starts  # no site, no promoter
    far = cassette(LAC_O, orf=40).replace(LAC_O, LAC_O + flank(300, 9))
    assert not analyse(far).starts  # the operator is more than 250 bases away
    assert not analyse(cassette(LAC_O, rbs="TTTTTT", orf=40)).starts  # no ribosome-binding site


def test_the_rbs_has_to_sit_at_a_sensible_distance_from_the_atg():
    assert not analyse(cassette(LAC_O, gap="A" * 14, orf=40)).starts
    assert analyse(cassette(LAC_O, gap="A" * 5, orf=40)).starts


def test_a_plasmid_stored_backwards_is_flipped_by_its_cassette():
    seq = cassette(LAC_O, orf=40)
    plus, minus = analyse(seq), analyse(rc(seq))
    assert not plus.flipped and minus.flipped
    assert minus.vector.seq == plus.vector.seq
    assert minus.best_start.position == plus.best_start.position
    assert any("minus strand" in n for n in minus.orientation_notes)


def test_the_cassette_does_not_outvote_a_t7_promoter_or_lacz_alpha():
    seq = cassette(T7, orf=40)
    votes = vectorio.orientation_votes(Vector.from_sequence(seq, name="t7"))
    assert votes[0] == 2 and votes[1] == 0  # the T7 promoter alone; the cassette adds nothing


# --------------------------------------------------------------------- real records


def record(name):
    if SAMPLES is None or not (SAMPLES / "unannotated" / name).is_file():
        pytest.skip(f"unannotated/{name} not available (set P2D_VECTOR_SAMPLES)")
    return analyze_vector(vectorio.read_vector_file(SAMPLES / "unannotated" / name))


@pytest.mark.parametrize(
    "name, start, flipped, site, confidence, n_terminal",
    [
        ("U13872.1.gb", 266, False, "NcoI", "high", "MEFELGTRGS"),  # pTrc99A
        ("X81837.1.gb", 1319, False, "NcoI", "medium", "MVPGDPLEST"),  # pBAD24
        ("LT727576.1.gb", 4089, True, "NdeI", "high", "MASMTGGQQMGR"),  # pET-3a
        ("LT727489.1.gb", 5317, True, "NdeI", "high", "MASMTGGQQMGR"),  # pET-11a
        ("PZ730019.1.gb", 4964, True, "", "high", "MHHHHHHSSGLVPR"),  # pET-30a
        ("AB186388.1.gb", 255, False, "", "high", "MNHKVHHHHHH"),  # pColdI
    ],
)
def test_real_unannotated_expression_vectors(name, start, flipped, site, confidence, n_terminal):
    a = record(name)
    best = a.best_start
    assert best.position == start and a.flipped == flipped and best.confidence == confidence
    assert a.vector.seq[start : start + 3] == "ATG" and best.n_terminal.startswith(n_terminal)
    if site:
        assert site in best.evidence[0]


@pytest.mark.parametrize(
    "name",
    [
        "M77749.1.gb",  # pKK223-3: the insert brings its own ATG
        "AY230150.1.gb",  # pT7-5: T7 promoter, no ATG in the vector
        "X81838.1.gb",  # pBAD18
        "LC760188.1.gb",  # pBAD33
    ],
)
def test_real_vectors_without_a_start_give_none_and_ask_for_one(name):
    a = record(name)
    assert not a.starts and not a.ready
    assert any("point to one" in w for w in a.warnings)
