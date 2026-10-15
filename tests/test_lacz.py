"""Plasmids that carry no annotation: the lacZ-alpha start is recognised from its sequence.

pUC19 in GenBank (L09137) has no features at all and lacZ-alpha lies on the minus strand, so
nothing but the sequence says where translation starts or which way the plasmid should read.
The real-plasmid tests need files the repository cannot hold and skip without them.
"""

from __future__ import annotations

import os
import random
from pathlib import Path

import pytest
from Bio.Restriction import SphI, XbaI
from Bio.Seq import Seq

import p2d
from p2d import vectorio
from p2d.analysis import analyze_vector, confirm
from p2d.recommend import group_by_outcome, rank_pairs
from p2d.vector import Vector

SAMPLES = Path(os.environ["P2D_VECTOR_SAMPLES"]) if os.environ.get("P2D_VECTOR_SAMPLES") else None
ALPHA = vectorio.LACZ_ALPHA_START


def rc(s):
    return str(Seq(s).reverse_complement())


def plasmid(strand, seed=1):
    """Random backbone with the start of lacZ-alpha and a short ORF after it."""
    rng = random.Random(seed)
    rand = lambda n: "".join(rng.choice("ACGT") for _ in range(n))  # noqa: E731
    orf = ALPHA + "TTGCATGCCTAGGAATTC" + "GCT" * 40 + "TAA"
    # nothing random may contain a stop in frame with the start, so keep the flank away from it
    seq = rand(400) + (orf if strand == 1 else rc(orf)) + rand(400)
    return Vector.from_sequence(seq, name="synthetic")


def test_the_motif_is_the_start_of_lacz_alpha():
    assert str(Seq(ALPHA).translate()) == "MTMITPS" and ALPHA.startswith("ATGACCATGATTACG")


def test_the_start_is_found_on_the_plus_strand():
    a = analyze_vector(plasmid(1))
    assert not a.flipped and a.ready
    best = a.best_start
    assert best.confidence == "high" and "lacZ-alpha" in best.evidence[0]
    assert a.vector.seq[best.position : best.position + 21] == ALPHA
    assert best.n_terminal.startswith("MTMITPSLHA")


def test_a_plasmid_stored_on_the_minus_strand_is_flipped():
    a = analyze_vector(plasmid(-1))
    assert a.flipped and a.ready
    assert any("lacZ-alpha" in n and "minus" in n for n in a.orientation_notes)
    assert a.vector.seq[a.best_start.position : a.best_start.position + 21] == ALPHA


def test_the_flipped_plasmid_is_the_reverse_complement_of_the_input():
    original = plasmid(-1)
    a = analyze_vector(original)
    assert a.vector.seq == rc(original.seq)


def test_a_sequence_without_lacz_is_still_reported_as_having_no_start():
    rng = random.Random(5)
    seq = "".join(rng.choice("ACGT") for _ in range(900))
    seq = seq.replace(ALPHA, "").replace(rc(ALPHA), "")
    a = analyze_vector(Vector.from_sequence(seq, name="plain"))
    assert not a.ready and not a.starts and not a.flipped


def test_a_one_base_change_is_not_recognised():
    bad = ALPHA[:10] + ("A" if ALPHA[10] != "A" else "C") + ALPHA[11:]
    seq = "ACGT" * 100 + bad + "TTGCATGCC" + "ACGT" * 100
    assert not analyze_vector(Vector.from_sequence(seq, name="x")).starts


# ------------------------------------------------------------------- the real pUC19


def sample(name):
    if SAMPLES is None or not (SAMPLES / name).is_file():
        pytest.skip(f"{name} not available (set P2D_VECTOR_SAMPLES)")
    return SAMPLES / name


@pytest.fixture(scope="module")
def puc19():
    return confirm(analyze_vector(vectorio.read_vector_file(sample("pUC19_L09137.gb"))))


@pytest.fixture(scope="module")
def gfp():
    lines = sample("GFP_pGLO.txt").read_text().split("\n")
    return "".join(lines[1:])


def test_real_puc19_is_oriented_and_started_from_sequence_alone():
    a = analyze_vector(vectorio.read_vector_file(sample("pUC19_L09137.gb")))
    assert a.flipped and a.ready and len(a.starts) == 1
    assert a.best_start.n_terminal.startswith("MTMITPSLHACRSTLED")


def test_real_puc19_lists_the_mcs_enzymes_in_lacz_order(puc19):
    from p2d.enzymes import list_cut_sites

    start = puc19.expression_start
    names = [c.enzyme for c in list_cut_sites(puc19, (start, start + 120))]
    assert names == [
        "HindIII", "SphI", "SbfI", "PstI", "SalI", "XbaI",
        "BamHI", "XmaI", "Acc65I", "KpnI", "SacI", "EcoRI",
    ]  # fmt: skip


def test_gfp_into_puc19_gives_the_fusion_a_hand_design_predicts(puc19, gfp):
    """SphI + XbaI, the pair chosen by hand: MTMITPSLHA in front of GFP, one frame base."""
    protein = str(Seq(gfp).translate(to_stop=True))
    host = p2d.load_host("ecoli_bl21")
    start = puc19.expression_start
    options, _ = rank_pairs(puc19, protein, host, window=(start, start + 120))
    chosen = [o for o in options if (o.up.enzyme, o.down.enzyme) == ("SphI", "XbaI")]
    assert len(chosen) == 1 and chosen[0].before == "MTMITPSLHA" and chosen[0].after == ""
    assert chosen[0].pad_bases == 1
    groups = {(g.before, g.pad_bases) for g in group_by_outcome(options)}
    assert ("MTMITPSL", 1) in groups, "HindIII leaves the fewest residues once the gene is recoded"


def test_a_pcr_made_with_the_hand_designed_primers_agrees_with_p2d(puc19, gfp):
    """Cut a PCR product and the plasmid with Biopython, ligate, and compare with p2d."""
    forward = "GCGCGCATGCCATGGCTAGCAAAGGAGA"
    reverse = "GCGCTCTAGATTATTTGTAGAGCTCATCCA"
    assert gfp.startswith(forward[-17:]) and gfp.endswith(rc(reverse[-20:]))
    product = forward[:-17] + gfp + rc(reverse[:-20])
    sph, xba = SphI.search(Seq(product), linear=True)[0], XbaI.search(Seq(product), linear=True)[0]
    v_sph, v_xba = SphI.search(Seq(puc19.seq))[0], XbaI.search(Seq(puc19.seq))[0]
    made = puc19.seq[: v_sph - 1] + product[sph - 1 : xba - 1] + puc19.seq[v_xba - 1 :]
    read = made[puc19.expression_start :]
    protein = str(Seq(read[: len(read) - len(read) % 3]).translate(to_stop=True))
    assert protein == "MTMITPSLHA" + str(Seq(gfp).translate(to_stop=True))

    from p2d.enzymes import list_cut_sites, plan_cut

    sites = {c.enzyme: c for c in list_cut_sites(puc19)}
    site = plan_cut(puc19, sites["SphI"], sites["XbaI"])
    protein = str(Seq(gfp).translate(to_stop=True))
    plan = p2d.ConstructPlan.single_chain(protein)
    host = p2d.load_host("ecoli_bl21")
    ours = p2d.assemble(p2d.design_insert(plan, puc19, site, host, cterm=p2d.CTerm.STOP))
    assert ours.ok and ours.protein == "MTMITPSLHA" + protein
    # p2d ends with a double stop (TAATGA); the hand design has a single one
    assert len(ours.seq) - len(made) == 3
