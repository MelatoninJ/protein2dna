"""lacZ-alpha start recognised from the sequence of an unannotated plasmid.

pUC19 (L09137) has no features and lacZ-alpha lies on the minus strand. Tests on the real
plasmid skip without ``P2D_VECTOR_SAMPLES``.
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


def test_motif_is_start_of_lacz_alpha():
    assert str(Seq(ALPHA).translate()) == "MTMITPS" and ALPHA.startswith("ATGACCATGATTACG")


def test_start_is_found_on_plus_strand():
    a = analyze_vector(plasmid(1))
    assert not a.flipped and a.ready
    best = a.best_start
    assert best.confidence == "high" and "lacZ-alpha" in best.evidence[0]
    assert a.vector.seq[best.position : best.position + 21] == ALPHA
    assert best.n_terminal.startswith("MTMITPSLHA")


def test_plasmid_stored_on_minus_strand_is_flipped():
    a = analyze_vector(plasmid(-1))
    assert a.flipped and a.ready
    assert any("lacZ-alpha" in n and "minus" in n for n in a.orientation_notes)
    assert a.vector.seq[a.best_start.position : a.best_start.position + 21] == ALPHA


def test_flipped_plasmid_is_reverse_complement_of_input():
    original = plasmid(-1)
    a = analyze_vector(original)
    assert a.vector.seq == rc(original.seq)


def test_sequence_without_lacz_has_no_start():
    rng = random.Random(5)
    seq = "".join(rng.choice("ACGT") for _ in range(900))
    seq = seq.replace(ALPHA, "").replace(rc(ALPHA), "")
    a = analyze_vector(Vector.from_sequence(seq, name="plain"))
    assert not a.ready and not a.starts and not a.flipped


def test_one_base_change_is_not_recognised():
    bad = ALPHA[:10] + ("A" if ALPHA[10] != "A" else "C") + ALPHA[11:]
    seq = "ACGT" * 100 + bad + "TTGCATGCC" + "ACGT" * 100
    assert not analyze_vector(Vector.from_sequence(seq, name="x")).starts


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


def test_real_puc19_oriented_from_sequence():
    a = analyze_vector(vectorio.read_vector_file(sample("pUC19_L09137.gb")))
    assert a.flipped and a.ready and len(a.starts) == 1
    assert a.best_start.n_terminal.startswith("MTMITPSLHACRSTLED")


def test_real_puc19_lists_mcs_enzymes_in_lacz_order(puc19):
    from p2d.enzymes import list_cut_sites

    start = puc19.expression_start
    names = [c.enzyme for c in list_cut_sites(puc19, (start, start + 120))]
    assert names == [
        "HindIII", "SphI", "SbfI", "PstI", "SalI", "XbaI",
        "BamHI", "XmaI", "Acc65I", "KpnI", "SacI", "EcoRI",
    ]  # fmt: skip


def test_gfp_into_puc19_gives_fusion_hand_design_predicts(puc19, gfp):
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


def test_pcr_made_with_hand_designed_primers_agrees_with_p2d(puc19, gfp):
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
    # one stop codon, like the hand design, so the two plasmids are the same size
    assert len(ours.seq) == len(made)


def kozak_plasmid(context="GCCACCATGG", orf_codons=120):
    rng = random.Random(7)
    rand = lambda n: "".join(rng.choice("ACGT") for _ in range(n))  # noqa: E731
    orf = context[6:] + "CT" + "GCT" * orf_codons + "TAA"  # ATG GCT GCT ... TAA, in frame
    return Vector.from_sequence(rand(300) + context[:6] + orf + rand(300), name="mammalian")


def test_kozak_start_is_medium():
    a = analyze_vector(kozak_plasmid())
    assert len(a.starts) == 1 and a.starts[0].confidence == "medium"
    assert "Kozak" in a.starts[0].evidence[0]
    assert a.vector.seq[a.starts[0].position : a.starts[0].position + 4] == "ATGG"
    assert not a.ready, "a medium start must be confirmed by the user"


def test_kozak_start_needs_full_consensus_and_real_orf():
    assert not analyze_vector(kozak_plasmid("GCCTCCATGG")).starts, "position -3 must be A or G"
    assert not analyze_vector(kozak_plasmid("GCCACCATGC")).starts, "G at +4 is part of the motif"
    assert not analyze_vector(kozak_plasmid(orf_codons=10)).starts, "too short to be a gene"


def test_vectors_own_bases_are_tried_first_as_pad():
    from p2d.frame import _pad_candidates

    pads = _pad_candidates(2, prefix="CTCGAG", suffix="", forbidden=[], frame_offset=0, native="CT")
    assert pads[0] == "CT" and len(pads) == len(set(pads))
    odd = _pad_candidates(2, prefix="CTCGAG", suffix="", forbidden=[], frame_offset=0, native="NN")
    assert "NN" not in odd, "only real bases count as native"


def test_altered_counts_changed_residues():
    from p2d.recommend import _altered

    v = p2d.load_bundled()
    native = str(Seq(v.seq[v.expression_start :]).translate(to_stop=True))
    assert _altered(v, native[:12]) == 0
    edited = native[:5] + ("A" if native[5] != "A" else "G") + native[6:12]
    assert _altered(v, edited) == 1


@pytest.fixture(scope="module")
def orange():
    a = analyze_vector(vectorio.read_vector_file(sample("mOrange2-C1.fa")))
    return a, confirm(a, start=a.starts[0].position)


@pytest.fixture(scope="module")
def mek1():
    from Bio import SeqIO

    rec = SeqIO.read(sample("MEK1_BC089772.gb"), "genbank")
    cds = next(f for f in rec.features if f.type == "CDS").location.extract(rec.seq)
    return str(cds)


def test_real_morange2_start_is_kozak_atg_of_morange2(orange):
    a, _ = orange
    assert [s.position for s in a.starts] == [645] and a.starts[0].confidence == "medium"
    assert a.starts[0].n_terminal.startswith("MVSKGEENNMAIIKEFMRFK")


def test_mek1_into_morange2_keeps_fusion_partner_intact(orange, mek1):
    _, v = orange
    protein = str(Seq(mek1).translate(to_stop=True))
    assert len(protein) == 393
    host = p2d.load_host("ecoli_bl21")
    options, _ = rank_pairs(v, protein, host, window=(1330, 1445))
    assert options and all(o.altered == 0 for o in options), "no pair may change mOrange2"
    by_pair = {(o.up.enzyme, o.down.enzyme): o for o in options}
    xho = by_pair[("XhoI", "BamHI")]
    assert xho.before.endswith("ELYKSGLRSRA") and xho.pad_bases == 2
    assert xho.before[-7:] == "SGLRSRA", "the seven residues the hand design predicts"
    groups = group_by_outcome(options)
    first = groups[0].best
    assert {first.up.enzyme, first.down.enzyme} == {"Bsp1407I", "AccIII"} or first.before.endswith(
        "ELYK"
    )
    assert groups[0].before.endswith("MDELYK"), "BsrGI keeps mOrange2 and adds nothing"


def test_hand_designed_mek1_primers_agree_with_p2d(orange, mek1):
    from Bio.Restriction import BamHI, BsrGI, XhoI

    _, v = orange
    mek = str(Seq(mek1).translate(to_stop=True))
    host = p2d.load_host("ecoli_bl21")

    def by_hand(fwd, rev, up, down, n_f, n_r):
        assert mek1.startswith(fwd[-n_f:]) and mek1.endswith(rc(rev[-n_r:]))
        product = fwd[:-n_f] + mek1 + rc(rev[:-n_r])
        cut = lambda enz, s: enz.search(Seq(s), linear=True)  # noqa: E731
        frag = product[cut(up, product)[0] - 1 : cut(down, product)[0] - 1]
        made = v.seq[: cut(up, v.seq)[0] - 1] + frag + v.seq[cut(down, v.seq)[0] - 1 :]
        read = made[v.expression_start :]
        return made, str(Seq(read[: len(read) - len(read) % 3]).translate(to_stop=True))

    from p2d.enzymes import list_cut_sites, plan_cut

    sites = {}
    for c in list_cut_sites(v, (1330, 1445)):
        for name in (c.enzyme, *c.aliases):
            sites.setdefault(name, c)

    cases = [
        ("GCGCCTCGAGCTATGCCCAAGAAGAAGCC", "GCGCGGATCCTCAGATGCTGGCAGCGT", XhoI, BamHI, 17, 17),
        ("GCGCTGTACAAGATGCCCAAGAAGAAGCCG", "GCGCGGATCCTCAGATGCTGGCAGCGT", BsrGI, BamHI, 18, 17),
    ]
    for fwd, rev, up, down, n_f, n_r in cases:
        made, protein = by_hand(fwd, rev, up, down, n_f, n_r)
        assert protein.endswith(mek) and len(protein) - len(mek) in (236 + 7, 236)
        site = plan_cut(v, sites[str(up)], sites[str(down)])
        ours = p2d.assemble(
            p2d.design_insert(
                p2d.ConstructPlan.single_chain(mek), v, site, host, cterm=p2d.CTerm.STOP
            )
        )
        assert ours.ok and ours.protein == protein
        assert len(ours.seq) == len(made), "one stop codon, like the hand design"
