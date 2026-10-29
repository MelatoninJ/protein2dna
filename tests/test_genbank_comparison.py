"""Rebuilding published constructs from their empty vector (tools/compare_genbank.py).

Synthetic constructs are built by hand with their own pad rule, so they share no frame code with
p2d. Tests on real records skip without ``P2D_VECTOR_SAMPLES`` (see VALIDATION.md).
"""

from __future__ import annotations

import os
import random
import sys
from pathlib import Path

import pytest
from Bio.Data import CodonTable
from Bio.Seq import Seq

import p2d
from p2d.enzymes import list_cut_sites, resolve_site

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "tools"))
import compare_genbank as cg  # noqa: E402

SAMPLES = Path(os.environ["P2D_VECTOR_SAMPLES"]) if os.environ.get("P2D_VECTOR_SAMPLES") else None
PROTEIN = "KVFLDWINEAYQRGTRVLAEMAKRGDEFVKRLIAEGHDPFEVLKELGYSE"
BACK = {}
for codon, aa in CodonTable.unambiguous_dna_by_id[1].forward_table.items():
    BACK.setdefault(aa, []).append(codon)


@pytest.fixture(scope="module")
def vector():
    return p2d.load_bundled()


@pytest.fixture(scope="module")
def host():
    return p2d.load_host("ecoli_bl21")


def coding(protein: str, seed: int = 1) -> str:
    rng = random.Random(seed)
    return "".join(rng.choice(BACK[a]) for a in protein)


def build(vector, up, down, protein, *, pad="GC", stop="TAA", seed=1):
    """A construct made by hand: sites kept, a pad that restores the frame, then the gene."""
    start = vector.expression_start
    n_pad = (3 - (up.end - start) % 3) % 3
    seq = vector.seq
    return seq[: up.end] + (pad * 3)[:n_pad] + coding(protein, seed) + stop + seq[down.start :]


def pair(vector, a="NdeI", b="XhoI"):
    return resolve_site(vector, a), resolve_site(vector, b)


def test_hand_built_construct_reproduced(vector, host):
    up, down = pair(vector)
    for seed in (1, 2, 3):
        r = cg.compare("hand", vector, build(vector, up, down, PROTEIN, seed=seed), host)
        assert r.status == "match", r
        assert r.up == up.id and r.down == down.id and r.cterm == "stop"
        assert r.real_protein == r.p2d_protein and r.real_protein.endswith(PROTEIN)


def test_construct_may_be_stored_backwards_or_rotated(vector, host):
    up, down = pair(vector)
    built = build(vector, up, down, PROTEIN)
    flipped = str(Seq(built).reverse_complement())
    rotated = built[400:] + built[:400]
    for text in (flipped, rotated, str(Seq(rotated).reverse_complement())):
        assert cg.compare("x", vector, text, host).status == "match"


@pytest.mark.parametrize(
    "a,b", [("NdeI", "XhoI"), ("NdeI", "EcoRI"), ("BamHI", "XhoI"), ("NcoI", "XhoI")]
)
def test_several_site_pairs_and_junction_residues_agree(vector, host, a, b):
    up, down = pair(vector, a, b)
    r = cg.compare("pair", vector, build(vector, up, down, PROTEIN), host)
    assert r.status == "match" and r.real_protein == r.p2d_protein


def test_wrong_junction_is_reported_as_mismatch(vector, host, monkeypatch):
    """Wrong frame-filling bases would change the protein; the tool must see that."""
    up, down = pair(vector, "NcoI", "XhoI")  # NcoI needs two pad bases in the test vector
    built = build(vector, up, down, PROTEIN, pad="GC")
    assert cg.compare("ok", vector, built, host).status == "match"

    import p2d.frame as frame

    monkeypatch.setattr(frame, "_pad_candidates", lambda n, *a, **k: ["T" * n] if n else [""])
    r = cg.compare("broken", vector, built, host)
    assert r.status == "mismatch", r.reason
    assert r.real_protein != r.p2d_protein
    assert cg.first_difference(r.real_protein, r.p2d_protein) is not None


def test_insert_that_is_not_read_through_is_not_comparable(vector, host):
    up, down = pair(vector)
    start = vector.expression_start
    n_pad = (3 - (up.end - start) % 3) % 3
    shifted = (
        vector.seq[: up.end]
        + "G" * (n_pad + 1)
        + coding(PROTEIN)
        + "TAA"
        + vector.seq[down.start :]
    )
    r = cg.compare("frameshift", vector, shifted, host)
    assert r.status == "not comparable" and "not read through" in r.reason


def test_construct_with_no_site_at_junction_is_not_comparable(vector, host):
    start = vector.expression_start
    seamless = vector.seq[: start + 3] + coding(PROTEIN) + "TAA" + vector.seq[start + 3 + 300 :]
    r = cg.compare("seamless", vector, seamless, host)
    assert r.status == "not comparable" and "no cloning site" in r.reason


def test_unrelated_or_identical_record_is_not_comparable(vector, host):
    rng = random.Random(9)
    unrelated = "".join(rng.choice("ACGT") for _ in range(3000))
    assert cg.compare("u", vector, unrelated, host).status == "not comparable"
    assert cg.compare("same", vector, vector.seq, host).status == "not comparable"


def test_pair_p2d_refuses_is_declined_not_called_failure(vector, host, monkeypatch):
    up, down = pair(vector)

    def refuse(*args):
        raise ValueError("the sites overlap")

    monkeypatch.setattr(cg, "plan_cut", refuse)
    r = cg.compare("refused", vector, build(vector, up, down, PROTEIN), host)
    assert r.status == "declined" and "overlap" in r.reason


def test_touching_sites_are_rebuilt_like_any_other_pair(vector, host):
    up, down = pair(vector, "NdeI", "BamHI")  # the two sites touch in the test vector
    r = cg.compare("touching", vector, build(vector, up, down, PROTEIN), host)
    assert r.status == "match" and r.real_protein == r.p2d_protein


def test_command_line_reports_and_sets_its_exit_status(vector, tmp_path, capsys):
    from Bio import SeqIO
    from Bio.Seq import Seq as S
    from Bio.SeqRecord import SeqRecord

    up, down = pair(vector)
    good = tmp_path / "good.gb"
    SeqIO.write(
        SeqRecord(
            S(build(vector, up, down, PROTEIN)),
            id="HAND1",
            name="HAND1",
            annotations={"molecule_type": "DNA"},
        ),
        str(good),
        "genbank",
    )
    vec = tmp_path / "vec.gb"
    SeqIO.write(vector.record, str(vec), "genbank")
    code = cg.main([str(vec), str(good), "--start", str(vector.expression_start + 1)])
    out = capsys.readouterr().out
    assert code == 0 and "HAND1" in out and "1 match" in out


# accession -> (5' site, 3' site, residues of the real fusion protein); all made by their authors
# in pET-28a(+) and taken from NCBI.  See VALIDATION.md.
MATCHING = {
    "MK847897.1": ("NcoI:5068", "BamHI:5166", 338),
    "MK847898.1": ("BmtI:5133", "NotI:5197", 317),
    "MK847902.1": ("NcoI:5068", "BamHI:5166", 278),
    "MK847903.1": ("BamHI:5166", "SacI:5178", 300),
    "MK847904.1": ("BamHI:5166", "SacI:5178", 280),
    "OK264409.1": ("NcoI:5068", "NdeI:5127", 210),
    "OL989372.1": ("NdeI:5127", "XhoI:5206", 533),
    "OL989373.1": ("NdeI:5127", "XhoI:5206", 533),
}


@pytest.fixture(scope="module")
def pet28a():
    from p2d.analysis import analyze_vector, confirm
    from p2d.vectorio import read_vector_file

    path = SAMPLES / "pET-28a_plus.dna" if SAMPLES else None
    if not path or not path.is_file():
        pytest.skip("pET-28a_plus.dna not available (set P2D_VECTOR_SAMPLES)")
    return confirm(analyze_vector(read_vector_file(path)))


def real(accession):
    path = SAMPLES / "constructs" / f"{accession}.gb"
    if not path.is_file():
        pytest.skip(f"{accession}.gb not available")
    return cg.read_construct(str(path))


@pytest.mark.parametrize("accession", sorted(MATCHING))
def test_p2d_rebuilds_protein_of_each_published_construct(pet28a, host, accession):
    up, down, residues = MATCHING[accession]
    name, seq = real(accession)
    r = cg.compare(name, pet28a, seq, host)
    assert r.status == "match", (r.reason, r.real_protein, r.p2d_protein)
    assert (r.up, r.down, len(r.real_protein)) == (up, down, residues)
    assert r.real_protein == r.p2d_protein and r.coverage >= cg.MIN_COVERAGE


def test_touching_pair_real_lab_used_is_rebuilt_too(pet28a, host):
    """BamHI and EcoRI sit side by side in pET-28a(+); p2d used to refuse that pair."""
    name, seq = real("MK847901.1")
    r = cg.compare(name, pet28a, seq, host)
    assert r.status == "match" and (r.up, r.down) == ("BamHI:5166", "EcoRI:5172")
    assert len(r.real_protein) == 304 and r.real_protein == r.p2d_protein


def test_records_that_cannot_test_p2d_say_why(pet28a, host):
    for accession, fragment in [
        ("PP098726.1", "no cloning site"),  # the gene sits on the start codon, cloned seamlessly
        ("MK816391.1", "not read through"),  # a gene out of frame with the His-tag start codon
    ]:
        name, seq = real(accession)
        r = cg.compare(name, pet28a, seq, host)
        assert r.status == "not comparable" and fragment in r.reason, accession


def test_published_constructs_use_several_site_pairs_and_genes(pet28a, host):
    pairs = {(u, d) for u, d, _ in MATCHING.values()}
    assert len(pairs) >= 5 and len({n for _, _, n in MATCHING.values()}) >= 7


def test_real_vector_lists_sites_constructs_used(pet28a):
    ids = {c.id for c in list_cut_sites(pet28a)}
    for up, down, _ in MATCHING.values():
        assert up in ids and down in ids
