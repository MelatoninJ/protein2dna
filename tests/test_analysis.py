"""Reading custom plasmids and working out their start codon, cloning region and enzymes.

Real vectors cannot live in the repository (their licences differ), so tests that
need them look in the folder named by ``P2D_VECTOR_SAMPLES`` and skip when it is not
set.  Everything else runs on synthetic plasmids, including a hand-built SnapGene
``.dna`` file whose expression cassette lies on the reverse strand, like the real
pET-28a(+).
"""

from __future__ import annotations

import os
from importlib import resources
from pathlib import Path

import pytest
from Bio.Seq import Seq

import p2d
from p2d import vectorio
from p2d.analysis import analyze_vector, confirm
from p2d.enzymes import find_cut_site, list_cut_sites, plan_cut
from synth import mirror, ptest1_features, snapgene_bytes

SAMPLES = Path(os.environ["P2D_VECTOR_SAMPLES"]) if os.environ.get("P2D_VECTOR_SAMPLES") else None
PAYLOAD = "SKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLKFICTTGKLPVPWPTL"


def sample(name: str) -> Path:
    if SAMPLES is None or not (SAMPLES / name).is_file():
        pytest.skip(f"{name} not available (set P2D_VECTOR_SAMPLES)")
    return SAMPLES / name


@pytest.fixture(scope="module")
def ptest1_genbank() -> str:
    return resources.files("p2d.data").joinpath("pTEST1.gb").read_text()


@pytest.fixture(scope="module")
def host():
    return p2d.load_host("ecoli_bl21")


@pytest.fixture(scope="module")
def forward_dna(ptest1_genbank):
    v = p2d.load_bundled()
    return snapgene_bytes(v.seq, ptest1_features(v))


@pytest.fixture(scope="module")
def reverse_dna():
    """pTEST1 stored the way SnapGene stores pET-28a(+): cassette on the minus strand."""
    v = p2d.load_bundled()
    rc = str(Seq(v.seq).reverse_complement())
    return snapgene_bytes(rc, mirror(ptest1_features(v), len(v.seq)))


# --------------------------------------------------------------------- reading


def test_genbank_text_bytes_and_path_agree(ptest1_genbank, tmp_path):
    path = tmp_path / "x.gb"
    path.write_text(ptest1_genbank)
    a = vectorio.read_vector_text(ptest1_genbank)
    b = vectorio.read_vector_bytes(ptest1_genbank.encode(), "whatever.gb")
    c = vectorio.read_vector_file(path)
    assert a.seq == b.seq == c.seq and len(a) == 1149 and a.circular


def test_formats_are_sniffed_without_help(ptest1_genbank):
    seq = p2d.load_bundled().seq
    assert vectorio.read_vector_bytes(ptest1_genbank.encode()).record.features  # LOCUS
    fasta = vectorio.read_vector_text(f">my plasmid\n{seq[:60]}\n{seq[60:]}\n")
    assert fasta.seq == seq and any("Sequence only" in n for n in fasta.notes)
    raw = vectorio.read_vector_text("\n".join(seq[i : i + 60] for i in range(0, len(seq), 60)))
    assert raw.seq == seq and raw.name == "pasted_vector"


def test_snapgene_file_is_read_with_features(forward_dna):
    v = vectorio.read_vector_bytes(forward_dna, "pTEST1.dna")
    labels = {f.qualifiers["label"][0] for f in v.record.features}
    assert {"T7 promoter", "RBS", "ATG", "MCS", "KanR"} <= labels
    assert v.circular and len(v) == 1149


@pytest.mark.parametrize(
    "data, filename, message",
    [
        (b"", "x.gb", "empty"),
        (b"   \n\n", "x.fa", "empty"),
        (b">a\nACGT\n", "x.fa", "bp"),
        (b">a\n" + b"ACGT" * 200 + b"\n>b\n" + b"ACGT" * 200, "x.fa", "2 records"),
        (b"\x00\x01\x02garbage" * 50, "x.bin", "binary"),
        (b">a\n" + b"ACGT" * 150 + b"XQZ" * 20 + b"\n", "x.fa", "Unsupported character"),
        (
            b">p\nMKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALPDAQFEVV\n",
            "p.fa",
            "protein",
        ),
        (b"just some prose that is not DNA at all", "x.txt", "format"),
        (b"\x09\x00\x00\x00\x0eSnapGene" + b"\xff" * 40, "x.dna", "SnapGene"),
        (b">a\n" + b"ACGTN" * 200, "x.fa", "ambiguous"),
    ],
)
def test_unusable_input_is_a_clean_error(data, filename, message):
    with pytest.raises(vectorio.VectorFileError, match=message):
        vectorio.read_vector_bytes(data, filename)


def test_size_limit_is_enforced(tmp_path):
    big = tmp_path / "big.fa"
    big.write_bytes(b">x\n" + b"A" * (vectorio.MAX_BYTES + 1))
    with pytest.raises(vectorio.VectorFileError, match="larger"):
        vectorio.read_vector_file(big)
    with pytest.raises(vectorio.VectorFileError, match="No such file"):
        vectorio.read_vector_file(tmp_path / "missing.gb")


def test_feature_text_strips_html_from_snapgene_notes(forward_dna):
    v = vectorio.read_vector_bytes(forward_dna, "p.dna")
    atg = next(f for f in v.record.features if f.qualifiers["label"][0] == "ATG")
    assert "<" not in vectorio.feature_text(atg) and "start codon" in vectorio.feature_text(atg)
    assert vectorio.is_start_codon_feature(atg, v.seq)


# ----------------------------------------------------------------- orientation


def test_forward_file_is_left_alone(forward_dna):
    v = vectorio.read_vector_bytes(forward_dna, "p.dna")
    same, flipped, notes = vectorio.normalize_orientation(v)
    assert not flipped and same is v and any("start codon" in n for n in notes)


def test_reverse_cassette_is_flipped_back_exactly(forward_dna, reverse_dna):
    fwd = vectorio.read_vector_bytes(forward_dna, "f.dna")
    rev = vectorio.read_vector_bytes(reverse_dna, "r.dna")
    assert rev.seq != fwd.seq
    oriented, flipped, _ = vectorio.normalize_orientation(rev)
    assert flipped and oriented.seq == fwd.seq

    def spans(v):
        return sorted(
            (
                f.qualifiers["label"][0],
                int(f.location.start),
                int(f.location.end),
                f.location.strand,
            )
            for f in v.record.features
        )

    assert spans(oriented) == spans(fwd), "features must land on the same bases after flipping"


def test_an_antibiotic_marker_cannot_outvote_the_expression_cassette(reverse_dna):
    """KanR points the other way in the file; only expression features may decide."""
    rev = vectorio.read_vector_bytes(reverse_dna, "r.dna")
    plus, minus, _ = vectorio.orientation_votes(rev)
    assert minus > plus


def test_promoter_sequence_alone_orients_a_bare_reverse_sequence():
    v = p2d.load_bundled()
    rc = str(Seq(v.seq).reverse_complement())
    bare = vectorio.read_vector_text(rc, "bare")
    oriented, flipped, notes = vectorio.normalize_orientation(bare)
    assert flipped and oriented.seq == v.seq and any("T7 promoter" in n for n in notes)


def test_no_signals_means_no_flip_and_a_note():
    seq = "ACGT" * 150
    v = vectorio.read_vector_text(seq, "plain")
    same, flipped, notes = vectorio.normalize_orientation(v)
    assert not flipped and "kept as given" in notes[0]


# -------------------------------------------------------------------- analysis


def test_reverse_stored_plasmid_analyses_like_the_forward_one(forward_dna, reverse_dna):
    a = analyze_vector(vectorio.read_vector_bytes(forward_dna, "f.dna"))
    b = analyze_vector(vectorio.read_vector_bytes(reverse_dna, "r.dna"))
    assert not a.flipped and b.flipped
    assert a.best_start.position == b.best_start.position == 87
    assert (a.region.start, a.region.end) == (b.region.start, b.region.end) == (144, 180)
    assert a.best_start.confidence == "high" and a.ready and b.ready
    assert a.best_start.n_terminal == b.best_start.n_terminal
    sites = lambda v: [(c.enzyme, c.start) for c in list_cut_sites(v)]  # noqa: E731
    assert sites(a.vector) == sites(b.vector), "the same sites at the same positions"


def test_annotated_start_is_checked_against_the_sequence(forward_dna):
    """An annotation that points at a base that is not ATG must not become a start codon."""
    v = vectorio.read_vector_bytes(forward_dna, "f.dna")
    wrong = snapgene_bytes(
        v.seq,
        [
            {
                "name": "ATG",
                "type": "CDS",
                "start": 89,
                "end": 91,
                "notes": {"product": "start codon"},
            }
        ],
    )
    a = analyze_vector(vectorio.read_vector_bytes(wrong, "w.dna"))
    assert all(s.position != 88 for s in a.starts)


def test_the_mcs_is_only_reported_when_annotated_never_guessed(ptest1_genbank):
    custom = analyze_vector(vectorio.read_vector_text(ptest1_genbank, "custom"))  # no MCS feature
    assert custom.region is None
    assert custom.best_start.position == 87 and custom.best_start.confidence == "high"
    assert custom.ready and custom.best_start.n_terminal.startswith("MGSSHHHHHHSSGLVPRGS")
    bare = analyze_vector(vectorio.read_vector_text(p2d.load_bundled().seq, "bare"))
    assert bare.region is None


def test_bare_sequence_still_gets_a_start_but_not_a_confident_one():
    a = analyze_vector(vectorio.read_vector_text(p2d.load_bundled().seq, "bare"))
    assert a.best_start.position == 87 and a.best_start.confidence == "medium"
    assert not a.ready and any("not certain" in w for w in a.warnings)


def test_nothing_to_go_on_gives_warnings_not_a_crash():
    a = analyze_vector(vectorio.read_vector_text("ACGTTGCA" * 100, "noise"))
    assert a.starts == [] and a.region is None and not a.ready
    assert any("start codon" in w for w in a.warnings)


# --------------------------------------------------------------------- confirm


def test_confirm_uses_the_best_start_only_when_it_is_certain(forward_dna):
    a = analyze_vector(vectorio.read_vector_bytes(forward_dna, "f.dna"))
    assert confirm(a).expression_start == 87
    guess = analyze_vector(vectorio.read_vector_text(p2d.load_bundled().seq, "bare"))
    with pytest.raises(ValueError, match="not certain"):
        confirm(guess)
    assert confirm(guess, start=87).expression_start == 87


@pytest.mark.parametrize("start", [5, -1, 10**6])
def test_confirm_rejects_a_start_that_is_not_an_atg(forward_dna, start):
    a = analyze_vector(vectorio.read_vector_bytes(forward_dna, "f.dna"))
    with pytest.raises(ValueError, match="not an ATG"):
        confirm(a, start=start)


def test_a_confirmed_custom_vector_designs_end_to_end(reverse_dna, host):
    """Upload (reverse-stored) -> analyse -> confirm -> pick two sites -> design -> verify."""
    v = confirm(analyze_vector(vectorio.read_vector_bytes(reverse_dna, "r.dna")))
    site = plan_cut(v, find_cut_site(v, "NdeI:144"), find_cut_site(v, "XhoI:174"))
    plan = p2d.ConstructPlan.single_chain(PAYLOAD)
    a = p2d.assemble(p2d.design_insert(plan, v, site, host, seed=1))
    assert a.ok, a.errors
    real = str(Seq(a.seq[v.expression_start :][: 3 * len(a.protein)]).translate())
    assert real == a.protein and a.protein.startswith("MGSSHHHHHHSSGLVPRGSHM" + PAYLOAD)


# --------------------------------------------- real vectors (when samples exist)


def test_real_pet28a_snapgene_file(host):
    """SnapGene's own pET-28a(+): cassette on the minus strand, start codon annotated."""
    a = analyze_vector(vectorio.read_vector_file(sample("pET-28a_plus.dna")))
    assert a.flipped, "SnapGene stores the pET-28a(+) cassette on the minus strand"
    assert a.best_start.confidence == "high" and a.ready
    assert a.best_start.n_terminal.startswith("MGSSHHHHHHSSGLVPRGSHMASMTGGQQMG")
    assert a.region.source == "annotation: MCS"
    # MCS annotations are loose (NdeI lies just outside this one), so sites are listed
    # across the whole stretch after the start codon and the user chooses
    window = (a.best_start.position, a.region.end + 20)
    names = {c.enzyme for c in list_cut_sites(a.vector, window)}
    assert {"NdeI", "NheI", "BamHI", "EcoRI", "SacI", "SalI", "HindIII", "NotI", "XhoI"} <= names


# the protease site each GST vector carries between GST and the MCS
PGEX_CLEAVAGE = {"U13853.gb": "LVPRGS", "U13852.gb": "IEGRGI", "U13850.gb": "LVPRGS"}


@pytest.mark.parametrize("name", sorted(PGEX_CLEAVAGE))
def test_real_pgex_backbones(name, host):
    a = analyze_vector(vectorio.read_vector_file(sample(name)))
    assert not a.flipped and a.ready
    assert a.best_start.position == 257 and a.best_start.confidence == "high"
    assert a.region.source == "annotation: MCS"
    v = confirm(a)
    inside = {c.enzyme for c in list_cut_sites(v, (a.region.start - 6, a.region.end + 6))}
    assert {"BamHI", "EcoRI"} <= inside
    plan = p2d.ConstructPlan.single_chain(PAYLOAD)
    bam = next(c for c in list_cut_sites(v) if c.enzyme == "BamHI" and c.count == 1)
    eco = next(c for c in list_cut_sites(v) if c.enzyme == "EcoRI" and c.count == 1)
    result = p2d.assemble(p2d.design_insert(plan, v, plan_cut(v, bam, eco), host))
    assert result.ok, result.errors
    real = str(Seq(result.seq[257:][: 3 * len(result.protein)]).translate())
    assert real == result.protein, "reported protein must match translation from the real start"
    assert result.protein.startswith(a.best_start.n_terminal), "GST stays on the N-terminus"
    assert PGEX_CLEAVAGE[name] in result.protein and PAYLOAD in result.protein
