"""Read a plasmid from whatever the user has, and put it in a standard orientation.

Accepted inputs, best first: a GenBank file or a SnapGene ``.dna`` file (both carry
the annotated features that make the analysis reliable), a FASTA file or pasted bare
sequence (features absent, so the user has to confirm more).  Everything is parsed
with Biopython and never executed; size, record count and alphabet are checked
before a sequence is accepted.

Orientation matters because files store a plasmid in whichever direction its source
used.  SnapGene's own pET-28a(+) keeps the T7 cassette on the reverse strand while
the antibiotic marker runs forward.  Everything downstream reads the plus strand, so
a plasmid whose expression cassette lies on the minus strand is reverse-complemented
as a whole (features included) before analysis.
"""

from __future__ import annotations

import io
import os
import re
from pathlib import Path

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from .vector import Vector

MAX_BYTES = 5_000_000
MIN_BP, MAX_BP = 500, 250_000
_IUPAC = set("ACGTRYSWKMBDHVN")

GENBANK_SUFFIXES = {".gb", ".gbk", ".genbank", ".gbff"}
FASTA_SUFFIXES = {".fa", ".fasta", ".fna", ".fas", ".seq", ".txt"}


class VectorFileError(ValueError):
    """The input cannot be used as a plasmid; the message is safe to show the user."""


# ------------------------------------------------------------------- reading


def read_vector_file(path: str | Path) -> Vector:
    path = Path(path)
    if not path.is_file():
        raise VectorFileError(f"No such file: {path.name}")
    if path.stat().st_size > MAX_BYTES:
        raise VectorFileError(f"{path.name} is larger than {MAX_BYTES // 1_000_000} MB.")
    return read_vector_bytes(path.read_bytes(), path.name)


def read_vector_text(text: str, name: str = "pasted_vector") -> Vector:
    return read_vector_bytes(text.encode("utf-8", errors="replace"), name)


def read_vector_bytes(data: bytes, filename: str = "vector") -> Vector:
    """Parse ``data`` as GenBank, SnapGene, FASTA or a bare sequence (sniffed)."""
    if not data or not data.strip():
        raise VectorFileError("The file is empty.")
    if len(data) > MAX_BYTES:
        raise VectorFileError(f"The input is larger than {MAX_BYTES // 1_000_000} MB.")
    fmt = _sniff(data, Path(filename).suffix.lower())
    stem = Path(filename).stem or "vector"
    try:
        record = _parse(data, fmt)
    except VectorFileError:
        raise
    except Exception as exc:  # Biopython raises many types on malformed input
        raise VectorFileError(
            f"Could not read this as a {_FORMAT_NAMES[fmt]} file ({type(exc).__name__}). "
            "Check that it holds one plasmid."
        ) from exc

    seq = str(record.seq).upper().replace("U", "T")
    _check_sequence(seq)
    record.seq = Seq(seq)
    topology = record.annotations.get("topology")
    notes: list[str] = []
    if topology is None:
        notes.append(
            f"No topology in the {_FORMAT_NAMES[fmt]} input; treated as a circular plasmid."
        )
    if fmt in {"fasta", "raw"}:
        notes.append(
            "Sequence only: no features, so the start codon and cloning site need confirming."
        )
    name = _clean_name(record.name if record.name not in ("", "<unknown name>") else stem) or stem
    v = Vector(record=record, name=name, circular=topology != "linear")
    v.notes = notes
    return v


_FORMAT_NAMES = {"genbank": "GenBank", "snapgene": "SnapGene", "fasta": "FASTA", "raw": "sequence"}


def _sniff(data: bytes, suffix: str) -> str:
    if suffix == ".dna" or data[:1] == b"\x09" and b"SnapGene" in data[:16]:
        return "snapgene"
    head = data[:2000].lstrip()
    if suffix in GENBANK_SUFFIXES or head.upper().startswith(b"LOCUS"):
        return "genbank"
    if head.startswith(b">") or suffix in {".fa", ".fasta", ".fna", ".fas"}:
        return "fasta"
    text = data[:5000].decode("utf-8", errors="ignore")
    if re.fullmatch(r"[\sACGTUNacgtun0-9]+", text):
        return "raw"
    if b"\x00" in data[:1000]:
        raise VectorFileError(
            "This looks like a binary file p2d does not know. Use GenBank, .dna or FASTA."
        )
    raise VectorFileError(
        "Could not tell the file format. Use GenBank, SnapGene .dna, FASTA, or a bare sequence."
    )


def _parse(data: bytes, fmt: str) -> SeqRecord:
    if fmt == "snapgene":
        return SeqIO.read(io.BytesIO(data), "snapgene")
    text = data.decode("utf-8", errors="replace")
    if fmt == "raw":
        seq = re.sub(r"[\s\d]+", "", text)
        return SeqRecord(Seq(seq), id="<unknown id>", name="<unknown name>")
    records = list(SeqIO.parse(io.StringIO(text), fmt))
    if len(records) != 1:
        raise VectorFileError(
            f"The file holds {len(records)} records; give p2d one plasmid at a time."
        )
    return records[0]


def _check_sequence(seq: str) -> None:
    bad = sorted(set(seq) - _IUPAC)
    if bad:
        if set(seq) <= set("ACDEFGHIKLMNPQRSTVWY") and len(set(seq)) > 8:
            raise VectorFileError(
                "This looks like a protein sequence. Put proteins in the design box; "
                "this input takes a plasmid as DNA."
            )
        raise VectorFileError(f"Unsupported character(s) in the sequence: {', '.join(bad[:5])}.")
    if not MIN_BP <= len(seq) <= MAX_BP:
        raise VectorFileError(
            f"The sequence is {len(seq)} bp; plasmids of {MIN_BP} to {MAX_BP} bp are supported."
        )
    if sum(1 for c in seq if c not in "ACGT") > 0.01 * len(seq):
        raise VectorFileError("More than 1% of the sequence is ambiguous (N or IUPAC codes).")


def _clean_name(name: str) -> str:
    return re.sub(r"[^\w.+()\-]+", "_", name.strip())[:60]


# ---------------------------------------------------------------- orientation

_HTML = re.compile(r"<[^>]+>")
T7_PROMOTER = "TAATACGACTCACTATAG"
_EXPRESSION_PROMOTER = re.compile(r"\b(t7|t5|tac|trc|lacuv5|arabad|pbad)\b", re.I)


def feature_text(feature) -> str:
    """Searchable lower-case text of a feature; SnapGene wraps notes in HTML tags."""
    parts = []
    for key in ("label", "name", "product", "gene", "standard_name", "note", "bound_moiety"):
        for value in feature.qualifiers.get(key, []):
            parts.append(_HTML.sub(" ", str(value)))
    return " ".join(parts).lower()


def is_start_codon_feature(feature, seq: str) -> bool:
    """A 3 bp ATG feature labelled as a start codon, read on its own strand."""
    start, end = int(feature.location.start), int(feature.location.end)
    if end - start != 3:
        return False
    codon = seq[start:end]
    if feature.location.strand == -1:
        codon = str(Seq(codon).reverse_complement())
    if codon != "ATG":
        return False
    text = feature_text(feature)
    return (
        "start codon" in text or feature.qualifiers.get("label", [""])[0].strip().upper() == "ATG"
    )


def orientation_votes(vector: Vector) -> tuple[int, int, list[str]]:
    """Evidence that the expression cassette is on the plus (first) or minus (second) strand."""
    seq = vector.seq
    plus = minus = 0
    reasons: list[str] = []

    def vote(strand, weight, why):
        nonlocal plus, minus
        if strand in (1, None, 0):
            plus += weight
        else:
            minus += weight
        reasons.append(f"{why} on the {'plus' if strand in (1, None, 0) else 'minus'} strand")

    # Only features whose direction a file states outright may vote: SnapGene leaves it
    # off RBS, MCS and operators, and Biopython then reports them as forward.
    for f in vector.record.features:
        text = feature_text(f)
        if is_start_codon_feature(f, seq):
            vote(f.location.strand, 3, "annotated start codon")
        elif (
            ("promoter" in text or f.type in {"promoter", "regulatory"})
            and "terminator" not in text
            and _EXPRESSION_PROMOTER.search(text)
        ):
            vote(f.location.strand, 3, "expression promoter")

    fwd, rev = seq.count(T7_PROMOTER), seq.count(str(Seq(T7_PROMOTER).reverse_complement()))
    if fwd:
        plus += 2 * fwd
        reasons.append(f"T7 promoter sequence on the plus strand ({fwd}x)")
    if rev:
        minus += 2 * rev
        reasons.append(f"T7 promoter sequence on the minus strand ({rev}x)")
    return plus, minus, reasons


def normalize_orientation(vector: Vector) -> tuple[Vector, bool, list[str]]:
    """Return ``(vector, flipped, notes)`` with the expression cassette on the plus strand."""
    plus, minus, reasons = orientation_votes(vector)
    if minus <= plus:
        notes = (
            reasons
            if plus
            else ["No expression promoter, RBS or start codon to orient by; kept as given."]
        )
        return vector, False, notes
    record = vector.record.reverse_complement(
        id=True, name=True, description=True, features=True, annotations=True
    )
    flipped = Vector(record=record, name=vector.name, circular=vector.circular)
    flipped.notes = list(getattr(vector, "notes", [])) + [
        "Reverse-complemented so the expression cassette reads forward."
    ]
    return flipped, True, reasons


# -------------------------------------------------------------------- library


def library_dir() -> Path:
    """Where the user keeps vectors that must not leave their computer.

    Files here (SnapGene maps under their own licence, lab plasmids) appear in the vector
    list without being uploaded each time or copied into the project.
    """
    return Path(os.environ.get("P2D_VECTOR_DIR") or Path.home() / ".p2d" / "vectors")


def list_library() -> list[str]:
    folder = library_dir()
    if not folder.is_dir():
        return []
    suffixes = GENBANK_SUFFIXES | FASTA_SUFFIXES | {".dna"}
    return sorted(p.name for p in folder.iterdir() if p.is_file() and p.suffix.lower() in suffixes)


def read_library_vector(filename: str) -> Vector:
    """Read ``filename`` from the library; only names the listing offers are accepted."""
    if filename not in list_library():
        raise VectorFileError("That file is not in your vector folder.")
    return read_vector_file(library_dir() / filename)
