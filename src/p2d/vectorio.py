"""Read a plasmid from GenBank, SnapGene .dna, FASTA or bare sequence, and orient it.

GenBank and SnapGene files carry annotations; FASTA and bare sequence do not, so more has to be
confirmed. Input is parsed with Biopython, never executed, and its size, record count and
alphabet are checked first. Everything downstream reads the plus strand, so a plasmid whose
expression cassette is on the minus strand (SnapGene's pET-28a(+), for one) is
reverse-complemented whole, features included.
"""

from __future__ import annotations

import io
import os
import re
from dataclasses import dataclass
from pathlib import Path

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord

from .vector import Vector

MAX_BYTES = 5_000_000
MIN_BP, MAX_BP = 500, 250_000
_IUPAC = set("ACGTRYSWKMBDHVN")

GENBANK_SUFFIXES = {".gb", ".gbk", ".genbank", ".gbff"}
FASTA_SUFFIXES = {".fa", ".fasta", ".fna", ".fas", ".seq", ".txt"}


class VectorFileError(ValueError):
    """The input cannot be used as a plasmid; the message is safe to show the user."""


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
    header_notes = _header_features(record) if fmt == "fasta" else []
    topology = record.annotations.get("topology")
    notes: list[str] = []
    if topology is None:
        notes.append(
            f"No topology in the {_FORMAT_NAMES[fmt]} input; treated as a circular plasmid."
        )
    if header_notes:
        notes.extend(header_notes)
    elif fmt in {"fasta", "raw"}:
        notes.append(
            "Sequence only: no features, so the start codon and cloning site need confirming."
        )
    name = _clean_name(record.name if record.name not in ("", "<unknown name>") else stem) or stem
    v = Vector(record=record, name=name, circular=topology != "linear")
    v.notes = notes
    return v


_HEADER_START = re.compile(r"\bstart\s*=\s*(\d+)\b", re.I)
_HEADER_MCS = re.compile(r"\bmcs\s*=\s*(\d+)\s*(?:-|\.\.|to)\s*(\d+)\b", re.I)


def _header_features(record) -> list[str]:
    """Features the user wrote in a FASTA header: ``>pET-28a start=5307 mcs=5305-5380``.

    Positions are 1-based, on the strand as written.  ``start`` is the A of the expression ATG
    and becomes an annotated start codon; ``mcs`` is the cloning region, a landmark only.
    """
    text, seq, notes = record.description, str(record.seq), []
    if m := _HEADER_START.search(text):
        pos = int(m.group(1)) - 1
        if not 0 <= pos <= len(seq) - 3 or seq[pos : pos + 3] != "ATG":
            raise VectorFileError(
                f"start={m.group(1)} in the FASTA header is not an ATG in this sequence."
            )
        record.features.append(
            SeqFeature(
                FeatureLocation(pos, pos + 3, strand=1),
                type="misc_feature",
                qualifiers={"label": ["start codon"]},
            )
        )
        notes.append(f"Start codon at {pos + 1} taken from the FASTA header.")
    if m := _HEADER_MCS.search(text):
        a, b = sorted((int(m.group(1)), int(m.group(2))))
        if a < 1 or b > len(seq):
            raise VectorFileError(f"mcs={a}-{b} in the FASTA header lies outside the sequence.")
        record.features.append(
            SeqFeature(
                FeatureLocation(a - 1, b, strand=1),
                type="misc_feature",
                qualifiers={"label": ["MCS"]},
            )
        )
        notes.append(f"Cloning region {a}-{b} taken from the FASTA header.")
    return notes


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


_HTML = re.compile(r"<[^>]+>")
T7_PROMOTER = "TAATACGACTCACTATAG"
# the first residues of lacZ-alpha (MTMITPS), shared by the pUC and pBluescript families
LACZ_ALPHA_START = "ATGACCATGATTACGCCAAGC"
# an annotated CDS that is the N-terminal partner of an expression vector (a tag or carrier protein)
EXPRESSED_PARTNER = re.compile(
    r"glutathione|\bgst\b|maltose|\bmal[ek]\b|thioredoxin|\btrx[a-z]?\b|\bsumo\b|\bnus[a-z]?\b",
    re.I,
)
# the 17 bp lac operator core, as found in the pTrc99A, pET and pCold records
LAC_OPERATOR = "TTGTGAGCGGATAACAA"
RBS_MOTIFS = ("AGGAGG", "AGGAG", "GGAGG", "AGGA", "GGAG", "GAGG")
_RBS_SPACING = range(4, 11)  # bases between the end of the motif and the ATG
_PROMOTER_REACH = 250  # an expression promoter or operator lies this close behind the start
# the vector's own ATG is often part of an NcoI (CCATGG) or NdeI (CATATG) site
_ATG_SITES = {"CCATGG": ("NcoI", 2), "CATATG": ("NdeI", 3)}
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


@dataclass(frozen=True)
class CassetteStart:
    position: int  # the A of ATG on the strand that was scanned
    site: str  # "NcoI" / "NdeI" when the ATG is part of the site, else ""
    rbs: int  # length of the ribosome-binding-site-like motif
    promoter: str  # "T7 promoter" / "lac operator" when one lies behind it, else ""


def _rbs_length(seq: str, atg: int) -> int:
    best = 0
    for gap in _RBS_SPACING:
        window = seq[max(0, atg - gap - 6) : max(0, atg - gap)]
        best = max([best] + [len(m) for m in RBS_MOTIFS if window.endswith(m)])
    return best


def cassette_starts(seq: str) -> list[CassetteStart]:
    """Start codons of an expression cassette, recognised from the sequence alone.

    An ATG counts when a ribosome-binding-site-like motif (AGGAGG and shorter cores) lies 4-10
    bases before it and either the ATG is part of an NcoI or NdeI site, or an expression promoter
    (a T7 promoter, a lac operator) lies within 250 bases before it.  Vectors with none of this,
    such as pKK223-3 whose insert brings its own ATG, give nothing.  The reading frame is not
    checked here: the caller decides how long an open frame it needs.
    """
    found = []
    for hit in re.finditer("(?=ATG)", seq):
        atg = hit.start()
        rbs = _rbs_length(seq, atg)
        if not rbs:
            continue
        site = next(
            (
                name
                for motif, (name, offset) in _ATG_SITES.items()
                if seq[atg - offset : atg - offset + 6] == motif
            ),
            "",
        )
        behind = seq[max(0, atg - _PROMOTER_REACH) : atg]
        promoter = (
            "T7 promoter"
            if T7_PROMOTER in behind
            else "lac operator"
            if LAC_OPERATOR in behind
            else ""
        )
        if site or promoter:
            found.append(CassetteStart(atg, site, rbs, promoter))
    return found


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
        elif (
            f.type == "CDS"
            and f.location.strand in (1, -1)
            and len(f) >= 300
            and EXPRESSED_PARTNER.search(text)
        ):
            vote(f.location.strand, 4, "annotated fusion-partner CDS")

    fwd, rev = seq.count(T7_PROMOTER), seq.count(str(Seq(T7_PROMOTER).reverse_complement()))
    if fwd:
        plus += 2 * fwd
        reasons.append(f"T7 promoter sequence on the plus strand ({fwd}x)")
    if rev:
        minus += 2 * rev
        reasons.append(f"T7 promoter sequence on the minus strand ({rev}x)")
    lac_f = seq.count(LACZ_ALPHA_START)
    lac_r = seq.count(str(Seq(LACZ_ALPHA_START).reverse_complement()))
    if lac_f:
        plus += 3 * lac_f
        reasons.append("start of lacZ-alpha (MTMITPS) on the plus strand")
    if lac_r:
        minus += 3 * lac_r
        reasons.append("start of lacZ-alpha (MTMITPS) on the minus strand")
    # the sequence-only cassette is a fallback: annotated features, a T7 promoter or lacZ-alpha
    # say more, and an unannotated vector (pTrc, pBAD, pET-3) has nothing else to go by
    if not (plus or minus):
        reverse = str(Seq(seq).reverse_complement())
        for strand, scanned in (("plus", seq), ("minus", reverse)):
            for c in cassette_starts(scanned):
                if strand == "plus":
                    plus += 1
                else:
                    minus += 1
                where = f" in the {c.site} site" if c.site else ""
                reasons.append(f"RBS-like motif and start codon{where} on the {strand} strand")
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
