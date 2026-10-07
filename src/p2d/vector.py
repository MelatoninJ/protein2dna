"""Vectors and insertion sites.

Design rule: **never trust the annotation labels.**  GenBank records in the
wild (Addgene included) routinely carry mislabelled or missing features.  The
ORF and the restriction sites are re-derived from the sequence itself, and a
mismatch against the record's own annotation is reported as a warning rather
than silently accepted.
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field
from enum import Enum
from pathlib import Path

from Bio import SeqIO
from Bio.Restriction import RestrictionBatch
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

IUPAC = {
    "A": "A",
    "C": "C",
    "G": "G",
    "T": "T",
    "R": "[AG]",
    "Y": "[CT]",
    "S": "[GC]",
    "W": "[AT]",
    "K": "[GT]",
    "M": "[AC]",
    "B": "[CGT]",
    "D": "[AGT]",
    "H": "[ACT]",
    "V": "[ACG]",
    "N": "[ACGT]",
}

STOP_CODONS = {"TAA", "TAG", "TGA"}


class StartSource(str, Enum):
    """Where the final construct's start codon comes from."""

    VECTOR = "vector"  # vector supplies ATG (and usually an N-terminal tag)
    UPSTREAM_SITE = "upstream_site"  # e.g. NdeI CAT-ATG donates the ATG itself


@dataclass(frozen=True)
class EnzymeHit:
    name: str
    site: str
    start: int  # 0-based, first base of the recognition sequence
    end: int  # exclusive

    def __len__(self) -> int:
        return self.end - self.start


@dataclass(frozen=True)
class InsertionSite:
    """A cloning strategy for one vector: which two enzymes, and where the
    reading frame begins.

    ``start_offset_in_site`` is the index of the A of ATG inside the upstream
    recognition sequence -- 2 for NdeI (CAT|ATG), 1 for NcoI (C|CATGG).
    """

    name: str
    upstream_enzyme: str
    downstream_enzyme: str
    start_source: StartSource = StartSource.VECTOR
    vector_orf_start: int | None = None
    start_offset_in_site: int = 0
    note: str = ""


def site_regex(site: str) -> re.Pattern:
    return re.compile("".join(IUPAC[b] for b in site.upper()))


@dataclass
class Vector:
    record: SeqRecord
    name: str = ""
    circular: bool = True
    sites: dict[str, InsertionSite] = field(default_factory=dict)

    def __post_init__(self) -> None:
        if not self.name:
            self.name = self.record.name or self.record.id or "unnamed_vector"
        self.seq = str(self.record.seq).upper()

    # ------------------------------------------------------------------ io
    @classmethod
    def from_genbank(cls, path: str | Path, **kw) -> Vector:
        record = SeqIO.read(str(path), "genbank")
        return cls(record=record, **kw)

    @classmethod
    def from_sequence(cls, seq: str, *, name: str = "vector", **kw) -> Vector:
        rec = SeqRecord(Seq(seq.upper()), id=name, name=name, description="")
        return cls(record=rec, name=name, **kw)

    def __len__(self) -> int:
        return len(self.seq)

    # --------------------------------------------------------- enzyme search
    def find_enzyme(self, enzyme: str, *, require_unique: bool = True) -> EnzymeHit:
        hits = self.find_all(enzyme)
        if not hits:
            raise ValueError(f"{enzyme} does not cut {self.name}")
        if require_unique and len(hits) > 1:
            positions = ", ".join(str(h.start) for h in hits)
            raise ValueError(
                f"{enzyme} is not unique in {self.name}: {len(hits)} sites at {positions}"
            )
        return hits[0]

    def find_all(self, enzyme: str) -> list[EnzymeHit]:
        enz = _enzyme(enzyme)
        pattern = site_regex(enz.site)
        search_space = self.seq + (self.seq[: len(enz.site) - 1] if self.circular else "")
        out = []
        for m in pattern.finditer(search_space):
            if m.start() >= len(self.seq):
                continue
            out.append(EnzymeHit(enzyme, enz.site, m.start(), m.start() + len(enz.site)))
        return out

    def cuts_once(self, enzyme: str) -> bool:
        return len(self.find_all(enzyme)) == 1

    # ------------------------------------------------------------ orf lookup
    def detect_orfs(self, *, min_aa: int = 20) -> list[tuple[int, int]]:
        """Re-derive ATG..stop ORFs on the plus strand from the sequence.

        Returns ``(start, stop_end)`` pairs, 0-based, longest first.
        """
        found: list[tuple[int, int]] = []
        for start in range(len(self.seq) - 2):
            if self.seq[start : start + 3] != "ATG":
                continue
            for i in range(start, len(self.seq) - 2, 3):
                if self.seq[i : i + 3] in STOP_CODONS:
                    if (i - start) // 3 >= min_aa:
                        found.append((start, i + 3))
                    break
        found.sort(key=lambda p: p[1] - p[0], reverse=True)
        return found

    def annotated_cds(self) -> list[tuple[int, int, str]]:
        out = []
        for f in self.record.features:
            if f.type == "CDS":
                label = (f.qualifiers.get("label") or f.qualifiers.get("gene") or ["CDS"])[0]
                out.append((int(f.location.start), int(f.location.end), label))
        return out

    def feature_label_at(self, pos: int) -> str | None:
        for f in self.record.features:
            if int(f.location.start) <= pos < int(f.location.end):
                lab = f.qualifiers.get("label") or f.qualifiers.get("note")
                if lab:
                    return lab[0]
        return None


def _enzyme(name: str):
    batch = RestrictionBatch([name])
    return list(batch)[0]


def enzyme_site(name: str) -> str:
    return _enzyme(name).site.upper()
