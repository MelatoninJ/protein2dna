"""Work out, from a plasmid alone, where translation starts.

The start codon is the one thing p2d cannot do without, and the one thing a file does
not always say, so it is answered with its evidence for a person to confirm.  An
annotated start-codon feature is the best evidence; an annotated CDS that runs through
an annotated cloning site is next; a ribosome binding site followed by an ATG whose ORF
reaches the cloning site, or sits behind a T7 promoter, is the fallback.

Where to cut is deliberately *not* guessed.  Restriction sites are listed in full by
``p2d.enzymes.list_cut_sites`` and the user chooses; an annotated MCS is reported here
only as a landmark.  Annotations are evidence, not truth (a core rule of this package),
so every annotated answer is checked against the sequence itself.
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field

from Bio.Seq import Seq

from .vector import STOP_CODONS, Vector
from .vectorio import (
    T7_PROMOTER,
    feature_text,
    is_start_codon_feature,
    normalize_orientation,
)

_MCS = re.compile(r"multiple cloning site|\bmcs\b|polylinker|poly-linker", re.I)
_NOT_EXPRESSED = re.compile(
    r"resistance|\b(kan|amp|cam|tet|bla|cat|lac[iz]|rop|rep[a-z]?|ori)\w*\b|replication", re.I
)
_SD = ("AGGAGG", "AGGAG", "GGAGG", "AGGA", "GGAG", "GAGG")

N_TERMINAL_PREVIEW = 30
CONFIDENCE_RANK = {"high": 0, "medium": 1, "low": 2}


@dataclass(frozen=True)
class StartCandidate:
    position: int  # 0-based index of the A of ATG on the plus strand
    confidence: str  # "high" | "medium" | "low"
    evidence: tuple[str, ...]
    orf_end: int  # end of the ORF (stop codon included), 0-based exclusive
    n_terminal: str  # residues encoded between the start codon and the cloning region


@dataclass(frozen=True)
class Region:
    start: int
    end: int
    source: str  # "annotation: MCS" or "guessed from enzyme sites"
    confident: bool


@dataclass
class VectorAnalysis:
    vector: Vector  # in the standard orientation
    flipped: bool
    orientation_notes: list[str]
    starts: list[StartCandidate]
    region: Region | None  # the annotated MCS, a landmark only
    warnings: list[str] = field(default_factory=list)

    @property
    def best_start(self) -> StartCandidate | None:
        return self.starts[0] if self.starts else None

    @property
    def ready(self) -> bool:
        """True when the start codon needs no confirmation."""
        start = self.best_start
        return bool(start and start.confidence == "high")


# --------------------------------------------------------------------- region


def find_region(vector: Vector) -> Region | None:
    """The annotated MCS, if the file has one.  Never guessed."""
    for f in vector.record.features:
        if f.type != "source" and _MCS.search(feature_text(f)):
            return Region(int(f.location.start), int(f.location.end), "annotation: MCS", True)
    return None


# ---------------------------------------------------------------------- starts


def _orf_end(seq: str, start: int) -> int | None:
    for i in range(start, len(seq) - 2, 3):
        if seq[i : i + 3] in STOP_CODONS:
            return i + 3
    return None


def _sd_score(seq: str, atg: int) -> int:
    window = seq[max(0, atg - 16) : max(0, atg - 3)]
    return max((len(m) for m in _SD if m in window), default=0)


def find_starts(vector: Vector, region: Region | None) -> list[StartCandidate]:
    seq = vector.seq
    found: dict[int, dict] = {}

    def add(pos, why, strength):
        entry = found.setdefault(pos, {"why": [], "strength": 0})
        entry["why"].append(why)
        entry["strength"] += strength

    for f in vector.record.features:
        if f.location.strand == -1:
            continue
        s, e = int(f.location.start), int(f.location.end)
        if is_start_codon_feature(f, seq):
            add(s, f"annotated as the start codon ({feature_label(f)})", 3)
        elif f.type == "CDS" and e - s >= 90 and seq[s : s + 3] == "ATG":
            text = feature_text(f)
            if _NOT_EXPRESSED.search(text):
                continue
            if region and s < region.start and (stop := _orf_end(seq, s)) and stop > region.start:
                add(
                    s,
                    f"annotated CDS '{feature_label(f)}' reads in frame into the cloning region",
                    2,
                )

    if region:
        for pos in range(len(seq) - 2):
            if seq[pos : pos + 3] != "ATG" or pos >= region.start:
                continue
            sd = _sd_score(seq, pos)
            stop = _orf_end(seq, pos)
            if sd >= 4 and stop and stop > region.start:
                add(
                    pos,
                    f"ribosome-binding-site-like motif ({sd} nt) upstream, ORF reaches the MCS",
                    1 + (sd >= 5),
                )

    if not region:
        for hit in re.finditer(T7_PROMOTER, seq):
            for pos in range(hit.end() + 8, min(hit.end() + 150, len(seq) - 2)):
                stop = _orf_end(seq, pos) if seq[pos : pos + 3] == "ATG" else None
                sd = _sd_score(seq, pos) if stop else 0
                if sd >= 4 and stop and stop - pos >= 90:
                    add(
                        pos,
                        f"RBS-like motif ({sd} nt) then ATG, downstream of a T7 promoter",
                        1 + (sd >= 5),
                    )

    out = []
    for pos, entry in found.items():
        stop = _orf_end(seq, pos)
        if stop is None:
            continue
        agree = len(entry["why"]) >= 2
        strength = entry["strength"]
        conf = (
            "high"
            if strength >= 3 and (agree or region is None or stop > region.start)
            else ("medium" if strength >= 2 else "low")
        )
        # residues the vector adds in front of the cloning site; if the file does not say
        # where that is, the first stretch is still enough to recognise a tag
        reach = (region.start - pos) // 3 if region and region.start > pos else N_TERMINAL_PREVIEW
        n = min(reach, (stop - pos) // 3 - 1, 60)
        n_term = str(Seq(seq[pos : pos + 3 * n]).translate())
        out.append(StartCandidate(pos, conf, tuple(entry["why"]), stop, n_term))
    out.sort(key=lambda c: (CONFIDENCE_RANK[c.confidence], -len(c.evidence), c.position))
    return out


def feature_label(feature) -> str:
    for key in ("label", "name", "product", "gene"):
        for value in feature.qualifiers.get(key, []):
            text = re.sub(r"<[^>]+>", " ", str(value)).strip()
            if text:
                return text[:40]
    return feature.type


# ---------------------------------------------------------------------- driver


def analyze_vector(vector: Vector) -> VectorAnalysis:
    """Orient the vector, then find its start codon candidates."""
    vector, flipped, notes = normalize_orientation(vector)
    region = find_region(vector)
    starts = find_starts(vector, region)

    warnings: list[str] = []
    if not starts:
        warnings.append("No start codon could be found; point to one to continue.")
    elif starts[0].confidence != "high":
        warnings.append("The start codon is not certain; confirm it before designing.")
    if len(starts) > 1 and starts[0].confidence == starts[1].confidence:
        warnings.append("More than one start codon is plausible; pick the one the host uses.")
    return VectorAnalysis(
        vector=vector,
        flipped=flipped,
        orientation_notes=notes,
        starts=starts,
        region=region,
        warnings=warnings,
    )


def confirm(analysis: VectorAnalysis, *, start: int | None = None) -> Vector:
    """Fix the start codon on the (oriented) vector and return it.

    With no argument the best candidate is used, but only if it is certain.
    """
    v = analysis.vector
    if start is None:
        best = analysis.best_start
        if not best or best.confidence != "high":
            raise ValueError("The start codon is not certain; pass start= explicitly.")
        start = best.position
    if not 0 <= start <= len(v.seq) - 3 or v.seq[start : start + 3] != "ATG":
        raise ValueError(f"Position {start + 1} is not an ATG.")
    v.expression_start = start
    return v
