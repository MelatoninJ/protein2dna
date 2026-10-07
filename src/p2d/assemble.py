"""Assembly and validation.

The hard gate of this whole package lives here: translate the assembled
plasmid *from the vector's own start codon* and check that what comes out is
what was intended.  Everything else in the package is advisory; this one is
pass/fail.
"""

from __future__ import annotations

from dataclasses import dataclass, field

from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord

from .frame import (
    InsertDesign,
    Issue,
    Segment,
    Severity,
    format_fusion,
    fusion_segments,
)

STOP_CODONS = {"TAA", "TAG", "TGA"}


@dataclass
class Assembly:
    design: InsertDesign
    seq: str
    protein: str
    segments: list[Segment]
    issues: list[Issue] = field(default_factory=list)

    @property
    def errors(self) -> list[Issue]:
        return [i for i in self.issues if i.severity is Severity.ERROR]

    @property
    def ok(self) -> bool:
        return not self.errors

    @property
    def insert_dna(self) -> str:
        return self.design.insert_dna

    def fusion_report(self) -> str:
        return format_fusion(self.segments)

    def to_record(self, name: str | None = None) -> SeqRecord:
        d = self.design
        name = name or f"{d.vector.name}_{d.chain.name}"
        rec = SeqRecord(
            Seq(self.seq),
            id=name[:16],
            name=name[:16],
            description=f"{name} (assembled by p2d)",
        )
        rec.annotations["molecule_type"] = "ds-DNA"
        rec.annotations["topology"] = "circular" if d.vector.circular else "linear"

        ins_start = d.up.start
        rec.features.append(
            SeqFeature(
                FeatureLocation(ins_start, ins_start + len(d.insert_dna)),
                type="misc_feature",
                qualifiers={"label": ["insert (p2d)"]},
            )
        )
        rec.features.append(
            SeqFeature(
                FeatureLocation(d.coding_start, d.coding_start + len(d.coding)),
                type="CDS",
                qualifiers={"label": [d.chain.name], "translation": [d.chain.insert_aa]},
            )
        )
        orf_end = d.orf_start + 3 * len(self.protein) + 3
        rec.features.append(
            SeqFeature(
                FeatureLocation(d.orf_start, min(orf_end, len(self.seq))),
                type="CDS",
                qualifiers={"label": ["fusion ORF"], "translation": [self.protein]},
            )
        )
        return rec


def assemble(design: InsertDesign) -> Assembly:
    """Splice the insert into the vector and validate the result."""
    d = design
    v = d.vector.seq
    final = v[: d.up.start] + d.insert_dna + v[d.down.end :]

    issues = list(d.issues)
    protein, stop_found = _translate_orf(final, d.orf_start)
    segments = fusion_segments(d, final)

    # ---- hard gate 1: the payload survives, contiguously and unmutated ------
    payload = d.chain.payload
    if payload not in protein:
        issues.append(
            Issue(
                Severity.ERROR,
                "payload_lost",
                "the payload sequence does not appear in the translated fusion "
                "protein -- frameshift or premature stop at a junction",
            )
        )
    else:
        idx = protein.index(payload)
        if protein.count(payload) > 1:
            issues.append(
                Issue(Severity.WARNING, "payload_repeated", "payload occurs more than once")
            )
        issues.append(
            Issue(
                Severity.INFO,
                "payload_at",
                f"payload begins at residue {idx + 1} of the fusion",
            )
        )

    # ---- hard gate 2: insert block translates exactly as designed -----------
    designed = d.chain.insert_aa
    observed = str(Seq(d.coding).translate())
    if observed != designed:
        issues.append(
            Issue(
                Severity.ERROR,
                "coding_mismatch",
                "reverse-translated coding region does not back-translate to the "
                f"designed sequence ({observed[:30]}... vs {designed[:30]}...)",
            )
        )

    # ---- hard gate 3: the ORF terminates -----------------------------------
    if not stop_found:
        issues.append(
            Issue(Severity.ERROR, "runaway_orf", "no in-frame stop codon found after the insert")
        )

    # ---- hard gate 4: no stop inside the designed region --------------------
    if "*" in observed:
        issues.append(
            Issue(Severity.ERROR, "internal_stop", "the coding region contains a stop codon")
        )

    # ---- frame bookkeeping, reported explicitly -----------------------------
    issues.append(
        Issue(
            Severity.INFO,
            "frame",
            f"ORF starts at {d.orf_start}; payload block at {d.coding_start} "
            f"(offset {(d.coding_start - d.orf_start) % 3} -- must be 0); "
            f"pads: 5'={len(d.pad_up)} nt, 3'={len(d.pad_down)} nt",
        )
    )
    if (d.coding_start - d.orf_start) % 3 != 0:
        issues.append(Issue(Severity.ERROR, "frameshift", "payload block is out of frame"))

    return Assembly(design=d, seq=final, protein=protein, segments=segments, issues=issues)


def _translate_orf(seq: str, start: int) -> tuple[str, bool]:
    aa = []
    i = start
    while i + 3 <= len(seq):
        codon = seq[i : i + 3]
        if codon in STOP_CODONS:
            return "".join(aa), True
        aa.append(str(Seq(codon).translate()))
        i += 3
    return "".join(aa), False
