"""Translate the assembled plasmid from the vector's own start codon and check that it is the
protein that was asked for. This is the pass/fail gate; the rest of the package only advises.
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
from .vector import translate_to_stop


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

        # carry the vector's annotations over: those before the replaced stretch stay,
        # those after it move by the change in length, those it overlaps are gone
        delta = len(self.seq) - len(d.vector.seq)
        for f in d.vector.record.features:
            s, e = int(f.location.start), int(f.location.end)
            if f.type == "source" or (s < d.down.end and e > d.up.start):
                continue
            try:
                loc = f.location if e <= d.up.start else f.location + delta
            except (TypeError, ValueError):
                continue
            rec.features.append(SeqFeature(loc, type=f.type, qualifiers=dict(f.qualifiers)))

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
        rec.features.sort(key=lambda f: (int(f.location.start), int(f.location.end)))
        return rec


def assemble(design: InsertDesign) -> Assembly:
    """Splice the insert into the vector and validate the result."""
    d = design
    v = d.vector.seq
    final = v[: d.up.start] + d.insert_dna + v[d.down.end :]

    issues = list(d.issues)
    protein, stop_found = _translate_orf(final, d.orf_start)
    segments = fusion_segments(d, final)

    payload = d.chain.payload
    if payload not in protein:
        issues.append(
            Issue(
                Severity.ERROR,
                "payload_lost",
                "the payload sequence does not appear in the translated fusion "
                "protein: frameshift or premature stop at a junction",
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

    if not stop_found:
        issues.append(
            Issue(Severity.ERROR, "runaway_orf", "no in-frame stop codon found after the insert")
        )

    if "*" in observed:
        issues.append(
            Issue(Severity.ERROR, "internal_stop", "the coding region contains a stop codon")
        )

    # A site definition can name a start codon of its own (e.g. an NdeI ATG), but if the
    # vector's real start codon lies upstream and in frame it still initiates first and
    # the tag it encodes stays on the product.  Translate from the real start and compare.
    real_start = d.vector.expression_start
    declared_atg = final[d.orf_start : d.orf_start + 3] == "ATG"
    if real_start is not None and declared_atg and d.orf_start != real_start:
        real, _ = _translate_orf(final, real_start)
        if real != protein:
            issues.append(
                Issue(
                    Severity.ERROR,
                    "wrong_start",
                    f"the design assumes translation starts at {d.orf_start}, but the "
                    f"vector's own start codon at {real_start} initiates first and gives a "
                    f"different protein (starting {real[:24]}... not {protein[:24]}...)",
                )
            )

    issues.append(
        Issue(
            Severity.INFO,
            "frame",
            f"ORF starts at {d.orf_start}; payload block at {d.coding_start} "
            f"(offset {(d.coding_start - d.orf_start) % 3}, must be 0); "
            f"pads: 5'={len(d.pad_up)} nt, 3'={len(d.pad_down)} nt",
        )
    )
    if (d.coding_start - d.orf_start) % 3 != 0:
        issues.append(Issue(Severity.ERROR, "frameshift", "payload block is out of frame"))

    return Assembly(design=d, seq=final, protein=protein, segments=segments, issues=issues)


def _translate_orf(seq: str, start: int) -> tuple[str, bool]:
    return translate_to_stop(seq, start)
