"""Numbers a bench scientist wants about a finished design: sizes, a digest to check it by,
what was added for the reading frame, and what the fusion protein is.

Lengths are top-strand lengths, the way a gel and Biopython's ``catalyse`` count them: a
cut falls ``fst5`` bases after the start of its recognition site (EcoRI cuts G^AATTC, so 1).
"""

from __future__ import annotations

import re
from dataclasses import dataclass

from Bio.SeqUtils import gc_fraction
from Bio.SeqUtils.ProtParam import ProteinAnalysis

from .assemble import Assembly
from .vector import Vector, _enzyme

# features whose loss from the plasmid should be pointed out to whoever picked the sites
_ESSENTIAL_TYPES = {"rep_origin", "CDS", "promoter", "terminator"}
_ESSENTIAL_LABEL = re.compile(
    r"\b(ori|origin|resistance|kan|amp|cam|tet|bla|lac[iz]|rop|promoter)", re.I
)


@dataclass(frozen=True)
class Fragment:
    start: int  # top-strand position of the cut that opens it, 0-based
    length: int
    role: str  # "backbone (kept)", "replaced piece", "insert", "fragment"


def digest(seq: str, circular: bool, enzyme_names: list[str], kept_start: int | None = None):
    """Complete digest of ``seq`` by ``enzyme_names``.  ``kept_start`` marks the fragment that
    opens at that cut position as the backbone and the rest as replaced."""
    vec = Vector.from_sequence(seq, name="digest", circular=circular)
    cuts = set()
    for name in enzyme_names:
        offset = _enzyme(name).fst5
        cuts.update(h.start + offset for h in vec.find_all(name))
    cuts = sorted(cuts)
    if not cuts:
        return [Fragment(0, len(seq), "fragment")]
    if not circular:
        edges = [0, *cuts, len(seq)]
        return [Fragment(a, b - a, "fragment") for a, b in zip(edges, edges[1:], strict=False)]
    out = []
    for i, cut in enumerate(cuts):
        nxt = cuts[(i + 1) % len(cuts)]
        length = (nxt - cut) % len(seq) or len(seq)
        role = "fragment"
        if kept_start is not None:
            role = "backbone (kept)" if cut == kept_start else "replaced piece"
        out.append(Fragment(cut, length, role))
    return out


def replaced_features(vector: Vector, lo: int, hi: int) -> list[dict]:
    """Annotated features inside the replaced stretch, so a careless choice is visible."""
    out = []
    for f in vector.record.features:
        s, e = int(f.location.start), int(f.location.end)
        if f.type == "source" or e <= lo or s >= hi:
            continue
        label = (f.qualifiers.get("label") or f.qualifiers.get("gene") or [f.type])[0]
        out.append(
            {
                "label": str(label)[:40],
                "type": f.type,
                "start": s,
                "end": e,
                "fully_removed": lo <= s and e <= hi,
                "essential": f.type in _ESSENTIAL_TYPES
                or bool(_ESSENTIAL_LABEL.search(str(label))),
            }
        )
    return out


def construct_report(a: Assembly) -> dict:
    d = a.design
    v = d.vector
    up_name, down_name = d.site.upstream_enzyme, d.site.downstream_enzyme
    cut_up = d.up.start + _enzyme(up_name).fst5
    cut_down = d.down.start + _enzyme(down_name).fst5
    up_len, down_len = len(d.up), len(d.down)

    replaced = cut_down - cut_up
    backbone = len(v.seq) - replaced
    insert_bp = len(d.insert_dna)
    cut_insert = (insert_bp - down_len + _enzyme(down_name).fst5) - _enzyme(up_name).fst5

    stop = "TAATGA" if d.stop_in_insert else ""
    pad_up, pad_down = d.pad_up, d.pad_down
    parts = [
        {"name": f"5' site ({up_name})", "dna": d.insert_dna[:up_len]},
        {"name": "frame pad", "dna": pad_up},
        {"name": "coding region", "bp": len(d.coding)},
        {"name": "stop codons", "dna": stop},
        {"name": "frame pad", "dna": pad_down},
        {"name": f"3' site ({down_name})", "dna": d.insert_dna[insert_bp - down_len :]},
    ]
    parts = [p for p in parts if p.get("dna") or p.get("bp")]

    names = sorted({up_name, down_name})
    original = digest(v.seq, v.circular, names, kept_start=cut_down)
    final = digest(a.seq, v.circular, names)
    final_cut_up = d.up.start + _enzyme(up_name).fst5
    final = [
        Fragment(
            f.start,
            f.length,
            "insert" if f.start == final_cut_up else "backbone (kept)",
        )
        for f in final
    ]

    protein = ProteinAnalysis(a.protein)
    return {
        "sizes": {
            "vector_bp": len(v.seq),
            "replaced_bp": replaced,
            "backbone_bp": backbone,
            "insert_bp": insert_bp,
            "insert_after_digest_bp": cut_insert,
            "plasmid_bp": len(a.seq),
            "coding_bp": len(d.coding),
            "net_change_bp": len(a.seq) - len(v.seq),
        },
        "insert_parts": parts,
        "frame": {
            "pad_5prime": pad_up,
            "pad_3prime": pad_down,
            "start_codon": d.orf_start,
        },
        "digest_check": {
            "enzymes": names,
            "vector": [f.__dict__ for f in sorted(original, key=lambda f: -f.length)],
            "final": [f.__dict__ for f in sorted(final, key=lambda f: -f.length)],
        },
        "protein": {
            "residues": len(a.protein),
            "payload_residues": len(d.chain.insert_aa),
            "mass_kda": round(protein.molecular_weight() / 1000, 2),
            "pi": round(protein.isoelectric_point(), 2),
        },
        "gc": {
            "plasmid": round(gc_fraction(a.seq), 4),
            "coding": round(gc_fraction(d.coding), 4),
        },
        "replaced_features": replaced_features(v, d.up.start, d.down.end),
        "ends": {
            "five_prime": f"{_enzyme(up_name).ovhgseq} ({up_name})",
            "three_prime": f"{_enzyme(down_name).ovhgseq} ({down_name})",
        },
    }
