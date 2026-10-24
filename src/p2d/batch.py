"""Design many proteins into one vector and one pair of sites, and write what was made.

Every protein gets the same vector, sites, host and options, and is designed on its own: a
protein that cannot be cloned (a site it forces, an unusable residue) is reported against its
name and the others carry on.  The outputs are a summary table, the DNA to order, the inserts,
and one annotated GenBank file per construct.
"""

from __future__ import annotations

import csv
import io
import re
from collections.abc import Iterable
from dataclasses import dataclass, field
from pathlib import Path

from Bio import SeqIO

from .assemble import Assembly, assemble
from .enzymes import CutSite, plan_cut
from .frame import CTerm, design_insert
from .host import HostProfile
from .model import AA_ALPHABET, Chain, ConstructPlan, Part
from .optimize import CodonStrategy, Limits
from .order import DEFAULT_CLAMP, Order, make_order
from .report import construct_report
from .vector import Vector

MAX_NAME = 60
SUMMARY_COLUMNS = [
    "name",
    "status",
    "residues",
    "fusion_residues",
    "insert_bp",
    "plasmid_bp",
    "order_bp",
    "order_gc",
    "cai",
    "longest_homopolymer",
    "warnings",
    "issue_codes",
    "error",
    "order_dna",
]


@dataclass
class Protein:
    """One input record: a name, its residues, and why it cannot be designed (if it cannot)."""

    name: str
    sequence: str
    error: str | None = None


@dataclass
class BatchEntry:
    name: str
    protein: str
    assembly: Assembly | None = None
    order: Order | None = None
    error: str | None = None

    @property
    def ok(self) -> bool:
        return (
            self.error is None
            and self.assembly is not None
            and self.assembly.ok
            and self.order is not None
            and self.order.ok
        )

    @property
    def issues(self):
        if self.assembly is None:
            return []
        return [*self.assembly.issues, *(self.order.issues if self.order else [])]

    @property
    def message(self) -> str | None:
        """Why this entry failed: its error, or the first error the checks found."""
        if self.error:
            return self.error
        errors = [i for i in self.issues if i.severity.value == "error"]
        return errors[0].message if errors else None


@dataclass
class BatchResult:
    vector: Vector
    entries: list[BatchEntry] = field(default_factory=list)

    @property
    def ok(self) -> bool:
        return all(e.ok for e in self.entries)

    @property
    def n_ok(self) -> int:
        return sum(e.ok for e in self.entries)


# --------------------------------------------------------------------------- reading


def safe_name(raw: str, fallback: str) -> str:
    """A name usable as a file name: first word of a FASTA header, plain characters only."""
    first = raw.strip().split()[0] if raw.strip() else ""
    name = re.sub(r"[^A-Za-z0-9_.\-]+", "_", first).strip("._-")[:MAX_NAME]
    return name or fallback


def read_proteins(text: str, *, max_residues: int | None = None) -> list[Protein]:
    """Parse FASTA text.  A record that is not a usable protein carries its error."""
    records: list[tuple[str, list[str]]] = []
    for line in text.replace("\r", "").split("\n"):
        if line.startswith(">"):
            records.append((line[1:], []))
        elif line.strip():
            if not records:  # a bare sequence with no header
                records.append(("", []))
            records[-1][1].append(line)
    taken: set[str] = set()
    out: list[Protein] = []
    for i, (header, lines) in enumerate(records, start=1):
        name = safe_name(header, f"protein_{i}")
        base, n = name, 1
        while name.lower() in taken:
            n += 1
            name = f"{base[: MAX_NAME - len(str(n)) - 1]}_{n}"
        taken.add(name.lower())
        seq = re.sub(r"[\s\d]+", "", "".join(lines)).upper().replace("*", "")
        error = None
        if not seq:
            error = "no sequence"
        elif bad := sorted(set(seq) - AA_ALPHABET):
            error = f"unsupported character(s): {', '.join(bad)}"
        elif max_residues and len(seq) > max_residues:
            error = f"{len(seq)} residues; the limit is {max_residues}"
        out.append(Protein(name, seq, error))
    return out


# --------------------------------------------------------------------------- designing


def design_batch(
    proteins: Iterable[Protein],
    vector: Vector,
    up: CutSite,
    down: CutSite,
    host: HostProfile,
    *,
    cterm: CTerm = CTerm.STOP,
    strategy: CodonStrategy = CodonStrategy.WEIGHTED_SAMPLE,
    seed: int | None = None,
    limits: Limits | None = None,
    clamp_5: int = DEFAULT_CLAMP,
    clamp_3: int = DEFAULT_CLAMP,
    n_parts: tuple[Part, ...] = (),
    c_parts: tuple[Part, ...] = (),
    stop_codons: int = 1,
) -> BatchResult:
    """Design each protein; a failure is recorded on its entry and never stops the others."""
    site = plan_cut(vector, up, down)  # a bad pair is the caller's error, not a protein's
    result = BatchResult(vector)
    for p in proteins:
        entry = BatchEntry(p.name, p.sequence, error=p.error)
        if p.error is None:
            try:
                chain = Chain(
                    p.name, p.sequence, n_terminal=list(n_parts), c_terminal=list(c_parts)
                )
                plan = ConstructPlan([chain], host=host.slug)
                entry.assembly = assemble(
                    design_insert(
                        plan,
                        vector,
                        site,
                        host,
                        strategy=strategy,
                        seed=seed,
                        cterm=cterm,
                        limits=limits,
                        stop_codons=stop_codons,
                    )
                )
                entry.order = make_order(entry.assembly, clamp_5, clamp_3)
            except ValueError as exc:
                entry.error = str(exc)
        result.entries.append(entry)
    return result


# --------------------------------------------------------------------------- writing


def _row(e: BatchEntry) -> dict:
    row = dict.fromkeys(SUMMARY_COLUMNS, "")
    row.update(name=e.name, status="ok" if e.ok else "failed", residues=len(e.protein))
    row["error"] = e.message or ""
    if e.assembly is not None:
        a, d = e.assembly, e.assembly.design
        sizes = construct_report(a)["sizes"]
        row.update(
            fusion_residues=len(a.protein),
            insert_bp=len(a.insert_dna),
            plasmid_bp=sizes["plasmid_bp"],
            cai=f"{d.metrics['cai']:.3f}" if d.metrics.get("cai") is not None else "",
            warnings=sum(i.severity.value == "warning" for i in e.issues),
            issue_codes=";".join(sorted({i.code for i in e.issues})),
        )
    if e.order is not None:
        row.update(
            order_bp=e.order.length,
            order_gc=f"{e.order.gc:.3f}",
            longest_homopolymer=e.order.longest_homopolymer,
            order_dna=e.order.dna,
        )
    return row


def summary_rows(result: BatchResult) -> list[dict]:
    return [_row(e) for e in result.entries]


def summary_csv(result: BatchResult) -> str:
    buf = io.StringIO()
    writer = csv.DictWriter(buf, fieldnames=SUMMARY_COLUMNS, lineterminator="\n")
    writer.writeheader()
    writer.writerows(summary_rows(result))
    return buf.getvalue()


def _fasta(result: BatchResult, suffix: str, get) -> str:
    return "".join(
        f">{e.name}_{suffix}\n{get(e)}\n" for e in result.entries if e.ok and get(e) is not None
    )


def orders_fasta(result: BatchResult) -> str:
    return _fasta(result, "order", lambda e: e.order.dna if e.order else None)


def inserts_fasta(result: BatchResult) -> str:
    return _fasta(result, "insert", lambda e: e.assembly.insert_dna if e.assembly else None)


def genbank_text(entry: BatchEntry) -> str:
    buf = io.StringIO()
    record = entry.assembly.to_record()
    record.id = record.name = entry.name[:16]
    SeqIO.write(record, buf, "genbank")
    return buf.getvalue()


def files(result: BatchResult) -> dict[str, str]:
    """Everything a batch makes, as ``{relative path: text}``."""
    out = {
        "summary.csv": summary_csv(result),
        "orders.fasta": orders_fasta(result),
        "inserts.fasta": inserts_fasta(result),
    }
    for e in result.entries:
        if e.ok:
            out[f"constructs/{e.name}.gb"] = genbank_text(e)
    failed = [f"{e.name}\t{e.message}" for e in result.entries if not e.ok]
    if failed:
        out["failed.txt"] = "\n".join(failed) + "\n"
    return out


def write_batch(result: BatchResult, directory: str | Path) -> list[Path]:
    """Write :func:`files` under ``directory`` (created if needed); returns the paths."""
    base = Path(directory)
    written = []
    for rel, text in files(result).items():
        path = base / rel
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(text)
        written.append(path)
    return written
