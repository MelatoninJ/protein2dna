"""JSON-friendly wrapper around the p2d library, used by the local web UI.

No logic lives in the UI: this module only validates input, calls the library,
and flattens the result into plain dicts.  It has no web-framework dependency so
the same functions can sit behind any server.
"""

from __future__ import annotations

import io
import re

from Bio import SeqIO

import p2d
from p2d import CodonStrategy, ConstructPlan, available_hosts, design_insert, load_host
from p2d.assemble import assemble
from p2d.enzymes import candidate_enzymes, expression_window, plan_site
from p2d.optimize import OptimizationRequest, optimize_coding
from p2d.vectors import available_bundled, load_bundled

MAX_RESIDUES = 5000
RESIDUES = set("ACDEFGHIKLMNPQRSTVWY")
STRATEGIES = {
    "weighted_sample": CodonStrategy.WEIGHTED_SAMPLE,
    "max_cai": CodonStrategy.MAX_CAI,
}


class InputError(ValueError):
    """The request is unusable; the message is safe to show to the user."""


def parse_protein(text: str) -> str:
    """Accept a bare sequence or FASTA text; return clean one-letter residues.

    Only the first FASTA record is used.  Digits and whitespace (GenPept-style
    numbering) are ignored; a trailing ``*`` is dropped like the library does.
    """
    lines = text.replace("\r", "").split("\n")
    seq_lines: list[str] = []
    seen_header = False
    for line in lines:
        if line.startswith(">"):
            if seen_header or seq_lines:
                break
            seen_header = True
            continue
        seq_lines.append(line)
    seq = re.sub(r"[\s\d]+", "", "".join(seq_lines)).upper().replace("*", "")
    if not seq:
        raise InputError("Paste a protein sequence to design.")
    bad = sorted(set(seq) - RESIDUES)
    if bad:
        raise InputError(
            f"Unsupported character(s): {', '.join(bad)}. "
            "Use the 20 standard one-letter amino acid codes."
        )
    if len(seq) > MAX_RESIDUES:
        raise InputError(f"Sequence is {len(seq)} residues; the limit is {MAX_RESIDUES}.")
    return seq


def meta() -> dict:
    vectors = []
    for name in available_bundled():
        v = load_bundled(name)
        vectors.append(
            {
                "name": name,
                "length": len(v),
                "picker": v.expression_start is not None,
                "sites": [{"name": k, "note": s.note} for k, s in v.sites.items()],
            }
        )
    hosts = []
    for slug in available_hosts():
        h = load_host(slug)
        hosts.append({"slug": slug, "name": h.name, "verified": h.codon_table_verified})
    return {
        "version": p2d.__version__,
        "vectors": vectors,
        "hosts": hosts,
        "strategies": [
            {
                "value": "weighted_sample",
                "label": "Balanced",
                "note": (
                    "Samples common codons, avoiding the tRNA depletion "
                    "of always taking the top one."
                ),
            },
            {
                "value": "max_cai",
                "label": "Maximum CAI",
                "note": "Always the most-used codon that keeps every site out.",
            },
        ],
        "max_residues": MAX_RESIDUES,
    }


def _protein_conflict(protein: str, site: str, host_slug: str) -> str | None:
    """Why ``site`` cannot be kept out of ``protein``'s coding DNA, or None if it can.

    A site is only truly unavoidable when every synonymous encoding creates it
    (a run of single-codon residues); this asks the optimiser, which knows.
    """
    result = optimize_coding(
        OptimizationRequest(
            aa=protein,
            host=load_host(host_slug),
            prefix="",
            suffix="",
            forbidden=[site],
            strategy=CodonStrategy.MAX_CAI,
        )
    )
    if result.ok:
        return None
    m = re.search(r"residue (\d+)", next((i.message for i in result.issues), ""))
    where = f" near residue {m.group(1)}" if m else ""
    return f"Your protein forces the {site} site{where}; no synonymous codons avoid it."


def _try_protein(text) -> str | None:
    """The cleaned protein, or None while the user is still typing something unusable."""
    try:
        return parse_protein(str(text or ""))
    except InputError:
        return None


def enzymes(request: dict) -> dict:
    """Single-cutter enzymes in the vector's expression ORF, and which can be used.

    Every enzyme is judged as a 5' choice.  When ``upstream`` is given each one is
    also judged as the 3' partner of it, by running the real design.
    """
    vector_name = str(request.get("vector") or "pTEST1")
    if vector_name not in available_bundled():
        raise InputError(f"Unknown vector {vector_name!r}.")
    vector = load_bundled(vector_name)
    host_slug = str(request.get("host") or "ecoli_bl21")
    if host_slug not in available_hosts():
        raise InputError(f"Unknown host {host_slug!r}.")
    try:
        window = expression_window(vector)
        candidates = candidate_enzymes(vector)
    except ValueError as exc:
        raise InputError(str(exc)) from exc

    protein = _try_protein(request.get("protein"))
    upstream = str(request.get("upstream") or "")
    host = load_host(host_slug)

    def partners(name: str) -> list:
        out = []
        for other in candidates:
            try:
                plan_site(vector, name, other.name)
            except ValueError:
                continue
            out.append(other)
        return out

    chosen = next((c for c in candidates if upstream in (c.name, *c.aliases)), None)
    rows = []
    for c in candidates:
        conflict = _protein_conflict(protein, c.site, host_slug) if protein else None
        if conflict:
            up = {"ok": False, "reason": conflict}
        elif not partners(c.name):
            up = {
                "ok": False,
                "reason": "No enzyme downstream of it can pair with it: the next sites "
                "touch it or leave the same ends.",
            }
        else:
            up = {"ok": True, "reason": None}
        site_in_start = c.start < window.start + 3
        up["starts_at"] = "site" if site_in_start else "vector"

        down = None
        if upstream and chosen is not None and c.start <= chosen.start:
            down = {"ok": False, "reason": None, "hidden": True}  # wrong side: not a candidate
        elif upstream:
            try:
                site = plan_site(vector, upstream, c.name)
            except ValueError as exc:
                down = {"ok": False, "reason": str(exc)}
            else:
                down = {"ok": True, "reason": None}
                if protein:
                    try:
                        a = assemble(
                            design_insert(
                                ConstructPlan.single_chain(protein, host=host_slug),
                                vector,
                                site,
                                host,
                                strategy=CodonStrategy.MAX_CAI,
                            )
                        )
                    except ValueError as exc:
                        down = {"ok": False, "reason": str(exc)}
                    else:
                        if not a.ok:
                            down = {"ok": False, "reason": a.errors[0].message}
        rows.append(
            {
                "name": c.name,
                "aliases": list(c.aliases),
                "site": c.site,
                "start": c.start,
                "end": c.end,
                "overhang_kind": c.overhang_kind,
                "overhang": c.overhang,
                "as_upstream": up,
                "as_downstream": down,
            }
        )
    return {
        "vector": vector.name,
        "window": {"start": window.start, "stop_end": window.stop_end},
        "protein_length": len(protein) if protein else None,
        "enzymes": rows,
    }


def design(request: dict) -> dict:
    """Run one design and return everything the UI renders."""
    protein = parse_protein(str(request.get("protein", "")))

    vector_name = str(request.get("vector") or "pTEST1")
    if vector_name not in available_bundled():
        raise InputError(f"Unknown vector {vector_name!r}.")
    vector = load_bundled(vector_name)

    upstream = str(request.get("upstream") or "")
    downstream = str(request.get("downstream") or "")
    if upstream or downstream:
        if not (upstream and downstream):
            raise InputError("Choose both a 5' and a 3' enzyme.")
        try:
            site = plan_site(vector, upstream, downstream)
        except ValueError as exc:
            raise InputError(str(exc)) from exc
        site_name = site.name
    else:
        site_name = str(request.get("site") or "")
        if not site_name:
            raise InputError(f"Choose a 5' and a 3' enzyme for {vector_name}.")
        if site_name not in vector.sites:
            raise InputError(f"{vector_name} has no insertion site {site_name!r}.")
        site = vector.sites[site_name]

    host_slug = str(request.get("host") or "ecoli_bl21")
    if host_slug not in available_hosts():
        raise InputError(f"Unknown host {host_slug!r}.")

    strategy_key = str(request.get("strategy") or "weighted_sample")
    if strategy_key not in STRATEGIES:
        raise InputError(f"Unknown strategy {strategy_key!r}.")

    seed = request.get("seed")
    if seed in (None, ""):
        seed = None
    else:
        try:
            seed = int(seed)
        except (TypeError, ValueError):
            raise InputError("Seed must be a whole number.") from None
        if not 0 <= seed < 2**32:
            raise InputError("Seed must be between 0 and 4294967295.")

    host = load_host(host_slug)
    plan = ConstructPlan.single_chain(protein, host=host_slug)
    try:
        a = assemble(
            design_insert(
                plan,
                vector,
                site,
                host,
                strategy=STRATEGIES[strategy_key],
                seed=seed,
            )
        )
    except ValueError as exc:
        # frame-engine refusals (e.g. no pad that avoids a stop) are user-facing
        raise InputError(str(exc)) from exc

    d = a.design
    coding_in_insert = d.coding_start - d.up.start
    segments, pos = [], 1
    for seg in a.segments:
        segments.append(
            {
                "label": seg.label,
                "aa": seg.aa,
                "kind": seg.kind.value,
                "source": seg.source.value,
                "start": pos,
                "end": pos + len(seg) - 1,
            }
        )
        pos += len(seg)

    gb = io.StringIO()
    SeqIO.write(a.to_record(), gb, "genbank")
    name = f"{vector.name}_{d.chain.name}"
    return {
        "ok": a.ok,
        "vector": vector.name,
        "site": site_name,
        "host": host_slug,
        "payload_length": len(protein),
        "protein": a.protein,
        "segments": segments,
        "issues": [
            {"severity": i.severity.value, "code": i.code, "message": i.message}
            for i in sorted(
                a.issues, key=lambda i: ("error", "warning", "info").index(i.severity.value)
            )
        ],
        "insert_dna": a.insert_dna,
        "coding_span": [coding_in_insert, coding_in_insert + len(d.coding)],
        "insert_bp": len(a.insert_dna),
        "plasmid_bp": len(a.seq),
        "metrics": {
            "cai": d.metrics.get("cai"),
            "gc": d.metrics.get("gc"),
            "n_rare": d.metrics.get("n_rare"),
            "strategy": d.metrics.get("strategy"),
        },
        "genbank": gb.getvalue(),
        "fasta": f">{name}_insert\n{a.insert_dna}\n",
        "filename": name,
    }
