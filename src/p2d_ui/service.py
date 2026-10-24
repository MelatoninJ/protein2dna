"""JSON-friendly wrapper around the library for the web UI.

Validates input, calls the library and flattens the results to dicts, with no logic of its own
and no web framework. A loaded plasmid is kept in memory under a random id so the browser sends
it once; nothing is written to disk.
"""

from __future__ import annotations

import base64
import binascii
import io
import re
import secrets
import threading
import zipfile
from collections import OrderedDict
from dataclasses import dataclass

from Bio import SeqIO

from p2d import CodonStrategy, ConstructPlan, CTerm, available_hosts, design_insert, load_host
from p2d import batch as batching
from p2d.analysis import VectorAnalysis, analyze_vector, feature_label
from p2d.assemble import assemble
from p2d.enzymes import CutSite, find_cut_site, list_cut_sites, plan_cut, upstream_problem
from p2d.frame import find_cterm_tag
from p2d.model import LINKERS, PROTEASE_SITES, TAGS, Chain, resolve_parts
from p2d.optimize import OptimizationRequest, optimize_coding
from p2d.order import DEFAULT_CLAMP, MAX_CLAMP, make_order
from p2d.recommend import group_by_outcome, rank_pairs
from p2d.report import construct_report
from p2d.vector import Vector
from p2d.vectorio import (
    VectorFileError,
    list_library,
    read_library_vector,
    read_vector_bytes,
    read_vector_text,
)
from p2d.vectors import available_bundled, load_bundled

MAX_RESIDUES = 5000
MAX_STORED = 16
MAX_FEATURES = 150
RESIDUES = set("ACDEFGHIKLMNPQRSTVWY")
STRATEGIES = {
    "weighted_sample": CodonStrategy.WEIGHTED_SAMPLE,
    "max_cai": CodonStrategy.MAX_CAI,
}
CTERMS = {"stop": CTerm.STOP, "tag": CTerm.TAG, "vector_frame": CTerm.VECTOR_FRAME}
MAX_GROUPS = 12  # distinct results returned by pairs()
MAX_PAIRS_PER_GROUP = 12
WINDOW_PAST_ORF = 150  # bp shown beyond the end of the expression ORF by default
HOLDS_START_SLACK = 10  # bp shown before the start codon (NcoI-style sites begin there)


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


def _try_protein(text) -> str | None:
    """The cleaned protein, or None while the user is still typing something unusable."""
    try:
        return parse_protein(str(text or ""))
    except InputError:
        return None


@dataclass
class _Entry:
    analysis: VectorAnalysis
    source: str
    sites: list[CutSite]


_STORE: OrderedDict[str, _Entry] = OrderedDict()
_LOCK = threading.Lock()


def _remember(analysis: VectorAnalysis, source: str) -> str:
    entry = _Entry(analysis, source, list_cut_sites(analysis.vector))
    vector_id = secrets.token_urlsafe(9)
    with _LOCK:
        _STORE[vector_id] = entry
        while len(_STORE) > MAX_STORED:
            _STORE.popitem(last=False)
    return vector_id


def _entry(vector_id) -> _Entry:
    with _LOCK:
        entry = _STORE.get(str(vector_id))
        if entry:
            _STORE.move_to_end(str(vector_id))
    if entry is None:
        raise InputError("That plasmid is no longer loaded. Load it again.")
    return entry


def _with_start(entry: _Entry, start) -> Vector:
    """A copy of the oriented vector with the confirmed start codon set."""
    if start in (None, ""):
        raise InputError("Choose the start codon first.")
    try:
        start = int(start)
    except (TypeError, ValueError):
        raise InputError("The start codon position must be a whole number.") from None
    base = entry.analysis.vector
    if not 0 <= start <= len(base.seq) - 3 or base.seq[start : start + 3] != "ATG":
        raise InputError(f"Position {start + 1} is not an ATG.")
    v = Vector(record=base.record, name=base.name, circular=base.circular)
    v.expression_start = start
    return v


MAX_PARTS = 4
MAX_BATCH = 100  # proteins per batch request
DEFAULT_HOST = "ecoli_bl21"


def _part_keys(request: dict, key: str) -> list[str]:
    raw = request.get(key) or []
    if (
        not isinstance(raw, list)
        or len(raw) > MAX_PARTS
        or not all(isinstance(k, str) for k in raw)
    ):
        raise InputError(f"{key} must be a short list of tag or linker names.")
    try:
        resolve_parts(raw)
    except ValueError as exc:
        raise InputError(str(exc).capitalize() + ".") from exc
    return raw


def _plan(protein: str, host_slug: str, request: dict) -> ConstructPlan:
    """The protein with any tags and linkers asked for, N-terminal ones first."""
    chain = Chain(
        "payload",
        protein,
        n_terminal=resolve_parts(_part_keys(request, "n_parts")),
        c_terminal=resolve_parts(_part_keys(request, "c_parts")),
    )
    return ConstructPlan([chain], host=host_slug)


LIMIT_FIELDS = {  # request key -> (Limits fields set, multiplier from the unit the user types)
    "gc_min": (("gc_min",), 0.01),
    "gc_max": (("gc_max",), 0.01),
    "window_min": (("window_min",), 0.01),
    "window_max": (("window_max",), 0.01),
    "window_size": (("window_size",), 1),
    "repeat_max": (("direct_max",), 1),
    "inverted_max": (("inverted_max",), 1),
}


def _limits(request: dict, host):
    """The host's GC and repeat limits with any the user changed (percent, bases)."""
    given = request.get("limits") or {}
    if not isinstance(given, dict):
        raise InputError("Limits must be a set of named numbers.")
    override = {}
    for key, raw in given.items():
        if key not in LIMIT_FIELDS:
            raise InputError(f"Unknown limit {key!r}.")
        if raw in (None, ""):
            continue
        if isinstance(raw, bool) or not isinstance(raw, int | float):
            raise InputError(f"{key} must be a number.")
        fields, unit = LIMIT_FIELDS[key]
        for field in fields:
            override[field] = int(raw) if unit == 1 else round(raw * unit, 4)
    try:
        return host.limits(**override)
    except ValueError as exc:
        raise InputError(str(exc)[:1].upper() + str(exc)[1:] + ".") from exc


def _stop_codons(request: dict) -> int:
    raw = request.get("stop_codons")
    if raw in (None, ""):
        return 1
    if isinstance(raw, bool) or str(raw) not in ("1", "2"):
        raise InputError("Stop codons must be 1 or 2.")
    return int(raw)


def _catalog() -> dict:
    def rows(table):
        return [{"key": k, "name": p.name, "aa": p.aa} for k, p in table.items()]

    return {"tags": rows(TAGS), "linkers": rows(LINKERS), "cleavage": rows(PROTEASE_SITES)}


def meta() -> dict:
    hosts = []
    for slug in sorted(available_hosts(), key=lambda s: (load_host(s).status != "released", s)):
        h = load_host(slug)
        lim = h.limits()
        hosts.append(
            {
                "slug": slug,
                "name": h.name,
                "verified": h.codon_table_verified,
                "status": h.status,
                "note": str(h.raw.get("codon_usage", {}).get("note", "")),
                "limits": {
                    "gc_min": round(lim.gc_min * 100),
                    "gc_max": round(lim.gc_max * 100),
                    "window_size": lim.window_size,
                    "window_min": round(lim.window_min * 100),
                    "window_max": round(lim.window_max * 100),
                    "repeat_max": lim.direct_max,
                },
            }
        )
    return {
        "bundled": [
            {"id": n, "name": n, "length": len(load_bundled(n))} for n in available_bundled()
        ],
        "library": list_library(),
        "hosts": hosts,
        "default_host": DEFAULT_HOST,
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
        "parts": _catalog(),
    }


def load_vector(request: dict) -> dict:
    """Load a plasmid from the bundled set, the user's folder, an upload or pasted text."""
    try:
        if request.get("bundled"):
            name = str(request["bundled"])
            if name not in available_bundled():
                raise InputError(f"Unknown bundled vector {name!r}.")
            vector, source = load_bundled(name), "bundled"
        elif request.get("library"):
            vector, source = read_library_vector(str(request["library"])), "your vector folder"
        elif request.get("data_b64"):
            try:
                data = base64.b64decode(str(request["data_b64"]), validate=True)
            except (binascii.Error, ValueError):
                raise InputError("The uploaded file could not be decoded.") from None
            vector = read_vector_bytes(data, str(request.get("filename") or "upload"))
            source = "uploaded file"
        elif request.get("text"):
            vector = read_vector_text(
                str(request["text"]), str(request.get("filename") or "pasted")
            )
            source = "pasted text"
        else:
            raise InputError("Choose a plasmid, upload a file, or paste a sequence.")
    except VectorFileError as exc:
        raise InputError(str(exc)) from exc

    analysis = analyze_vector(vector)
    v = analysis.vector
    features = []
    for f in sorted(v.record.features, key=lambda f: int(f.location.start)):
        if f.type == "source":
            continue
        features.append(
            {
                "label": feature_label(f),
                "type": f.type,
                "start": int(f.location.start),
                "end": int(f.location.end),
                "strand": f.location.strand or 0,
            }
        )
    sites = list_cut_sites(v)
    return {
        "id": _remember(analysis, source),
        "name": v.name,
        "source": source,
        "length": len(v),
        "circular": v.circular,
        "flipped": analysis.flipped,
        "orientation_notes": analysis.orientation_notes,
        "notes": list(v.notes),
        "starts": [
            {
                "position": s.position,
                "confidence": s.confidence,
                "evidence": list(s.evidence),
                "n_terminal": s.n_terminal,
                "orf_end": s.orf_end,
            }
            for s in analysis.starts
        ],
        "mcs": analysis.region and {"start": analysis.region.start, "end": analysis.region.end},
        "features": features[:MAX_FEATURES],
        "warnings": analysis.warnings,
        "site_count": len(sites),
        "single_cutters": sum(1 for s in sites if s.count == 1),
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


def _site_row(c: CutSite) -> dict:
    return {
        "id": c.id,
        "enzyme": c.enzyme,
        "aliases": list(c.aliases),
        "site": c.site,
        "start": c.start,
        "end": c.end,
        "number": c.number,
        "count": c.count,
        "overhang_kind": c.overhang_kind,
        "overhang": c.overhang,
    }


def _window(entry: _Entry, v: Vector | None, window) -> tuple[int, int]:
    """The stretch of plasmid to consider: the user's window, else the expression region."""
    base = entry.analysis.vector
    if window:
        try:
            lo, hi = int(window[0]), int(window[1])
        except (TypeError, ValueError, IndexError):
            raise InputError("The window must be two whole numbers.") from None
        lo, hi = max(0, lo), min(len(base.seq), hi)
        if lo >= hi:
            raise InputError("The window must end after it starts.")
        return lo, hi
    if v is None:
        return 0, len(base.seq)
    orf_end = next(
        (s.orf_end for s in entry.analysis.starts if s.position == v.expression_start), None
    )
    # a site that holds the start codon begins a few bases before it
    lo = max(0, v.expression_start - HOLDS_START_SLACK)
    return lo, min(len(base.seq), (orf_end or lo + 600) + WINDOW_PAST_ORF)


def sites(request: dict) -> dict:
    """Sites in a window, judged as 5' sites and, once one is chosen, as its 3' partner."""
    entry = _entry(request.get("vector_id"))
    base = entry.analysis.vector
    host_slug = str(request.get("host") or DEFAULT_HOST)
    if host_slug not in available_hosts():
        raise InputError(f"Unknown host {host_slug!r}.")
    cterm = CTERMS.get(str(request.get("cterm") or "stop"))
    if cterm is None:
        raise InputError("Unknown C-terminus choice.")

    start = request.get("start")
    v = _with_start(entry, start) if start not in (None, "") else None
    lo, hi = _window(entry, v, request.get("window"))

    shown = [c for c in entry.sites if lo <= c.start and c.end <= hi]
    if not request.get("include_multi", True):
        shown = [c for c in shown if c.count == 1]

    protein = _try_protein(request.get("protein"))
    conflicts: dict[str, str | None] = {}

    def conflict(c: CutSite) -> str | None:
        if protein and c.site not in conflicts:
            conflicts[c.site] = _protein_conflict(protein, c.site, host_slug)
        return conflicts.get(c.site)

    up_id = str(request.get("upstream") or "")
    chosen_up = None
    if up_id and v is not None:
        try:
            chosen_up = find_cut_site(base, up_id)
        except ValueError as exc:
            raise InputError(str(exc)) from exc

    host = load_host(host_slug)
    rows = []
    for c in shown:
        if v is None:
            up = {"ok": False, "reason": "Confirm the start codon first."}
        elif reason := (upstream_problem(v, c) or conflict(c)):
            up = {"ok": False, "reason": reason}
        elif not any(o.start > c.end for o in shown):
            up = {"ok": False, "reason": "No other site lies downstream of it in this window."}
        else:
            up = {"ok": True, "reason": None}

        down = None
        if chosen_up is not None:
            if c.start <= chosen_up.start:
                down = {"ok": False, "reason": None, "hidden": True}
            else:
                down = _downstream_verdict(
                    v, chosen_up, c, protein, host, host_slug, cterm, conflict, request
                )
        rows.append({**_site_row(c), "as_upstream": up, "as_downstream": down})

    return {
        "vector_id": request.get("vector_id"),
        "start": v.expression_start if v else None,
        "window": {"start": lo, "end": hi},
        "sites": rows,
        "protein_length": len(protein) if protein else None,
    }


def _downstream_verdict(v, up, c, protein, host, host_slug, cterm, conflict, request) -> dict:
    try:
        site = plan_cut(v, up, c)
    except ValueError as exc:
        return {"ok": False, "reason": str(exc), "has_tag": False}
    tag = find_cterm_tag(v, c.end)
    verdict = {
        "ok": True,
        "reason": None,
        "has_tag": tag is not None,
        "tag_name": tag.name if tag else None,
    }
    if reason := conflict(c):
        return {**verdict, "ok": False, "reason": reason}
    if protein:
        try:
            a = assemble(
                design_insert(
                    _plan(protein, host_slug, request),
                    v,
                    site,
                    host,
                    strategy=CodonStrategy.MAX_CAI,
                    cterm=cterm,
                    limits=_limits(request, host),
                    stop_codons=_stop_codons(request),
                )
            )
        except ValueError as exc:
            return {**verdict, "ok": False, "reason": str(exc)}
        if not a.ok:
            return {**verdict, "ok": False, "reason": a.errors[0].message}
    return verdict


def pairs(request: dict) -> dict:
    """Every workable pair of sites for the protein, grouped by the protein each one gives."""
    entry = _entry(request.get("vector_id"))
    v = _with_start(entry, request.get("start"))
    protein = parse_protein(str(request.get("protein", "")))
    host_slug = str(request.get("host") or DEFAULT_HOST)
    if host_slug not in available_hosts():
        raise InputError(f"Unknown host {host_slug!r}.")
    cterm = CTERMS.get(str(request.get("cterm") or "stop"))
    if cterm is None:
        raise InputError("Unknown C-terminus choice.")
    keep = {"any": None, "keep": True, "remove": False}.get(str(request.get("n_tag") or "any"))
    if str(request.get("n_tag") or "any") not in ("any", "keep", "remove"):
        raise InputError("N-terminal tag must be any, keep or remove.")
    lo, hi = _window(entry, v, request.get("window"))

    try:
        options, tried = rank_pairs(
            v,
            protein,
            load_host(host_slug),
            window=(lo, hi),
            cterm=cterm,
            n_tag=keep,
            limits=_limits(request, load_host(host_slug)),
            stop_codons=_stop_codons(request),
            n_parts=tuple(resolve_parts(_part_keys(request, "n_parts"))),
            c_parts=tuple(resolve_parts(_part_keys(request, "c_parts"))),
        )
    except ValueError as exc:
        raise InputError(str(exc)) from exc
    groups = group_by_outcome(options)

    def row(o):
        return {
            "up": o.up.id,
            "down": o.down.id,
            "up_label": o.up.label,
            "down_label": o.down.label,
            "single_cutters": o.single_cutters,
            "replaced_bp": o.replaced_bp,
            "plasmid_bp": o.plasmid_bp,
        }

    return {
        "window": {"start": lo, "end": hi},
        "considered": tried,
        "workable": len(options),
        "total_groups": len(groups),
        "groups": [
            {
                "before": g.before,
                "after": g.after,
                "n_tag": g.n_tag,
                "altered": g.altered,
                "pad_bases": g.pad_bases,
                "n_pairs": len(g.pairs),
                "pairs": [row(o) for o in g.pairs[:MAX_PAIRS_PER_GROUP]],
            }
            for g in groups[:MAX_GROUPS]
        ],
    }


@dataclass
class _Settings:
    vector: Vector
    host_slug: str
    host: object
    strategy_key: str
    cterm_key: str
    seed: int | None
    clamps: list[int]
    up: CutSite
    down: CutSite
    site: object
    limits: object
    stop_codons: int


def _settings(request: dict) -> _Settings:
    """Everything a design needs besides the protein, checked the same way for one or many."""
    entry = _entry(request.get("vector_id"))
    v = _with_start(entry, request.get("start"))

    host_slug = str(request.get("host") or DEFAULT_HOST)
    if host_slug not in available_hosts():
        raise InputError(f"Unknown host {host_slug!r}.")
    strategy_key = str(request.get("strategy") or "weighted_sample")
    if strategy_key not in STRATEGIES:
        raise InputError(f"Unknown strategy {strategy_key!r}.")
    cterm_key = str(request.get("cterm") or "stop")
    if cterm_key not in CTERMS:
        raise InputError(f"Unknown C-terminus choice {cterm_key!r}.")

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

    clamps = []
    for key in ("clamp_5", "clamp_3"):
        raw = request.get(key)
        if isinstance(raw, bool) or (isinstance(raw, float) and not raw.is_integer()):
            raise InputError("Extra bases must be whole numbers.")
        try:
            n = DEFAULT_CLAMP if raw in (None, "") else int(raw)
        except (TypeError, ValueError):
            raise InputError("Extra bases must be whole numbers.") from None
        if not 0 <= n <= MAX_CLAMP:
            raise InputError(f"Extra bases must be between 0 and {MAX_CLAMP}.")
        clamps.append(n)

    if not (request.get("upstream") and request.get("downstream")):
        raise InputError("Choose both a 5' and a 3' site.")
    try:
        up = find_cut_site(v, str(request["upstream"]))
        down = find_cut_site(v, str(request["downstream"]))
        site = plan_cut(v, up, down)
    except ValueError as exc:
        raise InputError(str(exc)) from exc

    host = load_host(host_slug)
    return _Settings(
        v,
        host_slug,
        host,
        strategy_key,
        cterm_key,
        seed,
        clamps,
        up,
        down,
        site,
        _limits(request, host),
        _stop_codons(request),
    )


def design(request: dict) -> dict:
    """Run one design and return everything the UI renders."""
    protein = parse_protein(str(request.get("protein", "")))
    s = _settings(request)
    v, host_slug, host, limits = s.vector, s.host_slug, s.host, s.limits
    strategy_key, cterm_key, seed, clamps = s.strategy_key, s.cterm_key, s.seed, s.clamps
    up, down, site = s.up, s.down, s.site
    plan = _plan(protein, host_slug, request)
    try:
        a = assemble(
            design_insert(
                plan,
                v,
                site,
                host,
                strategy=STRATEGIES[strategy_key],
                seed=seed,
                cterm=CTERMS[cterm_key],
                limits=limits,
                stop_codons=s.stop_codons,
            )
        )
    except ValueError as exc:
        # frame-engine refusals (e.g. no pad that avoids a stop) are user-facing
        raise InputError(str(exc)) from exc

    d = a.design
    order = make_order(a, *clamps)
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
    name = f"{v.name}_{up.enzyme}-{down.enzyme}"
    return {
        "ok": a.ok and order.ok,
        "vector": v.name,
        "site": site.name,
        "upstream": up.label,
        "downstream": down.label,
        "host": host_slug,
        "cterm": cterm_key,
        "payload_length": len(protein),
        "protein": a.protein,
        "segments": segments,
        "issues": [
            {"severity": i.severity.value, "code": i.code, "message": i.message}
            for i in sorted(
                [*a.issues, *order.issues],
                key=lambda i: ("error", "warning", "info").index(i.severity.value),
            )
        ],
        "insert_dna": a.insert_dna,
        "coding_span": [coding_in_insert, coding_in_insert + len(d.coding)],
        "order": {
            "dna": order.dna,
            "clamp_5": order.clamp_5,
            "clamp_3": order.clamp_3,
            "insert_span": [order.insert_start, order.insert_end],
            "coding_span": list(order.coding_span),
            "length": order.length,
            "gc": order.gc,
            "longest_homopolymer": order.longest_homopolymer,
            "ligation": {"ok": order.ligation.ok, "message": order.ligation.message},
        },
        "insert_bp": len(a.insert_dna),
        "plasmid_bp": len(a.seq),
        "metrics": {
            "cai": d.metrics.get("cai"),
            "gc": d.metrics.get("gc"),
            "n_rare": d.metrics.get("n_rare"),
            "gc_window_low": d.metrics.get("gc_window_low"),
            "gc_window_high": d.metrics.get("gc_window_high"),
            "repeat_bases": d.metrics.get("repeat_bases"),
            "refined_swaps": d.metrics.get("refined_swaps"),
            "window_size": limits.window_size,
            "strategy": d.metrics.get("strategy"),
        },
        "report": construct_report(a),
        "genbank": gb.getvalue(),
        "fasta": f">{name}_order\n{order.dna}\n",
        "insert_fasta": f">{name}_insert\n{a.insert_dna}\n",
        "filename": name,
    }


def _zip_b64(result) -> str:
    """The whole batch as a zip (summary, orders, inserts, one GenBank each), base64 encoded."""
    buf = io.BytesIO()
    with zipfile.ZipFile(buf, "w", zipfile.ZIP_DEFLATED) as z:
        for rel, text in batching.files(result).items():
            info = zipfile.ZipInfo(rel, date_time=(2020, 1, 1, 0, 0, 0))  # same bytes every time
            info.compress_type = zipfile.ZIP_DEFLATED
            z.writestr(info, text)
    return base64.b64encode(buf.getvalue()).decode("ascii")


def batch(request: dict) -> dict:
    """Design every protein in a FASTA text with one set of options; failures stay per protein."""
    text = request.get("fasta")
    if not isinstance(text, str) or not text.strip():
        raise InputError("Paste or load a FASTA file with at least one protein.")
    proteins = batching.read_proteins(text, max_residues=MAX_RESIDUES)
    if not proteins:
        raise InputError("No protein records found.")
    if len(proteins) > MAX_BATCH:
        raise InputError(
            f"{len(proteins)} proteins is too many; the limit is {MAX_BATCH} at a time."
        )
    s = _settings(request)
    parts = (
        tuple(resolve_parts(_part_keys(request, "n_parts"))),
        tuple(resolve_parts(_part_keys(request, "c_parts"))),
    )
    result = batching.design_batch(
        proteins,
        s.vector,
        s.up,
        s.down,
        s.host,
        cterm=CTERMS[s.cterm_key],
        strategy=STRATEGIES[s.strategy_key],
        seed=s.seed,
        limits=s.limits,
        stop_codons=s.stop_codons,
        clamp_5=s.clamps[0],
        clamp_3=s.clamps[1],
        n_parts=parts[0],
        c_parts=parts[1],
    )
    rows = [{k: v for k, v in r.items() if k != "order_dna"} for r in batching.summary_rows(result)]
    return {
        "vector": s.vector.name,
        "upstream": s.up.label,
        "downstream": s.down.label,
        "host": s.host_slug,
        "n": len(result.entries),
        "n_ok": result.n_ok,
        "rows": rows,
        "csv": batching.summary_csv(result),
        "orders_fasta": batching.orders_fasta(result),
        "inserts_fasta": batching.inserts_fasta(result),
        "zip_b64": _zip_b64(result),
        "filename": f"{s.vector.name}_{s.up.enzyme}-{s.down.enzyme}_batch",
    }
