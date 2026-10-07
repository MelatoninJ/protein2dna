"""Codon choice as a shortest path through (residue, automaton state).

Forbidden sites are an Aho-Corasick DFA; the DP carries the DFA state across
codon choices, so the result is provably free of every forbidden pattern
(either strand) and optimal for the per-codon objective.  When no path exists
that is a definitive answer -- the site is unavoidable by synonymous changes --
and is reported as such instead of being papered over.
"""

from __future__ import annotations

import hashlib
import math
import random
from dataclasses import dataclass, field
from enum import Enum

from ..host import HostProfile
from ..issues import Issue, Severity
from .automaton import Automaton

DEFAULT_USAGE_FLOOR = 0.10


class CodonStrategy(str, Enum):
    """What the DP minimises.  ``MAX_CAI`` is implemented but is not the default:
    always taking the best codon depletes charged tRNA and can jam ribosomes."""

    MAX_CAI = "max_cai"
    WEIGHTED_SAMPLE = "weighted_sample"
    HARMONIZED = "harmonized"


@dataclass
class OptimizationRequest:
    aa: str
    host: HostProfile
    prefix: str  # fixed DNA immediately 5' of the optimised region
    suffix: str  # fixed DNA immediately 3'
    forbidden: list[str]  # recognition sites, before RC expansion
    strategy: CodonStrategy = CodonStrategy.WEIGHTED_SAMPLE
    seed: int | None = None
    usage_floor: float = DEFAULT_USAGE_FLOOR  # WEIGHTED_SAMPLE only


@dataclass
class OptimizationResult:
    dna: str
    issues: list[Issue] = field(default_factory=list)
    metrics: dict = field(default_factory=dict)  # at minimum: cai, gc, n_rare

    @property
    def ok(self) -> bool:
        return not any(i.severity is Severity.ERROR for i in self.issues)


@dataclass
class _Failure:
    """Where and why the DP ran out of paths."""

    index: int  # residue index at which every state died; len(aa) means the 3' junction
    patterns: list[str]


# ------------------------------------------------------------------ weights


def _default_seed(req: OptimizationRequest) -> int:
    """Reproducible by default: the same design request always gives the same DNA."""
    key = "|".join([req.aa, req.prefix, req.suffix, req.host.slug, *sorted(req.forbidden)])
    return int.from_bytes(hashlib.sha256(key.encode()).digest()[:8], "big")


def _gumbel_noise(req: OptimizationRequest) -> list[dict[str, float]]:
    """One draw per (position, codon), made up-front from the seeded RNG.

    ``argmax(log f + Gumbel)`` samples codons in proportion to ``f``, so the DP
    over these perturbed scores is a constraint-respecting weighted sample.
    """
    seed = req.seed if req.seed is not None else _default_seed(req)
    rng = random.Random(seed)
    noise = []
    for residue in req.aa:
        noise.append(
            {c: -math.log(-math.log(1.0 - rng.random())) for c in sorted(_codons(req, residue))}
        )
    return noise


def _codons(req: OptimizationRequest, residue: str) -> dict[str, float]:
    try:
        return req.host.codons_for(residue)
    except KeyError as exc:
        raise ValueError(f"no codons for residue {residue!r} in host {req.host.slug!r}") from exc


def _options(
    req: OptimizationRequest, noise: list[dict[str, float]] | None, floor: float
) -> list[list[tuple[str, float]]]:
    """Per position: ``(codon, cost)`` pairs, cheapest-first not required."""
    out = []
    for i, residue in enumerate(req.aa):
        usage = _codons(req, residue)
        allowed = {c: f for c, f in usage.items() if f >= floor} or {
            max(usage, key=usage.get): max(usage.values())
        }
        if noise is None:
            out.append([(c, -math.log(f)) for c, f in sorted(allowed.items())])
        else:
            out.append([(c, -(math.log(f) + noise[i][c])) for c, f in sorted(allowed.items())])
    return out


# ---------------------------------------------------------------------- DP


def _solve(
    options: list[list[tuple[str, float]]], dfa: Automaton, prefix: str, suffix: str
) -> tuple[list[str] | None, _Failure | None]:
    """Min-cost codon path that never completes a forbidden match in the region.

    Matches wholly inside ``prefix`` or ``suffix`` are the caller's business and
    are ignored; a match counts only if it overlaps the optimised region.
    """
    suffix = suffix.upper()
    memo: dict[tuple[int, str], tuple[int, str | None]] = {}

    def advance(q: int, codon: str) -> tuple[int, str | None]:
        """Feed one codon; second item names the pattern that blocked it, if any."""
        key = (q, codon)
        if key not in memo:
            blocked = None
            state = q
            for base in codon:
                state = dfa.step(state, base)
                hit = dfa.match_at(state)
                if hit:
                    blocked = hit[1]
                    break
            memo[key] = (-1, blocked) if blocked else (state, None)
        return memo[key]

    tail_memo: dict[int, str | None] = {}

    def tail_blocker(q: int) -> str | None:
        """Pattern the fixed 3' context completes across the region boundary, if any."""
        if q not in tail_memo:
            blocker = None
            state = q
            for k, base in enumerate(suffix, start=1):
                state = dfa.step(state, base)
                hit = dfa.match_at(state)
                # longer than what we have consumed => it started inside the region
                if hit and hit[0] > k:
                    blocker = hit[1]
                    break
            tail_memo[q] = blocker
        return tail_memo[q]

    layer: dict[int, float] = {dfa.feed(prefix): 0.0}
    back: list[dict[int, tuple[int, str]]] = []
    for i, opts in enumerate(options):
        nxt_cost: dict[int, float] = {}
        nxt_back: dict[int, tuple[int, str]] = {}
        blocked_by: set[str] = set()
        for q, base_cost in layer.items():
            for codon, w in opts:
                q2, blocker = advance(q, codon)
                if blocker:
                    blocked_by.add(blocker)
                    continue
                c = base_cost + w
                if q2 not in nxt_cost or c < nxt_cost[q2]:
                    nxt_cost[q2] = c
                    nxt_back[q2] = (q, codon)
        if not nxt_cost:
            return None, _Failure(i, sorted(blocked_by))
        back.append(nxt_back)
        layer = nxt_cost

    survivors = {q: c for q, c in layer.items() if tail_blocker(q) is None}
    if not survivors:
        return None, _Failure(len(options), sorted({b for q in layer if (b := tail_blocker(q))}))

    q = min(survivors, key=lambda s: (survivors[s], s))
    codons: list[str] = []
    for nxt_back in reversed(back):
        q, codon = nxt_back[q]
        codons.append(codon)
    return codons[::-1], None


# ------------------------------------------------------------------ metrics


def _metrics(req: OptimizationRequest, codons: list[str]) -> dict:
    host = req.host
    logs = []
    for residue, codon in zip(req.aa, codons, strict=True):
        usage = host.codons_for(residue)
        logs.append(math.log(usage[codon] / max(usage.values())))
    dna = "".join(codons)
    return {
        "cai": math.exp(sum(logs) / len(logs)) if logs else 1.0,
        "gc": (dna.count("G") + dna.count("C")) / len(dna) if dna else 0.0,
        "n_rare": sum(1 for c in codons if c in host.rare_codons),
        "strategy": req.strategy.value,
        "n_codons": len(codons),
    }


# ------------------------------------------------------------------- public


def optimize_coding(req: OptimizationRequest) -> OptimizationResult:
    """Choose synonymous codons for ``req.aa`` so that no forbidden site appears
    anywhere in ``prefix + dna + suffix`` that overlaps ``dna``.

    On failure ``dna`` is the naive best-codon sequence and the result carries an
    ``unavoidable_site`` error; callers must check ``result.ok`` before using it.
    """
    if req.strategy is CodonStrategy.HARMONIZED:
        raise NotImplementedError(
            "HARMONIZED needs a source-organism usage table; planned after v0.2"
        )

    aa = req.aa.upper()
    if aa != req.aa:
        req = OptimizationRequest(**{**req.__dict__, "aa": aa})
    dfa = Automaton(req.forbidden)

    sampling = req.strategy is CodonStrategy.WEIGHTED_SAMPLE
    noise = _gumbel_noise(req) if sampling else None
    issues: list[Issue] = []

    # The floor keeps sampling away from rare codons.  If it leaves no path, drop
    # it before concluding a site is unavoidable: that is a statement about the
    # protein, not about codon preference.
    floors = [req.usage_floor, 0.0] if sampling and req.usage_floor > 0 else [0.0]
    codons: list[str] | None = None
    failure: _Failure | None = None
    for floor in floors:
        options = _options(req, noise, floor)
        codons, failure = _solve(options, dfa, req.prefix, req.suffix)
        if codons is not None:
            if floor != floors[0]:
                issues.append(
                    Issue(
                        Severity.INFO,
                        "usage_floor_relaxed",
                        f"no site-free path used only codons with usage >= {floors[0]:g}; "
                        "rarer codons were allowed",
                    )
                )
            break

    if codons is None:
        assert failure is not None
        names = ", ".join(failure.patterns) or "a forbidden site"
        where = (
            "at the 3' junction with the fixed downstream sequence"
            if failure.index >= len(aa)
            else f"at residue {failure.index + 1} ({aa[failure.index]}) of the optimised region"
        )
        issues.append(
            Issue(
                Severity.ERROR,
                "unavoidable_site",
                f"{names} cannot be avoided {where}: every synonymous codon choice creates it. "
                "Pick a different enzyme pair, or change the protein sequence there.",
            )
        )
        naive = [req.host.best_codon(r) for r in aa]
        return OptimizationResult(
            dna="".join(naive), issues=issues, metrics=_metrics(req, naive) | {"feasible": False}
        )

    # What the same weights would pick with no constraints; the gap is the price
    # of the forbidden sites and is worth telling the user about.
    free = [min(opts, key=lambda o: (o[1], o[0]))[0] for opts in _options(req, noise, floor)]
    n_forced = sum(1 for a, b in zip(free, codons, strict=True) if a != b)
    if n_forced:
        issues.append(
            Issue(
                Severity.INFO,
                "constraint_swaps",
                f"{n_forced} codon(s) differ from the unconstrained choice "
                "to avoid forbidden sites",
            )
        )

    return OptimizationResult(
        dna="".join(codons), issues=issues, metrics=_metrics(req, codons) | {"feasible": True}
    )
