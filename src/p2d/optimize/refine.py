"""Local search after the DP: pull GC and repeats inside limits without breaking the sites.

The DP picks codons one position at a time, so it cannot see properties of the whole
sequence.  This stage swaps single codons for synonymous ones, one at a time, whenever the
swap lowers a penalty for GC outside the allowed range or for repeated stretches.  It is
deterministic.  A swap is refused if it adds any occurrence of a forbidden pattern, so the
guarantee of the DP (no forbidden site on either strand, junctions included) is kept.

The penalty weights are internal and only decide which of several useful swaps goes first.
The limits themselves come from the host profile and can be set per design.
"""

from __future__ import annotations

import math
from collections import Counter
from dataclasses import dataclass, replace

from .automaton import Automaton, reverse_complement


@dataclass(frozen=True)
class Limits:
    gc_min: float = 0.35
    gc_max: float = 0.65
    window_size: int = 50
    window_min: float = 0.25
    window_max: float = 0.75
    direct_max: int = 20  # longest direct repeat allowed, in bases
    inverted_max: int = 20  # longest inverted repeat allowed

    def __post_init__(self) -> None:
        for lo, hi, what in (
            (self.gc_min, self.gc_max, "overall GC"),
            (self.window_min, self.window_max, "window GC"),
        ):
            if not 0 <= lo <= 1 or not 0 <= hi <= 1:
                raise ValueError(f"the {what} limits must be between 0% and 100%")
            if lo > hi:
                raise ValueError(f"the lower {what} limit is above the upper one")
        if not 10 <= self.window_size <= 500:
            raise ValueError("the GC window must be between 10 and 500 bases")
        if not 8 <= self.direct_max <= 200 or not 8 <= self.inverted_max <= 200:
            raise ValueError("repeat limits must be between 8 and 200 bases")

    @classmethod
    def from_profile(cls, raw: dict, **override) -> Limits:
        """Read the host YAML (``gc`` and ``repeats`` blocks); ``override`` wins."""
        gc, rep = raw.get("gc", {}), raw.get("repeats", {})
        win = gc.get("window", {})
        base = cls()
        values = {
            "gc_min": gc.get("global", [base.gc_min, base.gc_max])[0],
            "gc_max": gc.get("global", [base.gc_min, base.gc_max])[1],
            "window_size": win.get("size", base.window_size),
            "window_min": win.get("range", [base.window_min, base.window_max])[0],
            "window_max": win.get("range", [base.window_min, base.window_max])[1],
            "direct_max": rep.get("direct_max", base.direct_max),
            "inverted_max": rep.get("inverted_max", base.inverted_max),
        }
        values.update({k: v for k, v in override.items() if v is not None})
        return cls(**values)


@dataclass(frozen=True)
class Profile:
    """How a sequence stands against the limits."""

    gc: float
    window_low: float  # lowest and highest GC over any window
    window_high: float
    gc_excess: float  # percentage points outside the overall range
    window_excess: float  # percentage points outside the window range, summed over windows
    repeats: int  # bases that start a repeat longer than allowed (direct + inverted)

    @property
    def penalty(self) -> float:
        return 10 * self.gc_excess + self.window_excess + 5 * self.repeats

    @property
    def clean(self) -> bool:
        return self.penalty == 0


def measure(dna: str, limits: Limits) -> Profile:
    n = len(dna)
    gc = (dna.count("G") + dna.count("C")) / n if n else 0.0
    gc_excess = 100 * max(0.0, limits.gc_min - gc, gc - limits.gc_max) if n else 0.0

    low = high = gc
    window_excess = 0.0
    w = limits.window_size
    if n >= w:
        flags = [1 if c in "GC" else 0 for c in dna]
        total = sum(flags[:w])
        low = high = total / w
        window_excess += _outside(total / w, limits)
        for i in range(w, n):
            total += flags[i] - flags[i - w]
            share = total / w
            low, high = min(low, share), max(high, share)
            window_excess += _outside(share, limits)
    return Profile(gc, low, high, gc_excess, window_excess, _repeat_bases(dna, limits))


def _outside(share: float, limits: Limits) -> float:
    return 100 * max(0.0, limits.window_min - share, share - limits.window_max)


def _repeat_spans(dna: str, limits: Limits) -> list[tuple[int, int]]:
    """``(start, k)`` of every k-mer (k just over the limit) seen again elsewhere.

    A repeat longer than the limit contains such a k-mer, so this is empty exactly when no
    repeat exceeds it.  Copies must not overlap, or a homopolymer would repeat itself.
    """
    spans = []
    for limit, inverted in ((limits.direct_max, False), (limits.inverted_max, True)):
        k = limit + 1
        if len(dna) < 2 * k:
            continue
        seen: dict[str, list[int]] = {}
        for i in range(len(dna) - k + 1):
            seen.setdefault(dna[i : i + k], []).append(i)
        for i in range(len(dna) - k + 1):
            kmer = dna[i : i + k]
            other = seen.get(reverse_complement(kmer) if inverted else kmer, ())
            if any(abs(j - i) >= k for j in other):
                spans.append((i, k))
    return spans


def _repeat_bases(dna: str, limits: Limits) -> int:
    return len(_repeat_spans(dna, limits))


@dataclass
class Refinement:
    codons: list[str]
    before: Profile
    after: Profile
    swaps: int


def _touching(dfa: Automaton, text: str, lo: int, hi: int, margin: int) -> Counter:
    """Pattern matches in ``text`` that overlap ``[lo, hi)``, scanning only nearby."""
    a, b = max(0, lo - margin), min(len(text), hi + margin)
    _, found = dfa.scan(text[a:b])
    return Counter(m.pattern for m in found if a + m.start < hi and a + m.end > lo)


class _Gc:
    """GC counts per window, so a one-codon change is scored from the windows it touches."""

    def __init__(self, dna: str, limits: Limits) -> None:
        self.limits, self.w, self.n = limits, limits.window_size, len(dna)
        self.flags = [1 if c in "GC" else 0 for c in dna]
        self._rebuild()

    def _rebuild(self) -> None:
        self.total = sum(self.flags)
        n_windows = max(0, self.n - self.w + 1)
        self.sums = [0] * n_windows
        if n_windows:
            run = sum(self.flags[: self.w])
            self.sums[0] = run
            for s in range(1, n_windows):
                run += self.flags[s + self.w - 1] - self.flags[s - 1]
                self.sums[s] = run

    def _global(self, total: int) -> float:
        share = total / self.n
        return 100 * max(0.0, self.limits.gc_min - share, share - self.limits.gc_max)

    def penalty_gain(self, i: int, codon: str) -> float:
        """How much GC penalty replacing codon ``i`` by ``codon`` removes (positive = better)."""
        p, w = 3 * i, self.w
        delta = [(1 if c in "GC" else 0) - self.flags[p + k] for k, c in enumerate(codon)]
        if not any(delta):
            return 0.0
        gain = 10 * (self._global(self.total) - self._global(self.total + sum(delta)))
        for s in range(max(0, p - w + 1), min(p + 3, len(self.sums))):
            d = sum(delta[k] for k in range(3) if s <= p + k < s + w)
            if d:
                gain += _outside(self.sums[s] / w, self.limits) - _outside(
                    (self.sums[s] + d) / w, self.limits
                )
        return gain

    def apply(self, i: int, codon: str) -> None:
        p = 3 * i
        for k, c in enumerate(codon):
            self.flags[p + k] = 1 if c in "GC" else 0
        self._rebuild()

    def worst_codons(self) -> set[int]:
        """Codons inside a window that is out of range (all of them if overall GC is)."""
        if self._global(self.total):
            return set(range(self.n // 3))
        out: set[int] = set()
        for s, v in enumerate(self.sums):
            if _outside(v / self.w, self.limits):
                out.update(range(s // 3, (s + self.w - 1) // 3 + 1))
        return out


def _problem_codons(dna: str, limits: Limits) -> set[int]:
    """Codons that overlap a k-mer which repeats elsewhere."""
    out: set[int] = set()
    for i, k in _repeat_spans(dna, limits):
        out.update(range(i // 3, (i + k - 1) // 3 + 1))
    return out


MAX_REPEAT_TRIES = 600  # full re-measurements spent on repeats, so a hopeless one ends quickly


def refine(
    aa: str,
    codons: list[str],
    usage: dict[str, dict[str, float]],
    prefix: str,
    suffix: str,
    patterns: list[str],
    limits: Limits,
) -> Refinement:
    """Swap synonymous codons to bring ``codons`` inside ``limits``; sites stay out.

    GC first: each round ranks the useful single-codon swaps and applies them best-first,
    re-scoring each against what has already changed.  Then repeats, one swap at a time,
    never at the cost of GC.
    A swap is taken only if it adds no occurrence of any pattern near the codon.
    """
    codons = list(codons)
    before = measure("".join(codons), limits)
    dfa = Automaton(patterns)
    margin = max(dfa.max_len, 1)
    swaps = 0

    def safe(i: int, alt: str) -> bool:
        dna = "".join(codons)
        trial = dna[: 3 * i] + alt + dna[3 * i + 3 :]
        lo = len(prefix) + 3 * i
        old = _touching(dfa, prefix + dna + suffix, lo, lo + 3, margin)
        new = _touching(dfa, prefix + trial + suffix, lo, lo + 3, margin)
        return not any(new[p] > old[p] for p in new)

    def cost(i: int, alt: str) -> float:
        return math.log(usage[aa[i]][codons[i]] / usage[aa[i]][alt])

    for _ in range(60):
        gc = _Gc("".join(codons), limits)
        ranked = []
        for i in sorted(gc.worst_codons()):
            for alt in usage[aa[i]]:
                if alt != codons[i] and (gain := gc.penalty_gain(i, alt)) > 1e-9:
                    ranked.append((-gain, cost(i, alt), i, alt))
        if not ranked:
            break
        ranked.sort()
        moved = 0
        for _, _, i, alt in ranked:
            if codons[i] == alt or alt not in usage[aa[i]]:
                continue
            if gc.penalty_gain(i, alt) > 1e-9 and safe(i, alt):
                gc.apply(i, alt)
                codons[i] = alt
                moved += 1
        swaps += moved
        if not moved:
            break

    current = measure("".join(codons), limits)
    tries = 0
    while current.repeats and tries < MAX_REPEAT_TRIES:
        dna = "".join(codons)
        options = []
        for i in sorted(_problem_codons(dna, limits)):
            for alt in usage[aa[i]]:
                if alt != codons[i]:
                    options.append((cost(i, alt), i, alt))
        options.sort()
        chosen = None
        for _, i, alt in options:
            if tries >= MAX_REPEAT_TRIES:
                break
            tries += 1
            trial = measure(dna[: 3 * i] + alt + dna[3 * i + 3 :], limits)
            # GC comes first: a repeat is only fixed if that does not push GC further out
            better = (
                trial.repeats < current.repeats
                and trial.gc_excess <= current.gc_excess
                and trial.window_excess <= current.window_excess
            )
            if better and safe(i, alt):
                chosen = (i, alt, trial)
                break
        if chosen is None:
            break
        i, alt, current = chosen
        codons[i] = alt
        swaps += 1
    return Refinement(codons, before, measure("".join(codons), limits), swaps)


def with_limits(limits: Limits, **changes) -> Limits:
    return replace(limits, **changes)
