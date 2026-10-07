"""Aho-Corasick automaton over forbidden DNA patterns.

Patterns are IUPAC strings.  They are expanded to concrete ACGT strings and
their reverse complements are added, because a Type IIS site such as BsaI
(``GGTCTC``) is just as real on the opposite strand (``GAGACC``).

The automaton is built as a complete DFA: ``step(state, base)`` is one table
lookup with no failure-link chasing, which is what the DP in ``core.py`` needs
to run in ``O(n * |Q| * 6)``.
"""

from __future__ import annotations

from collections import deque
from collections.abc import Iterable
from dataclasses import dataclass
from itertools import product

from ..vector import IUPAC

ALPHABET = "ACGT"
_INDEX = {b: i for i, b in enumerate(ALPHABET)}
_COMPLEMENT = str.maketrans("ACGT", "TGCA")

# An ambiguity-heavy pattern (e.g. a long run of N) expands to 4**k strings.
# Real recognition sites are far below this; anything above it is a mistake.
MAX_EXPANSION = 4096


def reverse_complement(seq: str) -> str:
    return seq.translate(_COMPLEMENT)[::-1]


def _bases(code: str) -> str:
    try:
        cls = IUPAC[code]
    except KeyError as exc:
        raise ValueError(f"invalid IUPAC code {code!r} in forbidden pattern") from exc
    return cls.strip("[]")


def expand_pattern(pattern: str) -> list[str]:
    """Concrete ACGT strings matched by an IUPAC ``pattern``."""
    pattern = pattern.upper()
    if not pattern:
        raise ValueError("empty forbidden pattern")
    choices = [_bases(c) for c in pattern]
    size = 1
    for c in choices:
        size *= len(c)
        if size > MAX_EXPANSION:
            raise ValueError(
                f"forbidden pattern {pattern!r} expands to more than {MAX_EXPANSION} "
                "sequences; use a less ambiguous pattern"
            )
    return ["".join(p) for p in product(*choices)]


@dataclass(frozen=True)
class Match:
    end: int  # index one past the last base of the match, in the scanned text
    length: int
    pattern: str  # the caller's original pattern, not the expanded string

    @property
    def start(self) -> int:
        return self.end - self.length


class Automaton:
    """Multi-pattern matcher over ``ACGT``, both strands, IUPAC-aware.

    Any base outside ``ACGT`` (``N`` in a vector, say) resets to the root: it
    cannot be part of a concrete site, so no match may span it.
    """

    def __init__(self, patterns: Iterable[str]) -> None:
        sources: dict[str, str] = {}
        for original in patterns:
            for concrete in expand_pattern(original):
                sources.setdefault(concrete, original)
                sources.setdefault(reverse_complement(concrete), original)
        self.patterns = sorted(sources)
        self.max_len = max((len(p) for p in self.patterns), default=0)

        children: list[list[int]] = [[-1] * 4]
        depth = [0]
        hit_len = [0]
        hit_src: list[str | None] = [None]
        for concrete in self.patterns:
            node = 0
            for base in concrete:
                i = _INDEX[base]
                if children[node][i] == -1:
                    children.append([-1] * 4)
                    depth.append(depth[node] + 1)
                    hit_len.append(0)
                    hit_src.append(None)
                    children[node][i] = len(children) - 1
                node = children[node][i]
            hit_len[node] = len(concrete)
            hit_src[node] = sources[concrete]

        fail = [0] * len(children)
        queue: deque[int] = deque()
        for i in range(4):
            child = children[0][i]
            if child == -1:
                children[0][i] = 0
            else:
                queue.append(child)
        while queue:
            node = queue.popleft()
            # a node also ends every pattern that is a suffix of it: inherit the
            # longest from the failure state so a lookup is O(1)
            if hit_len[fail[node]] > hit_len[node]:
                hit_len[node] = hit_len[fail[node]]
                hit_src[node] = hit_src[fail[node]]
            for i in range(4):
                child = children[node][i]
                if child == -1:
                    children[node][i] = children[fail[node]][i]
                else:
                    fail[child] = children[fail[node]][i]
                    queue.append(child)

        self._next = children
        self._hit_len = hit_len
        self._hit_src = hit_src

    @property
    def n_states(self) -> int:
        return len(self._next)

    def step(self, state: int, base: str) -> int:
        i = _INDEX.get(base)
        return 0 if i is None else self._next[state][i]

    def match_at(self, state: int) -> tuple[int, str] | None:
        """Longest pattern ending at ``state`` as ``(length, original_pattern)``."""
        length = self._hit_len[state]
        return (length, self._hit_src[state]) if length else None  # type: ignore[arg-type]

    def feed(self, text: str, state: int = 0) -> int:
        """Advance through ``text``, ignoring matches (used for fixed context)."""
        for base in text.upper():
            state = self.step(state, base)
        return state

    def scan(self, text: str, state: int = 0) -> tuple[int, list[Match]]:
        """Advance through ``text`` and report the longest match ending at each base."""
        matches: list[Match] = []
        for i, base in enumerate(text.upper(), start=1):
            state = self.step(state, base)
            hit = self.match_at(state)
            if hit:
                matches.append(Match(end=i, length=hit[0], pattern=hit[1]))
        return state, matches

    def find_all(self, text: str) -> list[Match]:
        return self.scan(text)[1]
