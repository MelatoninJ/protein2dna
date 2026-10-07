"""Regex-based reference for forbidden-site checks.

Deliberately shares no code with ``p2d.optimize.automaton``, so a bug there
cannot hide itself from the tests that use this.
"""

import re

from p2d.vector import site_regex

_IUPAC_COMPLEMENT = str.maketrans("ACGTRYSWKMBDHVN", "TGCAYRSWMKVHDBN")


def reverse_complement_iupac(pattern: str) -> str:
    return pattern.upper().translate(_IUPAC_COMPLEMENT)[::-1]


def hits(text: str, patterns: list[str]) -> set[tuple[int, int]]:
    """``(start, end)`` of every overlapping match of any pattern on either strand."""
    found: set[tuple[int, int]] = set()
    for pat in patterns:
        for rx in {site_regex(pat), site_regex(reverse_complement_iupac(pat))}:
            for m in re.finditer(f"(?=({rx.pattern}))", text):
                found.add((m.start(), m.start() + len(m.group(1))))
    return found


def region_hits(prefix: str, dna: str, suffix: str, patterns: list[str]) -> set[tuple[int, int]]:
    """Matches in ``prefix + dna + suffix`` that overlap ``dna`` (the optimised region)."""
    lo, hi = len(prefix), len(prefix) + len(dna)
    return {(s, e) for s, e in hits(prefix + dna + suffix, patterns) if s < hi and e > lo}
