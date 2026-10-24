"""The DNA to order: the insert plus the extra bases a synthetic fragment needs.

An enzyme cuts poorly at the very end of a fragment, so a few bases (6 by default, NEB's general
guideline) go outside each site. They are in the order and the digest but not in the plasmid or
the reading frame. ``make_order`` picks bases that add no new site or host motif, then digests
the fragment and ligates it into the vector; the product must equal the planned plasmid.
"""

from __future__ import annotations

import hashlib
import random
import re
from dataclasses import dataclass, field

from Bio.Restriction import RestrictionBatch
from Bio.Seq import Seq
from Bio.SeqUtils import gc_fraction

from .assemble import Assembly
from .issues import Issue, Severity
from .optimize.automaton import Automaton, reverse_complement

DEFAULT_CLAMP = 6
MAX_CLAMP = 30
_TRIES = 600


@dataclass(frozen=True)
class LigationCheck:
    ok: bool
    message: str


@dataclass
class Order:
    dna: str  # clamp 5' + insert + clamp 3'
    clamp_5: str
    clamp_3: str
    insert_start: int  # where the insert begins inside ``dna``
    insert_end: int
    coding_span: tuple[int, int]  # coding region, in ``dna`` coordinates
    ligation: LigationCheck
    issues: list[Issue] = field(default_factory=list)

    @property
    def length(self) -> int:
        return len(self.dna)

    @property
    def gc(self) -> float:
        return gc_fraction(self.dna)

    @property
    def longest_homopolymer(self) -> int:
        return max(len(m.group()) for m in re.finditer(r"(.)\1*", self.dna))

    @property
    def ok(self) -> bool:
        return self.ligation.ok and not any(i.severity is Severity.ERROR for i in self.issues)


def _clamp_ok(clamp: str) -> bool:
    """Not a palindrome (it would pair with itself), no run of three, and not lopsided in GC."""
    if not clamp:
        return True
    if clamp == reverse_complement(clamp) or re.search(r"(.)\1\1", clamp):
        return False
    gc = sum(c in "GC" for c in clamp) / len(clamp)
    return len(clamp) < 6 or 0.34 <= gc <= 0.67


def _seed(insert: str, side: str) -> int:
    return int.from_bytes(hashlib.sha256(f"{side}|{insert}".encode()).digest()[:8], "big")


def _choose_clamp(n: int, insert: str, side: str, forbidden: list[str]) -> str:
    """The first acceptable ``n`` bases (from a seeded stream) that start no forbidden site.

    A match counts if it overlaps the clamp: the clamp must not create a recognition site or
    host motif together with the insert's end.  The insert's own sites do not overlap it.
    """
    if n == 0:
        return ""
    dfa = Automaton(forbidden)
    window = max(dfa.max_len, 1)
    rng = random.Random(_seed(insert, side))
    for _ in range(_TRIES):
        clamp = "".join(rng.choice("ACGT") for _ in range(n))
        if not _clamp_ok(clamp):
            continue
        text, lo, hi = (
            (clamp + insert[:window], 0, n)
            if side == "5"
            else (insert[-window:] + clamp, window, window + n)
        )
        _, matches = dfa.scan(text)
        if not any(m.start < hi and m.end > lo for m in matches):
            return clamp
    raise ValueError(
        f"could not find {n} extra bases that avoid every cloning site and host motif; "
        "try a different length"
    )


def _boundaries(seq: str, names: list[str], circular: bool) -> dict[str, list[int]]:
    """Cut positions as Biopython reports them, as 0-based top-strand boundaries."""
    found = RestrictionBatch(names).search(Seq(seq), linear=not circular)
    return {str(enzyme): sorted(p - 1 for p in cuts) for enzyme, cuts in found.items()}


def simulate_cloning(assembly: Assembly, order_dna: str) -> LigationCheck:
    """Digest ``order_dna`` and the vector, ligate, and compare with the planned plasmid.

    The cut positions come from Biopython, not from this package's own arithmetic, so a wrong
    cut offset or a mishandled clamp cannot agree with itself.
    """
    d = assembly.design
    up_name, down_name = d.site.upstream_enzyme, d.site.downstream_enzyme
    in_vector = _boundaries(d.vector.seq, [up_name, down_name], d.vector.circular)
    in_order = _boundaries(order_dna, [up_name, down_name], circular=False)

    for name in (up_name, down_name):
        n = len(in_order.get(name, []))
        if n != 1:
            return LigationCheck(
                False, f"{name} cuts the ordered fragment {n} times; it must cut exactly once"
            )
    cut_a = next((b for b in in_vector[up_name] if d.up.start <= b <= d.up.end), None)
    cut_b = next((b for b in in_vector[down_name] if d.down.start <= b <= d.down.end), None)
    if cut_a is None or cut_b is None:
        return LigationCheck(False, "the chosen sites are not cut in the vector")
    start, end = in_order[up_name][0], in_order[down_name][0]
    if not start < end:
        return LigationCheck(False, "the enzymes cut the ordered fragment in the wrong order")

    product = d.vector.seq[:cut_a] + order_dna[start:end] + d.vector.seq[cut_b:]
    if product != assembly.seq:
        return LigationCheck(
            False, "digesting and ligating the ordered fragment gives a different plasmid"
        )
    return LigationCheck(
        True,
        f"{up_name} and {down_name} digest of the order, ligated into the vector, gives the "
        f"planned {len(product):,} bp plasmid exactly",
    )


def make_order(
    assembly: Assembly, clamp_5: int = DEFAULT_CLAMP, clamp_3: int = DEFAULT_CLAMP
) -> Order:
    """Add ``clamp_5`` / ``clamp_3`` extra bases around the insert and verify the result."""
    for n in (clamp_5, clamp_3):
        if isinstance(n, bool) or not isinstance(n, int) or not 0 <= n <= MAX_CLAMP:
            raise ValueError(f"extra bases must be a whole number from 0 to {MAX_CLAMP}")
    d = assembly.design
    insert = d.insert_dna
    forbidden = [
        s
        for s in (
            d.vector.seq[d.up.start : d.up.end],
            d.vector.seq[d.down.start : d.down.end],
        )
    ] + list(d.host_motifs)
    c5 = _choose_clamp(clamp_5, insert, "5", forbidden)
    c3 = _choose_clamp(clamp_3, insert, "3", forbidden)
    dna = c5 + insert + c3

    coding = d.coding_start - d.up.start + len(c5)
    order = Order(
        dna=dna,
        clamp_5=c5,
        clamp_3=c3,
        insert_start=len(c5),
        insert_end=len(c5) + len(insert),
        coding_span=(coding, coding + len(d.coding)),
        ligation=simulate_cloning(assembly, dna),
    )
    if not order.ligation.ok:
        order.issues.append(Issue(Severity.ERROR, "order_does_not_clone", order.ligation.message))
    if not assembly.ok:
        order.issues.append(
            Issue(
                Severity.ERROR,
                "design_has_errors",
                "the design itself failed its checks, so this sequence should not be ordered",
            )
        )
    return order
