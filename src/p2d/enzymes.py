"""Which restriction enzymes can clone into a vector, and which pairs work together.

Everything here is derived from the vector sequence, never from its annotation
labels (see ``vector.py``).  The one curated input is ``Vector.expression_start``:
the ATG the ribosome actually uses.  A GenBank file cannot tell you that
reliably, and every frame decision below hangs on it.

Scope for now: palindromic Type II enzymes that cut inside their site and leave
a sticky end.  That covers the enzymes people clone with; Type IIS (Golden Gate)
arrives with v1.0 and blunt cutters are left out because they cannot give a
directional insert.
"""

from __future__ import annotations

from dataclasses import dataclass, field

from Bio.Restriction import CommOnly

from .vector import STOP_CODONS, InsertionSite, StartSource, Vector

_COMPLEMENT = str.maketrans("ACGT", "TGCA")


@dataclass(frozen=True)
class Window:
    """The expression ORF the vector gives us to clone into (0-based, stop included)."""

    start: int
    stop_end: int

    @property
    def stop_start(self) -> int:
        return self.stop_end - 3


@dataclass(frozen=True)
class EnzymeCandidate:
    name: str  # representative name; ``aliases`` cut the same site the same way
    site: str
    start: int
    end: int
    overhang_kind: str  # "5'" or "3'"
    overhang: str
    aliases: tuple[str, ...] = field(default_factory=tuple)

    @property
    def label(self) -> str:
        return f"{self.name} ({self.overhang_kind} {self.overhang})"


def expression_window(vector: Vector) -> Window:
    """The in-frame ORF from the curated start codon to its first stop."""
    start = vector.expression_start
    if start is None:
        raise ValueError(
            f"{vector.name} has no curated expression start codon, so p2d cannot say "
            "which reading frame to clone into. Set Vector.expression_start."
        )
    if vector.seq[start : start + 3] != "ATG":
        raise ValueError(f"expression_start {start} of {vector.name} is not an ATG")
    for i in range(start, len(vector.seq) - 2, 3):
        if vector.seq[i : i + 3] in STOP_CODONS:
            return Window(start, i + 3)
    raise ValueError(f"the expression ORF of {vector.name} never reaches a stop codon")


def _reverse_complement(seq: str) -> str:
    return seq.translate(_COMPLEMENT)[::-1]


def _usable_enzymes():
    """Commercial, palindromic, ACGT-only sticky cutters that cut inside the site."""
    for enzyme in CommOnly:
        site = enzyme.site
        if not set(site) <= set("ACGT") or enzyme.is_blunt() or not enzyme.is_palindromic():
            continue
        if not 0 < enzyme.fst5 < enzyme.size or "N" in enzyme.ovhgseq:
            continue
        yield enzyme


def candidate_enzymes(vector: Vector) -> list[EnzymeCandidate]:
    """Single cutters of ``vector`` that lie inside its expression ORF, by position.

    Enzymes with the same site and the same ends (isoschizomers) are merged into
    one entry; the most widely sold one represents them.
    """
    window = expression_window(vector)
    groups: dict[tuple, list] = {}
    for enzyme in _usable_enzymes():
        hits = vector.find_all(str(enzyme))
        if len(hits) != 1:
            continue
        hit = hits[0]
        if hit.end > window.stop_start:
            continue
        # a site must sit wholly after the start codon, unless it donates that ATG itself
        if hit.start < window.start + 3 and _atg_offset(enzyme.site, hit.start, window) is None:
            continue
        key = (enzyme.site, enzyme.is_5overhang(), enzyme.ovhgseq, hit.start)
        groups.setdefault(key, []).append((enzyme, hit))

    out = []
    for (site, five_prime, ovhg, _), members in groups.items():
        members.sort(key=lambda m: (-len(m[0].suppl), str(m[0])))
        enzyme, hit = members[0]
        out.append(
            EnzymeCandidate(
                name=str(enzyme),
                site=site,
                start=hit.start,
                end=hit.end,
                overhang_kind="5'" if five_prime else "3'",
                overhang=ovhg,
                aliases=tuple(sorted(str(m[0]) for m in members[1:])),
            )
        )
    return sorted(out, key=lambda c: (c.start, c.name))


def _atg_offset(site: str, site_start: int, window: Window) -> int | None:
    """Index of the expression ATG inside ``site`` if the site donates it, else None."""
    idx = window.start - site_start
    return idx if 0 <= idx <= len(site) - 3 and site[idx : idx + 3] == "ATG" else None


def ends_compatible(a: EnzymeCandidate, b: EnzymeCandidate) -> bool:
    """True if the two enzymes leave ends that ligate to each other.

    Compatible ends make directional cloning impossible: the insert can go in
    either way round, or the cut vector can close on itself.
    """
    if a.overhang_kind != b.overhang_kind:
        return False
    return a.overhang == b.overhang or a.overhang == _reverse_complement(b.overhang)


def plan_site(vector: Vector, upstream: str, downstream: str) -> InsertionSite:
    """Build the ``InsertionSite`` for an arbitrary enzyme pair, or say why it cannot work."""
    window = expression_window(vector)
    by_name = {}
    for c in candidate_enzymes(vector):
        by_name[c.name] = c
        for alias in c.aliases:
            by_name[alias] = c
    try:
        up, down = by_name[upstream], by_name[downstream]
    except KeyError as exc:
        raise ValueError(
            f"{exc.args[0]} is not a single-cutter inside the expression ORF of {vector.name}"
        ) from None
    if up.start >= down.start:
        raise ValueError(
            f"{up.name} (position {up.start}) must lie upstream of {down.name} "
            f"(position {down.start}) for the insert to go in the right way round"
        )
    if down.start <= up.end:
        raise ValueError(
            f"{up.name} and {down.name} sites touch or overlap: there is no spacer "
            "between them, so a double digest cannot cut both"
        )
    if ends_compatible(up, down):
        raise ValueError(
            f"{up.name} and {down.name} leave compatible {up.overhang_kind} {up.overhang} ends: "
            "the insert could ligate in either orientation"
        )

    offset = _atg_offset(up.site, up.start, window)
    if up.start >= window.start + 3:
        start = dict(start_source=StartSource.VECTOR, vector_orf_start=window.start)
        origin = "the vector's start codon"
    else:
        start = dict(start_source=StartSource.UPSTREAM_SITE, start_offset_in_site=offset)
        origin = f"the ATG inside {up.name}"
    return InsertionSite(
        name=f"{up.name}-{down.name}",
        upstream_enzyme=up.name,
        downstream_enzyme=down.name,
        note=f"translation starts at {origin}",
        **start,
    )
