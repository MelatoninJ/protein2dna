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
SLACK = 6  # bp of tolerance around an annotated cloning region
MIN_SITE = 6  # shortest recognition site considered for cloning


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
        if enzyme.size < MIN_SITE:
            continue  # 4-5 bp sites cut a plasmid dozens of times and are not cloning enzymes
        yield enzyme


def candidate_enzymes(vector: Vector) -> list[EnzymeCandidate]:
    """Single cutters of ``vector`` that lie inside its expression ORF, by position.

    When the vector has a ``cloning_region`` only enzymes inside it (plus one that
    holds the start codon) are offered, so a pair can never cut a gene the vector
    needs.  Without one the whole expression ORF is searched, which is only safe for
    vectors whose ORF is short.  Isoschizomers are merged into one entry; the most
    widely sold one represents them.
    """
    window = expression_window(vector)
    region = vector.cloning_region
    groups: dict[tuple, list] = {}
    for enzyme in _usable_enzymes():
        hits = vector.find_all(str(enzyme))
        if len(hits) != 1:
            continue
        hit = hits[0]
        if hit.end > window.stop_start:
            continue
        # outside the cloning region only if it donates the start codon itself
        outside = bool(region) and not (
            region[0] - SLACK <= hit.start and hit.end <= region[1] + SLACK
        )
        if outside and _atg_offset(enzyme.site, hit.start, window) is None:
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


# ============================================================ site-level picking
#
# The name-based API above answers "which enzymes cut this vector once".  Real
# plasmids often have an enzyme with several sites, and the right one to use is a
# matter of *where*.  Everything below works on individual sites, identified by
# enzyme name and position, and asks the only question that matters for a complete
# digest: does a second site of either enzyme lie in the part of the plasmid that is
# kept?  A second site inside the piece being replaced is harmless.


@dataclass(frozen=True)
class CutSite:
    """One occurrence of a recognition site."""

    enzyme: str  # representative name; ``aliases`` cut the same site the same way
    aliases: tuple[str, ...]
    site: str
    start: int  # 0-based
    end: int  # exclusive
    overhang_kind: str  # "5'" or "3'"
    overhang: str
    cut_offset: int  # top-strand cut, in bases after the site start (EcoRI: 1, G^AATTC)
    count: int  # occurrences of this site in the whole vector
    number: int  # which occurrence this is, 1-based by position

    @property
    def id(self) -> str:
        return f"{self.enzyme}:{self.start}"

    @property
    def label(self) -> str:
        which = f" ({self.number} of {self.count})" if self.count > 1 else ""
        return f"{self.enzyme} at {self.start + 1}{which}"


def list_cut_sites(vector: Vector, window: tuple[int, int] | None = None) -> list[CutSite]:
    """Every cloning-grade site in ``vector`` (6 bp or longer, palindromic, sticky), by position.

    ``window`` keeps only sites lying wholly inside ``(start, end)``.  Enzymes with the same
    site and ends are merged under the most widely sold name.
    """
    site_hits: dict[str, list] = {}
    groups: dict[tuple, list] = {}
    for enzyme in _usable_enzymes():
        if enzyme.site not in site_hits:
            site_hits[enzyme.site] = vector.find_all(str(enzyme))
        if site_hits[enzyme.site]:
            key = (enzyme.site, enzyme.is_5overhang(), enzyme.ovhgseq, enzyme.fst5)
            groups.setdefault(key, []).append(enzyme)

    out = []
    for (site, five_prime, overhang, cut), members in groups.items():
        members.sort(key=lambda e: (-len(e.suppl), str(e)))
        hits = sorted(site_hits[site], key=lambda h: h.start)
        for number, hit in enumerate(hits, start=1):
            if window and not (window[0] <= hit.start and hit.end <= window[1]):
                continue
            out.append(
                CutSite(
                    enzyme=str(members[0]),
                    aliases=tuple(sorted(str(m) for m in members[1:])),
                    site=site,
                    start=hit.start,
                    end=hit.end,
                    overhang_kind="5'" if five_prime else "3'",
                    overhang=overhang,
                    cut_offset=cut,
                    count=len(hits),
                    number=number,
                )
            )
    return sorted(out, key=lambda c: (c.start, c.enzyme))


def find_cut_site(vector: Vector, site_id: str) -> CutSite:
    """Resolve an id such as ``"BamHI:150"`` (name or alias, 0-based start) to a site."""
    name, _, pos = site_id.rpartition(":")
    if not name or not pos.isdigit():
        raise ValueError(f"{site_id!r} is not a site id like 'EcoRI:168'")
    for c in list_cut_sites(vector):
        if c.start == int(pos) and name in (c.enzyme, *c.aliases):
            return c
    raise ValueError(f"No {name} site starts at position {int(pos) + 1} of {vector.name}")


def backbone_conflicts(vector: Vector, up: CutSite, down: CutSite) -> list[tuple[str, int]]:
    """Other sites of either enzyme that a complete digest would cut in the kept backbone.

    The piece between the two chosen sites is replaced, so extra sites in it do not matter;
    one outside it (anywhere else on the circle) would cut the backbone.
    """
    lo, hi = up.start, down.end
    bad = []
    for cut in (up, down):
        for hit in vector.find_all(cut.enzyme):
            chosen = hit.start in (up.start, down.start)
            if not chosen and not (lo <= hit.start and hit.end <= hi):
                bad.append((cut.enzyme, hit.start))
    return sorted(set(bad), key=lambda b: b[1])


def _start_for(vector: Vector, up: CutSite) -> tuple[dict, str]:
    """How translation reaches ``up``: from the vector's start codon, or from the ATG it holds."""
    start = vector.expression_start
    if start is None:
        raise ValueError("The start codon has not been confirmed.")
    if up.start >= start + 3:
        for i in range(start, up.start - 2, 3):
            if vector.seq[i : i + 3] in STOP_CODONS:
                raise ValueError(
                    f"A stop codon at {i + 1} lies between the start codon and {up.label}, "
                    "so the insert would never be translated"
                )
        return (
            dict(start_source=StartSource.VECTOR, vector_orf_start=start),
            "the vector's start codon",
        )
    offset = start - up.start
    if not (0 <= offset <= len(up.site) - 3 and up.site[offset : offset + 3] == "ATG"):
        raise ValueError(f"{up.label} lies before the start codon at {start + 1}")
    return (
        dict(start_source=StartSource.UPSTREAM_SITE, start_offset_in_site=offset),
        f"the ATG inside {up.enzyme}",
    )


def upstream_problem(vector: Vector, up: CutSite) -> str | None:
    """Why ``up`` cannot be the 5' site (before a partner is chosen), or None."""
    try:
        _start_for(vector, up)
    except ValueError as exc:
        return str(exc)
    return None


def plan_cut(vector: Vector, up: CutSite, down: CutSite) -> InsertionSite:
    """The ``InsertionSite`` for two chosen sites, or a ``ValueError`` saying why not."""
    start = vector.expression_start
    if start is None:
        raise ValueError(
            f"{vector.name} has no confirmed start codon, so the reading frame is unknown. "
            "Choose the start codon first."
        )
    if up.start >= down.start:
        raise ValueError(
            f"{up.label} must lie upstream of {down.label} for the insert to go in the "
            "right way round"
        )
    if down.start <= up.end:
        raise ValueError(
            f"{up.enzyme} and {down.enzyme} sites touch or overlap: there is no spacer "
            "between them, so a double digest cannot cut both"
        )
    if ends_compatible(up, down):
        raise ValueError(
            f"{up.enzyme} and {down.enzyme} leave compatible {up.overhang_kind} {up.overhang} "
            "ends: the insert could ligate in either orientation"
        )
    bad = backbone_conflicts(vector, up, down)
    if bad:
        where = ", ".join(f"{name} at {pos + 1}" for name, pos in bad)
        raise ValueError(
            f"A complete digest would also cut the backbone you keep ({where}). Choose "
            "sites so that every other site of these enzymes lies between them, or use "
            "a partial digest."
        )

    begin, origin = _start_for(vector, up)
    return InsertionSite(
        name=f"{up.enzyme}-{down.enzyme}",
        upstream_enzyme=up.enzyme,
        downstream_enzyme=down.enzyme,
        upstream_pos=up.start,
        downstream_pos=down.start,
        note=f"translation starts at {origin}",
        **begin,
    )
