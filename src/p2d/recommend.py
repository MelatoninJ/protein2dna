"""Which two sites to use: every workable pair, designed for real, ranked by what it adds.

"Best" depends on what you want (keep the vector's N-terminal tag or lose it, end in a stop or read
into a C-terminal tag), so nothing is picked for you.  Every pair that can be cloned is designed
and checked, described by what it leaves on either end of your protein, and ordered by how little
it adds: first the residues put in front of and behind your protein, then the bases added only for
the reading frame, then whether the enzymes cut the vector once.  The order is a convenience, not a
verdict.
"""

from __future__ import annotations

from dataclasses import dataclass

from Bio.Seq import Seq

from .assemble import assemble
from .enzymes import CutSite, list_cut_sites, plan_cut
from .frame import CTerm, design_insert
from .host import HostProfile
from .model import Chain, ConstructPlan, Part, PartKind
from .optimize import CodonStrategy
from .optimize.refine import Limits
from .order import make_order
from .report import construct_report
from .vector import Vector

MAX_PAIRS = 600  # a window with more candidate pairs than this must be narrowed


@dataclass(frozen=True)
class PairOption:
    up: CutSite
    down: CutSite
    before: str  # residues the vector and the 5' site put in front of your protein
    after: str  # residues behind it (the 3' site and, if read, the C-terminal tag)
    n_tag: bool  # the vector's own N-terminal sequence stays on the protein
    pad_bases: int  # bases added only to keep the reading frame
    extra_residues: int
    altered: int  # residues of the vector's own protein this pair changes
    single_cutters: bool  # both enzymes cut the whole vector once
    replaced_bp: int
    plasmid_bp: int
    ligation_ok: bool

    @property
    def outcome(self) -> tuple:
        """What the protein ends up as; pairs with the same outcome are interchangeable."""
        return (self.before, self.after, self.n_tag, self.pad_bases)

    @property
    def sort_key(self) -> tuple:
        # least added to the protein first; then, among equal results, cleaner digests and the
        # smaller piece cut out of the vector
        return (
            self.altered,
            self.extra_residues,
            self.pad_bases,
            not self.single_cutters,
            self.replaced_bp,
            self.up.start,
            self.down.start,
        )


@dataclass(frozen=True)
class OutcomeGroup:
    """Pairs that give the same protein, best pair first."""

    before: str
    after: str
    n_tag: bool
    pad_bases: int
    pairs: tuple[PairOption, ...]

    @property
    def best(self) -> PairOption:
        return self.pairs[0]

    @property
    def altered(self) -> int:
        return self.best.altered


def _altered(vector: Vector, before: str) -> int:
    """How many residues in front of the protein differ from the vector's own reading.

    The vector's native protein (its tag, or a fusion partner such as a fluorescent protein) is
    read from its start codon.  A cut inside a codon can change a residue of it, for instance
    the last K of a tag, if the bases that complete the codon are not the vector's own.
    """
    start = vector.expression_start
    window = vector.seq[start : start + 3 * (len(before) + 1)]
    native = str(Seq(window[: len(window) - len(window) % 3]).translate(to_stop=True))
    return sum(1 for a, b in zip(before, native, strict=False) if a != b)


def _describe(vector, up, down, protein, host, cterm, strategy, n_parts, c_parts, limits):
    """Design one pair; returns ``(PairOption, None)`` or ``(None, reason)``."""
    try:
        site = plan_cut(vector, up, down)
        chain = Chain("payload", protein, n_terminal=list(n_parts), c_terminal=list(c_parts))
        plan = ConstructPlan([chain], host=host.slug)
        a = assemble(
            design_insert(plan, vector, site, host, strategy=strategy, cterm=cterm, limits=limits)
        )
    except ValueError as exc:
        return None, str(exc)
    if not a.ok:
        return None, a.errors[0].message
    order = make_order(a)
    if not order.ok:
        return None, order.ligation.message

    built = a.design.chain.insert_aa  # your protein with any tags and linkers you added
    at = a.protein.find(built)
    if at < 0:
        return None, "the protein is not found in the translation"
    before, after = a.protein[:at], a.protein[at + len(built) :]
    altered = _altered(vector, before)
    n_tag = bool(a.segments) and a.segments[0].kind is PartKind.TAG
    rep = construct_report(a)["sizes"]
    return (
        PairOption(
            up=up,
            down=down,
            before=before,
            after=after,
            n_tag=n_tag,
            pad_bases=len(a.design.pad_up) + len(a.design.pad_down),
            extra_residues=len(before) + len(after),
            altered=altered,
            single_cutters=up.count == 1 and down.count == 1,
            replaced_bp=rep["replaced_bp"],
            plasmid_bp=rep["plasmid_bp"],
            ligation_ok=order.ligation.ok,
        ),
        None,
    )


def rank_pairs(
    vector: Vector,
    protein: str,
    host: HostProfile,
    *,
    window: tuple[int, int] | None = None,
    cterm: CTerm = CTerm.STOP,
    n_tag: bool | None = None,
    strategy: CodonStrategy = CodonStrategy.MAX_CAI,
    limit: int | None = None,
    n_parts: tuple[Part, ...] = (),
    c_parts: tuple[Part, ...] = (),
    limits: Limits | None = None,
) -> tuple[list[PairOption], int]:
    """Every pair that clones ``protein`` into ``vector``, best first, and how many were tried.

    ``n_tag`` keeps only pairs that keep (True) or drop (False) the vector's N-terminal sequence;
    None keeps both.  ``window`` limits which sites are considered.  ``n_parts`` / ``c_parts``
    are tags and linkers built onto the protein; ``before`` / ``after`` still describe only what
    the vector and the sites add around that.
    """
    if vector.expression_start is None:
        raise ValueError("The start codon has not been confirmed.")
    sites = list_cut_sites(vector, window)
    candidates = [
        (up, down)
        for i, up in enumerate(sites)
        for down in sites[i + 1 :]
        if down.start > up.end and down.start != up.start
    ]
    if len(candidates) > MAX_PAIRS:
        raise ValueError(
            f"{len(candidates)} site pairs are too many to try; narrow the window to the "
            "region you want to clone into"
        )
    options: list[PairOption] = []
    for up, down in candidates:
        option, _ = _describe(
            vector, up, down, protein, host, cterm, strategy, n_parts, c_parts, limits
        )
        if option is not None and (n_tag is None or option.n_tag == n_tag):
            options.append(option)
    options.sort(key=lambda o: o.sort_key)
    return (options[:limit] if limit else options), len(candidates)


def group_by_outcome(options: list[PairOption]) -> list[OutcomeGroup]:
    """Collapse ranked pairs into the distinct results they give, in ranked order.

    Which 3' site you use rarely changes the protein when it ends in a stop codon, so a flat list
    would be dozens of near-identical rows.
    """
    groups: dict[tuple, list[PairOption]] = {}
    for option in options:  # already sorted, so the first member of each group is its best
        groups.setdefault(option.outcome, []).append(option)
    return [
        OutcomeGroup(before, after, n_tag, pads, tuple(members))
        for (before, after, n_tag, pads), members in groups.items()
    ]
