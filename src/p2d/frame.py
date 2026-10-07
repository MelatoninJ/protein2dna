"""The frame engine.

Everything here exists to answer one question precisely: *given this vector,
these two enzymes and this payload, what protein actually comes out?*

The arithmetic that matters:

    final = vector[:up.start] + insert + vector[down.end:]
    insert = up_site + pad_up + coding + pad_down + down_site

``pad_up`` makes the payload's first codon land in frame with the start codon.
``pad_down`` makes whatever the vector encodes downstream land in frame.  Both
pads are chosen by a small search that refuses bases creating a stop codon or a
new copy of either enzyme's site.  The coding region itself is chosen by
``p2d.optimize`` *given* those pads and the surrounding vector sequence, so a
site straddling a junction is avoided rather than discovered afterwards.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from enum import Enum

from Bio.Seq import Seq

from .host import HostProfile
from .issues import Issue, Severity
from .model import Chain, ConstructPlan, PartKind, PartSource, StopPolicy
from .optimize import Automaton, CodonStrategy, OptimizationRequest, optimize_coding
from .vector import EnzymeHit, InsertionSite, StartSource, Vector, enzyme_site, site_regex

PAD_BASES = ["GC", "GG", "CG", "GA", "AC", "CA", "TC", "CT", "AG", "GT"]
STOP_CODONS = {"TAA", "TAG", "TGA"}
DOWNSTREAM_SCAN_NT = 600


@dataclass
class Segment:
    """A contiguous run of codons with one provenance, for the fusion preview."""

    label: str
    aa: str
    kind: PartKind
    source: PartSource

    def __len__(self) -> int:
        return len(self.aa)


@dataclass
class InsertDesign:
    vector: Vector
    site: InsertionSite
    chain: Chain
    up: EnzymeHit
    down: EnzymeHit
    orf_start: int
    pad_up: str
    pad_down: str
    coding: str
    insert_dna: str
    coding_start: int  # index of payload-block first base in final coords
    issues: list[Issue] = field(default_factory=list)
    stop_in_insert: bool = False
    metrics: dict = field(default_factory=dict)  # codon metrics from the optimiser

    @property
    def errors(self) -> list[Issue]:
        return [i for i in self.issues if i.severity is Severity.ERROR]

    @property
    def ok(self) -> bool:
        return not self.errors


# --------------------------------------------------------------------- helpers


def reverse_translate(aa: str, host: HostProfile) -> str:
    """Naive best-codon baseline.

    No longer used by ``design_insert`` (``p2d.optimize`` replaced it in v0.2).
    Kept as the reference the optimiser is tested against, and for isolating
    whether a failure is in the optimiser or in the frame engine.
    """
    return "".join(host.best_codon(r) for r in aa)


def _pad_candidates(
    n: int,
    prefix: str,
    suffix: str,
    forbidden: list[str],
    frame_offset: int,
) -> list[str]:
    """All ``n``-base pads that introduce neither a stop codon nor a new site.

    ``frame_offset`` is the position of ``prefix``'s first base within its
    codon, so stop codons can be read in the correct frame.  Sites running into
    the coding region are not judged here -- the optimiser sees the pad as fixed
    context and avoids them -- so more than one candidate is returned and the
    caller tries them in order.
    """
    if n == 0:
        return [""]
    candidates = [p[:n] for p in PAD_BASES] if n <= 2 else ["GCG", "GGC"]
    pad_lo, pad_hi = len(prefix), len(prefix) + n
    usable = []
    for pad in dict.fromkeys(candidates):
        window = prefix + pad + suffix
        codons = [window[i : i + 3] for i in range(-frame_offset % 3, len(window) - 2, 3)]
        if any(c in STOP_CODONS for c in codons):
            continue
        # a forbidden site already present in prefix/suffix is intentional;
        # only reject matches the pad itself creates
        if any(
            m.start() < pad_hi and m.end() > pad_lo
            for site in forbidden
            for m in site_regex(site).finditer(window)
        ):
            continue
        usable.append(pad)
    if not usable:
        raise ValueError(
            f"could not find {n} pad base(s) free of stop codons and {forbidden}; "
            "pick a different enzyme pair"
        )
    return usable


def _vector_supplies_stop(vector: Vector, down_end: int, frame_offset: int) -> tuple[bool, str]:
    """Walk the vector downstream of the insert, in frame, looking for a stop."""
    seq = vector.seq
    i = down_end + frame_offset
    aa = []
    limit = min(len(seq), i + DOWNSTREAM_SCAN_NT)
    while i + 3 <= limit:
        codon = seq[i : i + 3]
        if codon in STOP_CODONS:
            return True, "".join(aa)
        aa.append(str(Seq(codon).translate()))
        i += 3
    return False, "".join(aa)


class CTerm(str, Enum):
    """What happens at the C-terminus of the insert.

    The vector's own reading frame after the 3' site is not always the frame of its
    C-terminal tag (pET-28a(+): XhoI and the His6 tag sit one base off the start codon's
    frame), so the choice changes which bases are added in front of the 3' site.
    """

    VECTOR_FRAME = "vector_frame"  # keep the vector's frame after the site, as it was
    STOP = "stop"  # end the protein with a stop codon; nothing from the vector is read
    TAG = "tag"  # read into the vector's C-terminal His tag, whatever frame it is in


@dataclass(frozen=True)
class CTermTag:
    position: int  # 0-based start of the first tag codon in the vector
    residues: str  # the tag as translated, e.g. "HHHHHH"


def find_cterm_tag(vector: Vector, after: int, window: int = 240) -> CTermTag | None:
    """A His6 run in one of the three frames after ``after``, before that frame's stop."""
    best = None
    for frame in range(3):
        seg = vector.seq[after + frame : after + frame + window]
        seg = seg[: len(seg) - len(seg) % 3]
        aa = str(Seq(seg).translate(to_stop=True))
        k = aa.find("HHHHHH")
        if k >= 0 and (best is None or k < best[0]):
            run = len(aa[k:]) - len(aa[k:].lstrip("H"))
            best = (k, CTermTag(after + frame + 3 * k, "H" * run))
    return best[1] if best else None


# ----------------------------------------------------------------- main entry


def design_insert(
    plan: ConstructPlan,
    vector: Vector,
    site: InsertionSite,
    host: HostProfile,
    *,
    strategy: CodonStrategy = CodonStrategy.WEIGHTED_SAMPLE,
    seed: int | None = None,
    cterm: CTerm = CTerm.VECTOR_FRAME,
) -> InsertDesign:
    if len(plan.chains) != 1:
        raise NotImplementedError("v0.1 handles single-chain constructs only")
    chain = plan.chains[0]

    up = (
        vector.hit_at(site.upstream_enzyme, site.upstream_pos)
        if site.upstream_pos is not None
        else vector.find_enzyme(site.upstream_enzyme)
    )
    down = (
        vector.hit_at(site.downstream_enzyme, site.downstream_pos)
        if site.downstream_pos is not None
        else vector.find_enzyme(site.downstream_enzyme)
    )
    issues: list[Issue] = []

    if down.start <= up.end:
        raise ValueError(
            f"{site.downstream_enzyme} (pos {down.start}) is not downstream of "
            f"{site.upstream_enzyme} (pos {up.start}) -- v0.1 requires a non-wrapping MCS"
        )

    up_site = enzyme_site(site.upstream_enzyme)
    down_site = enzyme_site(site.downstream_enzyme)
    forbidden = [up_site, down_site]

    # --- where does translation start, in final-plasmid coordinates? ---------
    if site.start_source is StartSource.VECTOR:
        if site.vector_orf_start is None:
            raise ValueError("start_source=VECTOR requires vector_orf_start")
        orf_start = site.vector_orf_start
        if orf_start >= up.start:
            raise ValueError("vector_orf_start must lie upstream of the upstream enzyme site")
        if vector.seq[orf_start : orf_start + 3] != "ATG":
            issues.append(
                Issue(Severity.WARNING, "no_atg", f"vector_orf_start {orf_start} is not an ATG")
            )
    else:
        orf_start = up.start + site.start_offset_in_site
        codon = vector.seq[orf_start : orf_start + 3]
        if codon != "ATG":
            raise ValueError(
                f"{site.upstream_enzyme} site offset {site.start_offset_in_site} "
                f"gives {codon!r}, not ATG"
            )

    # --- frame arithmetic ---------------------------------------------------
    offset = (up.start + len(up_site) - orf_start) % 3
    n_pad_up = (3 - offset) % 3
    # The vector's own downstream frame must survive the splice, or its C-terminal tag
    # is read out of frame.  That frame is set by the vector's real start codon, which
    # may differ from ``orf_start`` (an NdeI-style site restarts translation inside the
    # tag).  This fixes the 3' pad for any enzyme pair, not only ones whose downstream
    # site happens to start in frame.
    frame_grid = vector.expression_start if vector.expression_start is not None else orf_start
    if cterm is CTerm.TAG:
        tag = find_cterm_tag(vector, down.end)
        if tag is None:
            raise ValueError(
                f"No C-terminal His tag follows the {site.downstream_enzyme} site; "
                "choose a stop codon or the vector's own frame instead"
            )
        frame_grid = tag.position
    # with a stop codon straight after the insert nothing downstream is read, so no pad
    n_pad_down = 0 if cterm is CTerm.STOP else (down.start - frame_grid) % 3

    insert_aa = chain.insert_aa

    # --- stop codon policy --------------------------------------------------
    down_frame_offset = (n_pad_down + len(down_site)) % 3
    has_downstream_stop, downstream_aa = _vector_supplies_stop(
        vector, down.end, (frame_grid - down.end) % 3
    )
    vector_encodes_cterm = bool(downstream_aa)

    if chain.stop is StopPolicy.FORCE or cterm is CTerm.STOP:
        want_stop = True
    elif chain.stop is StopPolicy.OMIT:
        want_stop = False
    else:
        want_stop = not has_downstream_stop

    if want_stop and vector_encodes_cterm:
        issues.append(
            Issue(
                # asking for a stop outright is a choice, not something to warn about
                Severity.INFO if cterm is CTerm.STOP else Severity.WARNING,
                "cterm_lost",
                f"insert carries a stop but the vector encodes {len(downstream_aa)} "
                f"residue(s) downstream ({downstream_aa[:20]}...) -- they will not be translated",
            )
        )
    if not want_stop and not has_downstream_stop:
        issues.append(
            Issue(
                Severity.ERROR,
                "no_stop",
                "no stop codon in the insert and none found in frame within "
                f"{DOWNSTREAM_SCAN_NT} nt of the vector -- the ORF would run on",
            )
        )

    stop_dna = host.best_codon("*") if want_stop else ""
    if want_stop:
        # double stop is standard practice; the second one is belt-and-braces
        stop_dna = "TAA" + "TGA"

    # --- coding region: optimised against the real junction context ---------
    # The pads are fixed DNA from the optimiser's point of view, and a site can
    # run across pad and coding region, so try pad combinations until the DP
    # finds a path.  If none does, the first combination's error is reported.
    # include up to two vector bases on each side so a stop codon straddling the site
    # boundary is seen, but never bases before the start codon (they are not translated)
    lead = max(0, min(2, up.start - orf_start))
    up_pads = _pad_candidates(
        n_pad_up,
        prefix=vector.seq[up.start - lead : up.start] + up_site,
        suffix="",
        forbidden=forbidden,
        frame_offset=(up.start - lead - orf_start) % 3,
    )
    down_pads = _pad_candidates(
        n_pad_down,
        prefix=stop_dna,
        suffix=down_site + vector.seq[down.end : down.end + 2],
        forbidden=forbidden,
        frame_offset=0,
    )
    context = Automaton(forbidden).max_len
    vector_before = vector.seq[max(0, up.start - context) : up.start]
    vector_after = vector.seq[down.end : down.end + context]

    first = None
    for pad_up, pad_down in ((u, d) for u in up_pads for d in down_pads):
        result = optimize_coding(
            OptimizationRequest(
                aa=insert_aa,
                host=host,
                prefix=vector_before + up_site + pad_up,
                suffix=stop_dna + pad_down + down_site + vector_after,
                forbidden=forbidden,
                strategy=strategy,
                seed=seed,
            )
        )
        first = first or (result, pad_up, pad_down)
        if result.ok:
            break
    else:
        result, pad_up, pad_down = first
    coding = result.dna
    issues += result.issues

    insert_dna = up_site + pad_up + coding + stop_dna + pad_down + down_site
    coding_start = up.start + len(up_site) + len(pad_up)

    # --- duplicate initiator Met --------------------------------------------
    if site.start_source is StartSource.UPSTREAM_SITE and insert_aa.startswith("M"):
        issues.append(
            Issue(
                Severity.WARNING,
                "duplicate_met",
                f"{site.upstream_enzyme} donates the initiator Met and the insert also "
                "starts with M -- the product begins MM; drop the leading M from the "
                "payload unless that is intended",
            )
        )

    # --- duplicate-tag check -------------------------------------------------
    for part in chain.insert_parts():
        if part.kind is PartKind.TAG and part.aa and part.aa in downstream_aa:
            issues.append(
                Issue(
                    Severity.WARNING,
                    "duplicate_tag",
                    f"{part.name} is encoded by the insert but the vector already "
                    "supplies the same sequence downstream",
                )
            )

    design = InsertDesign(
        vector=vector,
        site=site,
        chain=chain,
        up=up,
        down=down,
        orf_start=orf_start,
        pad_up=pad_up,
        pad_down=pad_down,
        coding=coding,
        insert_dna=insert_dna,
        coding_start=coding_start,
        issues=issues,
        stop_in_insert=want_stop,
        metrics=result.metrics,
    )
    _ = down_frame_offset
    return design


# ------------------------------------------------------------ fusion preview


def fusion_segments(design: InsertDesign, final_seq: str) -> list[Segment]:
    """Walk the final ORF codon by codon and label each by provenance.

    Codon-level provenance is used rather than slicing by part boundaries,
    because a junction codon can straddle vector and insert sequence -- which
    is exactly the case users get wrong by hand.
    """
    d = design
    insert_start = d.up.start
    insert_end = insert_start + len(d.insert_dna)
    coding_start = d.coding_start
    coding_end = coding_start + len(d.coding)

    boundaries = [
        (insert_start, "vector (N-terminal)", PartKind.TAG, PartSource.VECTOR),
        (coding_start, "scar: 5' junction", PartKind.SCAR, PartSource.INSERT),
        (coding_end, "insert", PartKind.PAYLOAD, PartSource.INSERT),
        (insert_end, "scar: 3' junction", PartKind.SCAR, PartSource.INSERT),
        (len(final_seq), "vector (C-terminal)", PartKind.TAG, PartSource.VECTOR),
    ]

    def provenance(pos: int):
        for limit, label, kind, source in boundaries:
            if pos < limit:
                return label, kind, source
        return "vector (C-terminal)", PartKind.TAG, PartSource.VECTOR

    segments: list[Segment] = []
    i = d.orf_start
    while i + 3 <= len(final_seq):
        codon = final_seq[i : i + 3]
        if codon in STOP_CODONS:
            break
        aa = str(Seq(codon).translate())
        label, kind, source = provenance(i)
        if segments and segments[-1].label == label:
            segments[-1].aa += aa
        else:
            segments.append(Segment(label, aa, kind, source))
        i += 3
    return segments


def format_fusion(segments: list[Segment]) -> str:
    lines = []
    pos = 1
    for seg in segments:
        end = pos + len(seg) - 1
        lines.append(f"  {pos:>5}-{end:<5} {seg.label:<24} {seg.aa}")
        pos = end + 1
    lines.append(f"  total: {pos - 1} aa")
    return "\n".join(lines)
