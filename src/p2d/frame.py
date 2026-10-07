"""The frame engine.

Everything here exists to answer one question precisely: *given this vector,
these two enzymes and this payload, what protein actually comes out?*

The arithmetic that matters:

    final = vector[:up.start] + insert + vector[down.end:]
    insert = up_site + pad_up + coding + pad_down + down_site

``pad_up`` makes the payload's first codon land in frame with the start codon.
``pad_down`` makes whatever the vector encodes downstream land in frame.  Both
pads are chosen by a small search that refuses bases creating a stop codon or a
new copy of either enzyme's site.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from enum import Enum

from Bio.Seq import Seq

from .host import HostProfile
from .model import Chain, ConstructPlan, PartKind, PartSource, StopPolicy
from .vector import EnzymeHit, InsertionSite, StartSource, Vector, enzyme_site, site_regex

PAD_BASES = ["GC", "GG", "CG", "GA", "AC", "CA", "TC", "CT", "AG", "GT"]
STOP_CODONS = {"TAA", "TAG", "TGA"}
DOWNSTREAM_SCAN_NT = 600


class Severity(str, Enum):
    ERROR = "error"
    WARNING = "warning"
    INFO = "info"


@dataclass(frozen=True)
class Issue:
    severity: Severity
    code: str
    message: str

    def __str__(self) -> str:
        return f"[{self.severity.value}] {self.code}: {self.message}"


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

    @property
    def errors(self) -> list[Issue]:
        return [i for i in self.issues if i.severity is Severity.ERROR]

    @property
    def ok(self) -> bool:
        return not self.errors


# --------------------------------------------------------------------- helpers


def reverse_translate(aa: str, host: HostProfile) -> str:
    """Naive codon choice (v0.1).

    This is deliberately the simplest thing that works, so the frame engine can
    be tested on its own.  The constrained optimiser replaces it in v0.2 -- the
    signature is what the rest of the package depends on.
    """
    return "".join(host.best_codon(r) for r in aa)


def _creates_problem(seq: str, forbidden: list[str]) -> bool:
    return any(site_regex(site).search(seq) for site in forbidden)


def _choose_pad(
    n: int,
    prefix: str,
    suffix: str,
    forbidden: list[str],
    frame_offset: int,
) -> str:
    """Pick ``n`` pad bases that introduce neither a stop codon nor a new site.

    ``frame_offset`` is the position of ``prefix``'s first base within its
    codon, so stop codons can be read in the correct frame.
    """
    if n == 0:
        return ""
    candidates = [p[:n] for p in PAD_BASES] if n <= 2 else ["GCG", "GGC"]
    pad_lo, pad_hi = len(prefix), len(prefix) + n
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
        return pad
    raise ValueError(
        f"could not find {n} pad base(s) free of stop codons and {forbidden}; "
        "pick a different enzyme pair"
    )


def _repair_internal_sites(
    coding: str, aa: str, host: HostProfile, forbidden: list[str]
) -> tuple[str, list[Issue]]:
    """Remove forbidden sites from the coding region by synonymous swap.

    Greedy and local -- good enough for v0.1.  The DFA/DP optimiser in v0.2
    makes this guarantee global instead of best-effort.
    """
    issues: list[Issue] = []
    for site in forbidden:
        rx = site_regex(site)
        for _ in range(50):
            m = rx.search(coding)
            if not m:
                break
            fixed = False
            for pos in range(m.start() // 3, (m.end() + 2) // 3):
                if pos >= len(aa):
                    break
                residue = aa[pos]
                options = sorted(host.codons_for(residue).items(), key=lambda kv: -kv[1])
                current = coding[pos * 3 : pos * 3 + 3]
                for codon, _freq in options:
                    if codon == current:
                        continue
                    trial = coding[: pos * 3] + codon + coding[pos * 3 + 3 :]
                    if not rx.search(trial[max(0, m.start() - 6) : m.end() + 6]):
                        coding = trial
                        issues.append(
                            Issue(
                                Severity.INFO,
                                "silent_change",
                                f"codon {pos + 1} ({residue}) {current}->{codon} to remove {site}",
                            )
                        )
                        fixed = True
                        break
                if fixed:
                    break
            if not fixed:
                issues.append(
                    Issue(
                        Severity.ERROR,
                        "unremovable_site",
                        f"{site} at coding position {m.start()} cannot be removed synonymously",
                    )
                )
                break
    return coding, issues


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


# ----------------------------------------------------------------- main entry


def design_insert(
    plan: ConstructPlan,
    vector: Vector,
    site: InsertionSite,
    host: HostProfile,
) -> InsertDesign:
    if len(plan.chains) != 1:
        raise NotImplementedError("v0.1 handles single-chain constructs only")
    chain = plan.chains[0]

    up = vector.find_enzyme(site.upstream_enzyme)
    down = vector.find_enzyme(site.downstream_enzyme)
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
    n_pad_down = (-len(down_site)) % 3

    insert_aa = chain.insert_aa
    coding = reverse_translate(insert_aa, host)
    coding, repair_issues = _repair_internal_sites(coding, insert_aa, host, forbidden)
    issues += repair_issues

    pad_up = _choose_pad(
        n_pad_up,
        prefix=up_site,
        suffix=coding[:6],
        forbidden=forbidden,
        frame_offset=(up.start - orf_start) % 3,
    )

    # --- stop codon policy --------------------------------------------------
    down_frame_offset = (n_pad_down + len(down_site)) % 3
    has_downstream_stop, downstream_aa = _vector_supplies_stop(vector, down.end, 0)
    vector_encodes_cterm = bool(downstream_aa)

    if chain.stop is StopPolicy.FORCE:
        want_stop = True
    elif chain.stop is StopPolicy.OMIT:
        want_stop = False
    else:
        want_stop = not has_downstream_stop

    if want_stop and vector_encodes_cterm:
        issues.append(
            Issue(
                Severity.WARNING,
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

    pad_down = _choose_pad(
        n_pad_down,
        prefix=coding[-6:] + stop_dna,
        suffix=down_site,
        forbidden=forbidden,
        frame_offset=0,
    )

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
