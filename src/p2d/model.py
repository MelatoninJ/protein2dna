"""Construct data model.

The central idea: a construct is a *list of chains*, even when there is only
one.  Single-chain is just ``len(chains) == 1``.  This keeps the multi-chain
extension (Duet vectors, 2A peptides, co-transfection) from requiring a
rewrite of every call site later.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from enum import Enum

AA_ALPHABET = set("ACDEFGHIKLMNPQRSTVWY")


class PartKind(str, Enum):
    """What role a stretch of protein plays in the construct."""

    SIGNAL_PEPTIDE = "signal_peptide"
    TAG = "tag"
    LINKER = "linker"
    PROTEASE_SITE = "protease_site"
    PAYLOAD = "payload"
    SCAR = "scar"  # residues forced by a restriction site or a frame pad


class PartSource(str, Enum):
    """Who encodes this part.

    VECTOR parts already exist in the plasmid; the insert must *not* encode
    them again.  This single field is what catches the most common construct
    bug there is: adding a His-tag to a protein going into pET-28a, which
    already supplies one.
    """

    VECTOR = "vector"
    INSERT = "insert"


@dataclass(frozen=True)
class Part:
    name: str
    aa: str
    kind: PartKind
    source: PartSource = PartSource.INSERT

    def __post_init__(self) -> None:
        bad = set(self.aa.upper()) - AA_ALPHABET
        if bad:
            raise ValueError(f"part {self.name!r} has non-standard residues: {sorted(bad)}")
        object.__setattr__(self, "aa", self.aa.upper())

    def __len__(self) -> int:
        return len(self.aa)


class StopPolicy(str, Enum):
    AUTO = "auto"  # stop iff nothing in-frame follows in the vector
    FORCE = "force"
    OMIT = "omit"


@dataclass
class Chain:
    """One polypeptide the insert must encode."""

    name: str
    payload: str
    n_terminal: list[Part] = field(default_factory=list)
    c_terminal: list[Part] = field(default_factory=list)
    stop: StopPolicy = StopPolicy.AUTO

    def __post_init__(self) -> None:
        self.payload = self.payload.upper().strip().replace("*", "")
        bad = set(self.payload) - AA_ALPHABET
        if bad:
            raise ValueError(
                f"chain {self.name!r} payload has non-standard residues: {sorted(bad)}"
            )

    @property
    def payload_part(self) -> Part:
        return Part(self.name, self.payload, PartKind.PAYLOAD, PartSource.INSERT)

    def insert_parts(self) -> list[Part]:
        """Parts the insert itself must encode, N-to-C, excluding scars."""
        parts = [p for p in self.n_terminal if p.source is PartSource.INSERT]
        parts.append(self.payload_part)
        parts += [p for p in self.c_terminal if p.source is PartSource.INSERT]
        return parts

    @property
    def insert_aa(self) -> str:
        return "".join(p.aa for p in self.insert_parts())


class LinkageStrategy(str, Enum):
    """How multiple chains are connected.

    Only SINGLE is implemented in v0.1.  The others are registered so the API
    shape is settled now; selecting one raises NotImplementedError with a
    pointer to the roadmap rather than silently doing the wrong thing.
    """

    SINGLE = "single"
    DUAL_MCS = "dual_mcs"  # pETDuet / pACYCDuet style
    POLYCISTRONIC = "polycistronic"  # one promoter, one RBS per ORF
    SELF_CLEAVING_2A = "self_cleaving_2a"  # P2A/T2A polyprotein
    SEPARATE_PLASMIDS = "separate_plasmids"


@dataclass
class ConstructPlan:
    chains: list[Chain]
    host: str = "ecoli_bl21"
    linkage: LinkageStrategy = LinkageStrategy.SINGLE

    def __post_init__(self) -> None:
        if not self.chains:
            raise ValueError("a construct needs at least one chain")
        if self.linkage is LinkageStrategy.SINGLE and len(self.chains) > 1:
            raise ValueError(
                f"{len(self.chains)} chains given but linkage is 'single'; "
                "pick a multi-chain strategy (not implemented in v0.1)"
            )
        if self.linkage is not LinkageStrategy.SINGLE:
            raise NotImplementedError(
                f"linkage strategy {self.linkage.value!r} is planned for v1.0. "
                "v0.1 implements 'single' only."
            )

    @classmethod
    def single_chain(cls, payload: str, *, name: str = "payload", **kw) -> ConstructPlan:
        return cls(chains=[Chain(name=name, payload=payload)], **kw)


# ---------------------------------------------------------------- tag library

TAGS: dict[str, Part] = {
    "his6": Part("His6", "HHHHHH", PartKind.TAG),
    "his8": Part("His8", "HHHHHHHH", PartKind.TAG),
    "flag": Part("FLAG", "DYKDDDDK", PartKind.TAG),
    "myc": Part("Myc", "EQKLISEEDL", PartKind.TAG),
    "strep2": Part("Strep-II", "WSHPQFEK", PartKind.TAG),
    "ha": Part("HA", "YPYDVPDYA", PartKind.TAG),
}

LINKERS: dict[str, Part] = {
    "gs": Part("GS", "GS", PartKind.LINKER),
    "g4s": Part("G4S", "GGGGS", PartKind.LINKER),
    "g4s3": Part("(G4S)3", "GGGGSGGGGSGGGGS", PartKind.LINKER),
    "eaaak": Part("EAAAK", "EAAAK", PartKind.LINKER),
}

PROTEASE_SITES: dict[str, Part] = {
    # cleavage is C-terminal to the last residue listed unless noted
    "tev": Part("TEV", "ENLYFQG", PartKind.PROTEASE_SITE),
    "3c": Part("HRV 3C", "LEVLFQGP", PartKind.PROTEASE_SITE),
    "thrombin": Part("thrombin", "LVPRGS", PartKind.PROTEASE_SITE),
    "enterokinase": Part("enterokinase", "DDDDK", PartKind.PROTEASE_SITE),
}

SIGNAL_PEPTIDES: dict[str, Part] = {
    # periplasmic export in E. coli
    "pelb": Part("pelB", "MKYLLPTAAAGLLLLAAQPAMA", PartKind.SIGNAL_PEPTIDE),
    "ompa": Part("OmpA", "MKKTAIAIAVALAGFATVAQA", PartKind.SIGNAL_PEPTIDE),
    "dsba": Part("DsbA", "MKKIWLALAGLVLAFSASA", PartKind.SIGNAL_PEPTIDE),
}
