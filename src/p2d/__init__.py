"""p2d -- protein-to-DNA reverse translation that knows about your vector.

Status: v0.1. The frame engine and round-trip validation are implemented; the
constrained codon optimiser is not (see README roadmap).  Reverse translation
currently picks the most-used codon per residue, which is explicitly *not* a
good production strategy -- treat v0.1 output as a frame check, not as a
sequence to order.
"""

from .assemble import Assembly, assemble
from .frame import InsertDesign, Issue, Severity, design_insert
from .host import HostProfile, available_hosts, load_host
from .model import (
    LINKERS,
    PROTEASE_SITES,
    SIGNAL_PEPTIDES,
    TAGS,
    Chain,
    ConstructPlan,
    LinkageStrategy,
    Part,
    PartKind,
    PartSource,
    StopPolicy,
)
from .vector import InsertionSite, StartSource, Vector
from .vectors import available_bundled, load_bundled

__version__ = "0.1.0"

__all__ = [
    "Assembly",
    "Chain",
    "ConstructPlan",
    "HostProfile",
    "InsertDesign",
    "InsertionSite",
    "Issue",
    "LINKERS",
    "LinkageStrategy",
    "PROTEASE_SITES",
    "Part",
    "PartKind",
    "PartSource",
    "SIGNAL_PEPTIDES",
    "Severity",
    "StartSource",
    "StopPolicy",
    "TAGS",
    "Vector",
    "__version__",
    "assemble",
    "available_bundled",
    "available_hosts",
    "design_insert",
    "load_bundled",
    "load_host",
    "run",
]


def run(payload: str, *, vector=None, site: str = "NdeI-XhoI", host: str = "ecoli_bl21", **kw):
    """One-call convenience path: payload in, validated assembly out."""
    v = vector if vector is not None else load_bundled()
    plan = ConstructPlan.single_chain(payload, host=host, **kw)
    design = design_insert(plan, v, v.sites[site], load_host(host))
    return assemble(design)
