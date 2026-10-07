"""p2d -- protein-to-DNA reverse translation that knows about your vector.

Status: v0.2. The frame engine, round-trip validation and the DFA-constrained
codon optimiser are implemented.  GC windows, mRNA folding and synthesis checks
are not (see README roadmap) -- review a design before ordering it.
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
from .optimize import CodonStrategy
from .vector import InsertionSite, StartSource, Vector
from .vectors import available_bundled, load_bundled

__version__ = "0.1.0"

__all__ = [
    "Assembly",
    "Chain",
    "CodonStrategy",
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
