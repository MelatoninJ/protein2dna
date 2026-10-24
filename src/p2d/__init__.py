"""p2d: reverse-translate a protein into DNA that fits your vector.

It does not fold mRNA or predict expression. Check a design before you order it.
"""

from .assemble import Assembly, assemble
from .frame import CTerm, InsertDesign, Issue, Severity, design_insert
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
from .optimize import CodonStrategy, Limits
from .order import Order, make_order
from .vector import InsertionSite, StartSource, Vector
from .vectors import available_bundled, load_bundled

__version__ = "0.3.0.dev1"

__all__ = [
    "Assembly",
    "Chain",
    "CodonStrategy",
    "ConstructPlan",
    "CTerm",
    "HostProfile",
    "InsertDesign",
    "InsertionSite",
    "Issue",
    "LINKERS",
    "LinkageStrategy",
    "Limits",
    "PROTEASE_SITES",
    "Order",
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
    "make_order",
    "run",
]


def run(payload: str, *, vector=None, site: str = "NdeI-XhoI", host: str = "ecoli_bl21", **kw):
    """One-call convenience path: payload in, validated assembly out."""
    v = vector if vector is not None else load_bundled()
    plan = ConstructPlan.single_chain(payload, host=host, **kw)
    design = design_insert(plan, v, v.sites[site], load_host(host))
    return assemble(design)
