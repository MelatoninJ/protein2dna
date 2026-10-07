"""Bundled vectors.

v0.1 ships exactly one: ``pTEST1``, a small synthetic pET-like vector used by
the test suite.  It is **not** a real plasmid and must not be ordered or used
at the bench.  Real vectors come from :func:`p2d.vector.Vector.from_genbank`
(user upload) or, from v0.5, the Addgene developers API.

Adding a vector here means supplying its ``InsertionSite`` definitions -- the
enzyme pair plus where the reading frame starts.  That is the only part a
GenBank file cannot tell you reliably.
"""

from __future__ import annotations

from importlib import resources

from .vector import InsertionSite, StartSource, Vector

DATA_PACKAGE = "p2d.data"

PTEST1_SITES = {
    # NdeI donates the ATG: the vector's own N-terminal tag is excised.
    "NdeI-XhoI": InsertionSite(
        name="NdeI-XhoI",
        upstream_enzyme="NdeI",
        downstream_enzyme="XhoI",
        start_source=StartSource.UPSTREAM_SITE,
        start_offset_in_site=3,  # CAT|ATG
        note="untagged N-terminus; C-terminal His6 retained",
    ),
    "NdeI-EcoRI": InsertionSite(
        name="NdeI-EcoRI",
        upstream_enzyme="NdeI",
        downstream_enzyme="EcoRI",
        start_source=StartSource.UPSTREAM_SITE,
        start_offset_in_site=3,
        note="untagged N-terminus; XhoI scar and C-terminal His6 retained",
    ),
    # BamHI lies downstream of the vector ATG, so the N-terminal tag is kept.
    "BamHI-XhoI": InsertionSite(
        name="BamHI-XhoI",
        upstream_enzyme="BamHI",
        downstream_enzyme="XhoI",
        start_source=StartSource.VECTOR,
        vector_orf_start=87,
        note="N-terminal His6 + thrombin retained; adds HM scar from the NdeI site",
    ),
}


def load_bundled(name: str = "pTEST1") -> Vector:
    if name != "pTEST1":
        raise ValueError(f"no bundled vector named {name!r}; available: pTEST1")
    path = resources.files(DATA_PACKAGE).joinpath("pTEST1.gb")
    with resources.as_file(path) as p:
        v = Vector.from_genbank(p, name="pTEST1")
    v.sites = dict(PTEST1_SITES)
    return v


def available_bundled() -> list[str]:
    return ["pTEST1"]
