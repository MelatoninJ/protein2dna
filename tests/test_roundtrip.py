"""The hard gate.

Every case here asserts that translating the *assembled plasmid* from the
vector's own start codon reproduces the intended fusion protein.  When this
file goes red, nothing else in the package matters.
"""

from __future__ import annotations

import pytest
from Bio.Seq import Seq

import p2d
from p2d.vector import InsertionSite, StartSource

PAYLOAD = "KVFLDWINEAYQRGTRVLAEMAKRGDEFVKRLIAEGHDPFEVLKELGYSE"
NANOBODY_LIKE = "QVQLVESGGGLVQAGGSLRLSCAASGRTFSSYAMGWFRQAPGKEREFVA"


@pytest.fixture(scope="module")
def vector():
    return p2d.load_bundled()


@pytest.fixture(scope="module")
def host():
    return p2d.load_host("ecoli_bl21")


def build(vector, host, site_key, payload, **kw):
    plan = p2d.ConstructPlan.single_chain(payload, **kw)
    return p2d.assemble(p2d.design_insert(plan, vector, vector.sites[site_key], host))


# ------------------------------------------------------------------ golden set

GOLDEN = [
    ("NdeI-XhoI", PAYLOAD, "M" + PAYLOAD + "LE" + "HHHHHH"),
    ("NdeI-EcoRI", PAYLOAD, "M" + PAYLOAD + "EF" + "LEHHHHHH"),
    (
        "BamHI-XhoI",
        PAYLOAD,
        "MGSSHHHHHHSSGLVPRGSHM" + "GS" + PAYLOAD + "LE" + "HHHHHH",
    ),
    ("NdeI-XhoI", NANOBODY_LIKE, "M" + NANOBODY_LIKE + "LE" + "HHHHHH"),
]


@pytest.mark.parametrize("site_key,payload,expected", GOLDEN)
def test_golden_fusion_protein(vector, host, site_key, payload, expected):
    a = build(vector, host, site_key, payload)
    assert a.ok, a.errors
    assert a.protein == expected


@pytest.mark.parametrize("site_key,payload,expected", GOLDEN)
def test_translation_matches_segments(vector, host, site_key, payload, expected):
    """The fusion preview must agree with the actual DNA."""
    a = build(vector, host, site_key, payload)
    assert "".join(s.aa for s in a.segments) == a.protein


@pytest.mark.parametrize("site_key,payload,expected", GOLDEN)
def test_orf_terminates(vector, host, site_key, payload, expected):
    a = build(vector, host, site_key, payload)
    tail = a.seq[a.design.orf_start + 3 * len(a.protein) :][:3]
    assert tail in {"TAA", "TAG", "TGA"}


# -------------------------------------------------------------- frame rescue


@pytest.mark.parametrize("shift", [0, 1, 2])
def test_pad_rescues_out_of_frame_site(vector, host, shift):
    """A site that is out of frame with the ORF must be rescued by padding."""
    site = InsertionSite(
        "shifted", "BamHI", "XhoI", StartSource.VECTOR, vector_orf_start=87 - shift
    )
    plan = p2d.ConstructPlan.single_chain(PAYLOAD)
    d = p2d.design_insert(plan, vector, site, host)
    a = p2d.assemble(d)
    assert (d.coding_start - d.orf_start) % 3 == 0
    assert len(d.pad_up) == (3 - shift) % 3 if shift else len(d.pad_up) == 0
    assert a.ok, a.errors
    assert PAYLOAD in a.protein
    assert a.protein.endswith("LEHHHHHH")


def test_insert_is_multiple_of_three_from_start(vector, host):
    a = build(vector, host, "NdeI-XhoI", PAYLOAD)
    d = a.design
    assert (d.coding_start - d.orf_start) % 3 == 0


# ------------------------------------------------------------- coding fidelity


def test_coding_backtranslates(vector, host):
    a = build(vector, host, "NdeI-XhoI", PAYLOAD)
    assert str(Seq(a.design.coding).translate()) == PAYLOAD


def test_no_internal_cloning_sites(vector, host):
    """Synonymous repair must clear the cloning enzymes from the coding region."""
    # a payload whose naive reverse translation contains CTCGAG (XhoI: Leu-Glu)
    payload = "AALEAA" + PAYLOAD
    a = build(vector, host, "NdeI-XhoI", payload)
    assert a.ok, a.errors
    assert "CTCGAG" not in a.design.coding
    assert "CATATG" not in a.design.coding
    assert (
        str(Seq(a.design.coding).translate()) == "M" + payload
        or str(Seq(a.design.coding).translate()) == payload
    )


def test_enzyme_sites_unique_in_final_plasmid(vector, host):
    from Bio.Restriction import RestrictionBatch

    a = build(vector, host, "NdeI-XhoI", PAYLOAD)
    for enzyme in ("NdeI", "XhoI"):
        n = len(list(RestrictionBatch([enzyme]))[0].search(Seq(a.seq), linear=False))
        assert n == 1, f"{enzyme} appears {n} times in the assembled plasmid"


# ----------------------------------------------------------------- behaviour


def test_duplicate_met_warned(vector, host):
    a = build(vector, host, "NdeI-XhoI", "M" + PAYLOAD)
    assert any(i.code == "duplicate_met" for i in a.issues)


def test_duplicate_tag_warned(vector, host):
    from p2d.model import Chain, Part, PartKind, PartSource

    chain = Chain(
        name="tagged",
        payload=PAYLOAD,
        c_terminal=[Part("His6", "HHHHHH", PartKind.TAG, PartSource.INSERT)],
    )
    plan = p2d.ConstructPlan(chains=[chain])
    a = p2d.assemble(p2d.design_insert(plan, vector, vector.sites["NdeI-XhoI"], host))
    assert any(i.code == "duplicate_tag" for i in a.issues)


def test_stop_policy_auto_omits_when_vector_has_cterm(vector, host):
    a = build(vector, host, "NdeI-XhoI", PAYLOAD)
    assert not a.design.stop_in_insert
    assert a.protein.endswith("HHHHHH")


def test_forced_stop_truncates_and_warns(vector, host):
    from p2d.model import Chain, StopPolicy

    plan = p2d.ConstructPlan(chains=[Chain("x", PAYLOAD, stop=StopPolicy.FORCE)])
    a = p2d.assemble(p2d.design_insert(plan, vector, vector.sites["NdeI-XhoI"], host))
    assert any(i.code == "cterm_lost" for i in a.issues)
    assert not a.protein.endswith("HHHHHH")


def test_multichain_rejected_in_v01():
    from p2d.model import Chain

    with pytest.raises(ValueError):
        p2d.ConstructPlan(chains=[Chain("a", PAYLOAD), Chain("b", PAYLOAD)])


def test_unimplemented_linkage_is_explicit():
    from p2d.model import Chain, LinkageStrategy

    with pytest.raises(NotImplementedError):
        p2d.ConstructPlan(chains=[Chain("a", PAYLOAD)], linkage=LinkageStrategy.SELF_CLEAVING_2A)
