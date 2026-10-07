"""Tests for the vector layer, host profiles, and the architecture boundary."""

from __future__ import annotations

import ast
from pathlib import Path

import pytest

import p2d

SRC = Path(p2d.__file__).parent
UI_MODULES = {"streamlit", "gradio", "flask", "fastapi", "dash", "django"}


def test_library_imports_no_ui_framework():
    """The core must stay usable from a batch script.

    This is the rule that keeps the eventual web UI a thin shell instead of
    the place the logic ends up living.
    """
    offenders = []
    for path in SRC.rglob("*.py"):
        tree = ast.parse(path.read_text())
        for node in ast.walk(tree):
            if isinstance(node, ast.Import):
                names = [a.name.split(".")[0] for a in node.names]
            elif isinstance(node, ast.ImportFrom):
                names = [(node.module or "").split(".")[0]]
            else:
                continue
            for n in names:
                if n in UI_MODULES:
                    offenders.append(f"{path.relative_to(SRC)}: {n}")
    assert not offenders, f"UI framework imported in the library: {offenders}"


# ------------------------------------------------------------------ vectors


def test_bundled_vector_loads():
    v = p2d.load_bundled()
    assert len(v) == 1149
    assert v.circular


@pytest.mark.parametrize("enzyme", ["NdeI", "BamHI", "EcoRI", "XhoI"])
def test_cloning_enzymes_are_unique(enzyme):
    assert p2d.load_bundled().cuts_once(enzyme)


def test_find_enzyme_rejects_non_unique():
    v = p2d.Vector.from_sequence("GGATCC" * 2 + "A" * 50)
    with pytest.raises(ValueError, match="not unique"):
        v.find_enzyme("BamHI")


def test_find_enzyme_rejects_absent():
    v = p2d.Vector.from_sequence("A" * 60)
    with pytest.raises(ValueError, match="does not cut"):
        v.find_enzyme("BamHI")


def test_orf_detected_from_sequence_not_labels():
    """ORF detection must not depend on the record's CDS annotation."""
    v = p2d.load_bundled()
    orfs = v.detect_orfs(min_aa=30)
    assert (87, 87 + 37 * 3 + 3) in orfs


def test_empty_vector_orf_translates_as_expected():
    from Bio.Seq import Seq

    v = p2d.load_bundled()
    assert str(Seq(v.seq[87:]).translate(to_stop=True)) == "MGSSHHHHHHSSGLVPRGSHMGSGTGSEFLEHHHHHH"


# -------------------------------------------------------------------- hosts


def test_host_loads_and_flags_placeholder_table():
    with pytest.warns(UserWarning, match="placeholder"):
        h = p2d.load_host("ecoli_bl21")
    assert not h.codon_table_verified
    assert h.best_codon("M") == "ATG"
    assert h.best_codon("*") == "TAA"


def test_every_residue_has_a_codon():
    h = p2d.load_host("ecoli_bl21")
    for aa in "ACDEFGHIKLMNPQRSTVWY":
        assert h.best_codon(aa)


def test_codon_table_is_self_consistent():
    """Each amino acid's codons must actually encode that amino acid."""
    from Bio.Seq import Seq

    h = p2d.load_host("ecoli_bl21")
    for aa, codons in h.codon_usage.items():
        for codon in codons:
            assert str(Seq(codon).translate()) == aa, f"{codon} is not {aa}"


def test_codon_table_is_complete():
    """All 61 sense codons must be present exactly once."""
    seen = [c for codons in p2d.load_host("ecoli_bl21").codon_usage.values() for c in codons]
    assert len(seen) == len(set(seen)) == 61


def test_unknown_host_raises():
    with pytest.raises(ValueError, match="unknown host"):
        p2d.load_host("nonexistent_organism")


# -------------------------------------------------------------------- output


def test_genbank_record_roundtrips(tmp_path):
    from Bio import SeqIO

    a = p2d.run("KVFLDWINEAYQRGTRVLAEM")
    path = tmp_path / "out.gb"
    SeqIO.write(a.to_record(), str(path), "genbank")
    back = SeqIO.read(str(path), "genbank")
    assert str(back.seq) == a.seq
    assert any("insert" in f.qualifiers.get("label", [""])[0] for f in back.features)
