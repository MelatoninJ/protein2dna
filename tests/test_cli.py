"""The command line is a thin shell; these tests keep it honest about output and exit codes."""

from Bio import SeqIO
from typer.testing import CliRunner

import p2d
from p2d.cli import app

PROTEIN = "KVFLDWINEAYQRGTRVLAEMAKRGDEFVKRLIAEGHDPFEVLKELGYSE"
runner = CliRunner()


def test_version_matches_the_package():
    r = runner.invoke(app, ["version"])
    assert r.exit_code == 0 and r.output.strip() == f"p2d {p2d.__version__}"


def test_hosts_and_vectors_are_listed():
    assert "ecoli_bl21" in runner.invoke(app, ["hosts"]).output
    out = runner.invoke(app, ["vectors"]).output
    assert "pTEST1" in out and "NdeI-XhoI" in out


def test_design_reports_the_fusion_protein(tmp_path):
    gb, fa = tmp_path / "out.gb", tmp_path / "insert.fa"
    r = runner.invoke(app, ["design", PROTEIN, "--gb", str(gb), "--fasta", str(fa)])
    assert r.exit_code == 0, r.output
    assert "pTEST1 / NdeI-XhoI" in r.output and "total: 79 aa" in r.output
    record = next(SeqIO.parse(str(gb), "genbank"))
    insert = next(SeqIO.parse(str(fa), "fasta"))
    assert str(insert.seq) in str(record.seq)
    assert str(insert.seq) == p2d.run(PROTEIN, site="NdeI-XhoI").insert_dna


def test_a_fasta_file_can_be_given_instead_of_a_sequence(tmp_path):
    path = tmp_path / "p.fa"
    path.write_text(f">binder\n{PROTEIN}\n")
    from_file = runner.invoke(app, ["design", str(path)])
    inline = runner.invoke(app, ["design", PROTEIN])
    assert from_file.exit_code == 0 and from_file.output == inline.output


def test_another_site_of_the_bundled_vector_can_be_chosen():
    sites = list(p2d.load_bundled().sites)
    other = next(s for s in sites if s != "NdeI-XhoI")
    r = runner.invoke(app, ["design", PROTEIN, "--site", other])
    assert f"pTEST1 / {other}" in r.output


def test_a_vector_without_the_named_site_is_refused(tmp_path):
    vec = tmp_path / "v.gb"
    SeqIO.write(p2d.load_bundled().record, str(vec), "genbank")
    r = runner.invoke(app, ["design", PROTEIN, "--vector", str(vec), "--site", "Nope-Nope"])
    assert r.exit_code == 2 and "no registered insertion site" in r.output


def test_a_residue_the_codon_table_lacks_fails_cleanly():
    r = runner.invoke(app, ["design", "MKV#Z"])
    assert r.exit_code != 0


def test_ui_starts_the_local_server(monkeypatch):
    import p2d_ui

    calls = []
    monkeypatch.setattr(p2d_ui, "serve", lambda **kw: calls.append(kw))
    r = runner.invoke(app, ["ui", "--port", "9123", "--no-browser"])
    assert r.exit_code == 0 and calls == [{"port": 9123, "open_browser": False}]
