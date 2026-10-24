"""Designing many proteins at once: reading, isolating failures, and what is written."""

from __future__ import annotations

import csv
import io
import re

import pytest
from Bio import SeqIO
from hypothesis import given
from hypothesis import strategies as st
from typer.testing import CliRunner

import p2d
from p2d import batch as b
from p2d.cli import app
from p2d.enzymes import list_cut_sites, resolve_site

A = "KVFLDWINEAYQRGTRVLAEMAKRGDEFVKRLIAEGHDPFEVLKELGYSE"
B = "MKVLAAGIVGLLLAAQPAMAHHHHHH"
C = "SKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLKFICTTGKLPVPWPTL"
FASTA = f">binder_A first design\n{A[:30]}\n{A[30:]}\n>binder_B\n{B}\n>binder_C\n{C}\n"


@pytest.fixture(scope="module")
def vector():
    return p2d.load_bundled()


@pytest.fixture(scope="module")
def host():
    return p2d.load_host("ecoli_bl21")


@pytest.fixture(scope="module")
def sites(vector):
    return resolve_site(vector, "NdeI"), resolve_site(vector, "XhoI")


def test_records_are_read_with_wrapped_lines_and_headers():
    got = b.read_proteins(FASTA)
    assert [p.name for p in got] == ["binder_A", "binder_B", "binder_C"]
    assert [p.sequence for p in got] == [A, B, C] and not any(p.error for p in got)


def test_numbers_spaces_and_a_final_stop_are_ignored():
    (p,) = b.read_proteins(">x\n 1 MKVL  AAGI\n  11 VGLL*\n")
    assert p.sequence == "MKVLAAGIVGLL" and p.error is None


def test_a_sequence_without_a_header_still_counts():
    (p,) = b.read_proteins("MKVLAAGI\n")
    assert p.name == "protein_1" and p.sequence == "MKVLAAGI"


def test_duplicate_names_are_made_unique_whatever_the_case():
    names = [p.name for p in b.read_proteins(">a\nMK\n>A\nMK\n>a\nMK\n")]
    assert names == ["a", "A_2", "a_3"]
    assert len({n.lower() for n in names}) == 3


def test_a_bad_record_carries_its_error_and_the_rest_are_unaffected():
    got = b.read_proteins(">ok\nMKV\n>bad\nMK#Z\n>empty\n>fine\nMKVL\n", max_residues=100)
    assert [bool(p.error) for p in got] == [False, True, True, False]
    assert "#" in got[1].error and got[2].error == "no sequence"
    assert b.read_proteins(">long\n" + "M" * 11, max_residues=10)[0].error.startswith("11 residues")


def test_no_records_gives_an_empty_list():
    assert b.read_proteins("") == [] and b.read_proteins("\n  \n") == []


@given(st.text(max_size=120))
def test_names_are_always_safe_file_names(header):
    name = b.safe_name(header, "fallback")
    assert re.fullmatch(r"[A-Za-z0-9_.\-]+", name) and len(name) <= b.MAX_NAME
    assert "/" not in name and "\\" not in name and not name.startswith(".")
    assert name not in (".", ".."), "never a directory reference"


def test_unsafe_headers_cannot_escape_the_output_folder():
    assert b.safe_name("../../etc/passwd x", "f") == "etc_passwd"
    assert b.safe_name("   ", "f") == "f" and b.safe_name("###", "f") == "f"


def test_each_protein_matches_a_design_made_on_its_own(vector, host, sites):
    up, down = sites
    result = b.design_batch(
        b.read_proteins(FASTA), vector, up, down, host, strategy=p2d.CodonStrategy.MAX_CAI
    )
    assert result.ok and result.n_ok == 3
    for entry in result.entries:
        single = p2d.assemble(
            p2d.design_insert(
                p2d.ConstructPlan.single_chain(entry.protein, name=entry.name),
                vector,
                p2d.enzymes.plan_cut(vector, up, down),
                host,
                strategy=p2d.CodonStrategy.MAX_CAI,
                cterm=p2d.CTerm.STOP,
            )
        )
        assert entry.assembly.insert_dna == single.insert_dna
        assert entry.assembly.seq == single.seq and entry.order.ok


def test_one_failure_does_not_stop_the_others(vector, host, sites):
    got = b.read_proteins(">one\nMKV\n>two\nMK#\n>three\nMKVL\n")
    result = b.design_batch(got, vector, *sites, host)
    assert [e.ok for e in result.entries] == [True, False, True]
    assert result.entries[1].message and not result.ok and result.n_ok == 2


def test_a_design_level_failure_is_recorded_per_protein(vector, host, sites, monkeypatch):
    """A refusal from the frame engine (no tag after the site) lands on every protein."""
    monkeypatch.setattr("p2d.frame.find_cterm_tag", lambda *a, **k: None)
    result = b.design_batch(b.read_proteins(FASTA), vector, *sites, host, cterm=p2d.CTerm.TAG)
    assert result.n_ok == 0 and len(result.entries) == 3
    assert all("No C-terminal tag" in e.message for e in result.entries)


def test_tags_and_extra_bases_reach_every_entry(vector, host, sites):
    result = b.design_batch(
        b.read_proteins(FASTA),
        vector,
        *sites,
        host,
        clamp_5=8,
        clamp_3=0,
        n_parts=(p2d.TAGS["flag"],),
    )
    for e in result.entries:
        assert "DYKDDDDK" + e.protein in e.assembly.protein
        assert len(e.order.clamp_5) == 8 and e.order.clamp_3 == ""


def test_a_bad_pair_of_sites_is_the_callers_error_not_a_protein_failure(vector, host, sites):
    up, down = sites
    with pytest.raises(ValueError):
        b.design_batch(b.read_proteins(FASTA), vector, down, up, host)


@pytest.fixture(scope="module")
def made(vector, host, sites):
    got = b.read_proteins(FASTA + ">oops\nMK#\n")
    return b.design_batch(got, vector, *sites, host, strategy=p2d.CodonStrategy.MAX_CAI)


def test_the_summary_has_one_row_per_protein(made):
    rows = list(csv.DictReader(io.StringIO(b.summary_csv(made))))
    assert [r["name"] for r in rows] == ["binder_A", "binder_B", "binder_C", "oops"]
    assert [r["status"] for r in rows] == ["ok", "ok", "ok", "failed"]
    assert list(rows[0]) == b.SUMMARY_COLUMNS
    ok = rows[0]
    assert int(ok["residues"]) == len(A) and int(ok["plasmid_bp"]) > int(ok["insert_bp"])
    assert ok["order_dna"] == made.entries[0].order.dna and ok["error"] == ""
    assert "#" in rows[3]["error"] and rows[3]["order_dna"] == ""


def test_orders_and_inserts_list_only_what_was_made(made):
    orders = list(SeqIO.parse(io.StringIO(b.orders_fasta(made)), "fasta"))
    assert [r.id for r in orders] == ["binder_A_order", "binder_B_order", "binder_C_order"]
    assert [str(r.seq) for r in orders] == [e.order.dna for e in made.entries[:3]]
    inserts = list(SeqIO.parse(io.StringIO(b.inserts_fasta(made)), "fasta"))
    assert [str(r.seq) for r in inserts] == [e.assembly.insert_dna for e in made.entries[:3]]


def test_every_genbank_file_holds_its_plasmid(made):
    files = b.files(made)
    assert sorted(k for k in files if k.startswith("constructs/")) == [
        f"constructs/binder_{x}.gb" for x in "ABC"
    ]
    for entry in made.entries[:3]:
        record = SeqIO.read(io.StringIO(files[f"constructs/{entry.name}.gb"]), "genbank")
        assert str(record.seq) == entry.assembly.seq and len(record.features) > 0
    assert files["failed.txt"].startswith("oops\t")


def test_a_clean_batch_writes_no_failure_file(vector, host, sites):
    ok = b.design_batch(b.read_proteins(FASTA), vector, *sites, host)
    assert "failed.txt" not in b.files(ok)


def test_write_batch_creates_the_files(made, tmp_path):
    paths = b.write_batch(made, tmp_path / "deep" / "out")
    assert {p.name for p in paths} >= {"summary.csv", "orders.fasta", "binder_A.gb", "failed.txt"}
    assert all(p.is_file() and tmp_path in p.parents for p in paths)


def test_every_ok_order_proves_it_clones(made):
    for e in made.entries[:3]:
        assert e.order.ligation.ok and "exactly" in e.order.ligation.message


def test_a_plain_enzyme_name_resolves_when_it_cuts_once(vector):
    assert resolve_site(vector, "xhoi").enzyme == "XhoI"
    assert resolve_site(vector, " NdeI ").enzyme == "NdeI"
    site = resolve_site(vector, "BamHI")
    assert resolve_site(vector, site.id) == site


def test_an_unknown_or_repeated_enzyme_is_explained(vector):
    with pytest.raises(ValueError, match="no Nope site"):
        resolve_site(vector, "Nope")
    seq = vector.seq + "GGATCC" + "ACGT" * 10
    twice = p2d.Vector.from_sequence(seq, name="twice")
    twice.expression_start = vector.expression_start
    with pytest.raises(ValueError, match="2 times; choose one of"):
        resolve_site(twice, "BamHI")


runner = CliRunner()


def test_the_cli_writes_a_batch_and_exits_zero(tmp_path):
    fa = tmp_path / "p.fa"
    fa.write_text(FASTA)
    r = runner.invoke(app, ["batch", str(fa), "--out", str(tmp_path / "o"), "--cterm", "stop"])
    assert r.exit_code == 0, r.output
    assert "3 of 3 designed" in r.output
    assert (tmp_path / "o" / "summary.csv").is_file()
    assert len(list((tmp_path / "o" / "constructs").glob("*.gb"))) == 3


def test_the_cli_exits_one_when_a_protein_fails_but_still_writes_the_rest(tmp_path):
    fa = tmp_path / "p.fa"
    fa.write_text(FASTA + ">bad\nMK#\n")
    r = runner.invoke(app, ["batch", str(fa), "-o", str(tmp_path / "o")])
    assert r.exit_code == 1 and "failed  bad" in r.output and "3 of 4" in r.output
    assert (tmp_path / "o" / "failed.txt").read_text().startswith("bad\t")


@pytest.mark.parametrize(
    "args,message",
    [
        (["--host", "nowhere"], "unknown host"),
        (["--upstream", "Nope"], "no Nope site"),
        (["--cterm", "sideways"], "sideways"),
        (["--n-part", "nothing"], "unknown tag"),
    ],
)
def test_the_cli_refuses_bad_options_with_exit_two(tmp_path, args, message):
    fa = tmp_path / "p.fa"
    fa.write_text(FASTA)
    r = runner.invoke(app, ["batch", str(fa), "-o", str(tmp_path / "o"), *args])
    assert r.exit_code == 2 and message.lower() in r.output.lower()
    assert not (tmp_path / "o").exists()


def test_the_cli_needs_a_start_unless_the_plasmid_gives_one(tmp_path, vector):
    gb = tmp_path / "v.gb"
    SeqIO.write(vector.record, str(gb), "genbank")
    fa = tmp_path / "p.fa"
    fa.write_text(FASTA)
    out = tmp_path / "o"
    base = ["batch", str(fa), "-o", str(out), "--vector", str(gb)]
    no_start = runner.invoke(app, [*base, "--start", str(vector.expression_start + 1)])
    assert no_start.exit_code == 0, no_start.output
    assert (out / "summary.csv").is_file()


def test_the_cli_reports_an_empty_or_missing_file(tmp_path):
    empty = tmp_path / "e.fa"
    empty.write_text("")
    r = runner.invoke(app, ["batch", str(empty), "-o", str(tmp_path / "o")])
    assert r.exit_code == 2 and "no protein records" in r.output
    r = runner.invoke(app, ["batch", str(tmp_path / "missing.fa"), "-o", str(tmp_path / "o")])
    assert r.exit_code == 2


def test_the_vector_bundled_sites_are_cut_once():
    v = p2d.load_bundled()
    assert {c.enzyme for c in list_cut_sites(v)} >= {"NdeI", "XhoI", "BamHI"}


def test_an_entry_that_fails_its_checks_is_never_written_as_something_to_order(made):
    import dataclasses

    from p2d.order import LigationCheck

    good = made.entries[0]
    broken_order = dataclasses.replace(good.order, ligation=LigationCheck(False, "does not clone"))
    broken = dataclasses.replace(good, order=broken_order)
    assert good.ok and not broken.ok
    result = b.BatchResult(made.vector, [broken])
    assert b.orders_fasta(result) == "" and b.inserts_fasta(result) == ""
    assert not any(k.startswith("constructs/") for k in b.files(result))
    assert b.summary_rows(result)[0]["status"] == "failed"


def test_a_batch_into_the_real_pet28a_from_the_command_line(tmp_path):
    import os
    from pathlib import Path

    samples = os.environ.get("P2D_VECTOR_SAMPLES")
    path = Path(samples) / "pET-28a_plus.dna" if samples else None
    if not path or not path.is_file():
        pytest.skip("pET-28a_plus.dna not available (set P2D_VECTOR_SAMPLES)")
    fa = tmp_path / "p.fa"
    fa.write_text(FASTA)
    r = runner.invoke(
        app,
        [
            "batch",
            str(fa),
            "-o",
            str(tmp_path / "o"),
            "--vector",
            str(path),
            "-u",
            "NdeI",
            "-d",
            "XhoI",
        ],
    )
    assert r.exit_code == 0, r.output
    rows = list(csv.DictReader((tmp_path / "o" / "summary.csv").open()))
    assert [x["status"] for x in rows] == ["ok"] * 3
    for x in rows:
        record = SeqIO.read(tmp_path / "o" / "constructs" / f"{x['name']}.gb", "genbank")
        assert len(record.seq) == int(x["plasmid_bp"])
        seq = str(record.seq)
        start = seq.find("ATGGGCAGCAGCCATCATCATCATCATCAC")  # the His6-thrombin start of pET-28a
        assert start >= 0
        protein = str(record.seq[start : start + 63].translate())
        assert protein.startswith("MGSSHHHHHHSSGLVPRGSHM")
