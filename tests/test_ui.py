"""The web UI layer: input handling, the JSON service, and the local server."""

import http.client
import json
import threading
from importlib import resources

import pytest
from Bio import SeqIO

from p2d_ui import service
from p2d_ui.server import make_server

PEPTIDE = "SKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLKFICTTGKLPVPWPTL"


# ------------------------------------------------------------ parse_protein


def test_parse_plain_and_fasta_and_numbered():
    assert service.parse_protein("mkta yiak") == "MKTAYIAK"
    assert service.parse_protein(">sp|X|NAME desc\nMKTA\nYIAK\n") == "MKTAYIAK"
    assert service.parse_protein("  1 MKTAYIAKQR\n 11 LQ") == "MKTAYIAKQRLQ"
    assert service.parse_protein("MKTA*") == "MKTA"


def test_parse_uses_only_the_first_fasta_record():
    assert service.parse_protein(">a\nMKTA\n>b\nGGGG\n") == "MKTA"


@pytest.mark.parametrize("text", ["", "   \n", ">only a header\n"])
def test_parse_rejects_empty(text):
    with pytest.raises(service.InputError, match="Paste a protein"):
        service.parse_protein(text)


def test_parse_rejects_unknown_residues_and_reports_them():
    with pytest.raises(service.InputError, match="B, Z"):
        service.parse_protein("MKBZA")


def test_parse_enforces_length_limit():
    with pytest.raises(service.InputError, match="limit"):
        service.parse_protein("A" * (service.MAX_RESIDUES + 1))
    assert len(service.parse_protein("A" * service.MAX_RESIDUES)) == service.MAX_RESIDUES


# ------------------------------------------------------------------ service


def test_meta_lists_vectors_hosts_strategies():
    m = service.meta()
    assert m["vectors"][0]["name"] == "pTEST1"
    assert {"NdeI-XhoI", "BamHI-XhoI"} <= {s["name"] for s in m["vectors"][0]["sites"]}
    assert any(h["slug"] == "ecoli_bl21" and h["verified"] for h in m["hosts"])
    assert [s["value"] for s in m["strategies"]] == ["weighted_sample", "max_cai"]


def test_design_returns_everything_the_ui_renders():
    r = service.design({"protein": PEPTIDE, "vector": "pTEST1", "site": "NdeI-XhoI"})
    assert r["ok"] is True
    assert r["payload_length"] == len(PEPTIDE)
    assert [s["start"] for s in r["segments"]][0] == 1
    assert r["segments"][-1]["end"] == len(r["protein"])
    payload_seg = next(s for s in r["segments"] if s["kind"] == "payload")
    assert payload_seg["aa"] == PEPTIDE  # the payload block, scar residues excluded
    assert 0 < r["metrics"]["cai"] <= 1
    a, b = r["coding_span"]
    assert b - a == 3 * len(PEPTIDE)
    assert r["insert_dna"][:6] == "CATATG"


def test_design_outputs_parse_back():
    r = service.design({"protein": PEPTIDE, "site": "NdeI-XhoI"})
    from io import StringIO

    record = SeqIO.read(StringIO(r["genbank"]), "genbank")
    assert len(record.seq) == r["plasmid_bp"]
    assert r["fasta"].startswith(">") and r["insert_dna"] in r["fasta"]


def test_issues_are_sorted_errors_first():
    r = service.design(
        {"protein": "M" + PEPTIDE, "site": "NcoI-XhoI"}
    )  # leading M duplicates the Met NcoI donates
    order = [i["severity"] for i in r["issues"]]
    assert any(i["code"] == "duplicate_met" for i in r["issues"])
    assert order == sorted(order, key=("error", "warning", "info").index)


def test_seed_and_strategy_are_honoured():
    base = {"protein": PEPTIDE * 2, "site": "NdeI-XhoI"}
    a = service.design({**base, "seed": 1})["insert_dna"]
    assert a == service.design({**base, "seed": "1"})["insert_dna"]
    assert a != service.design({**base, "seed": 2})["insert_dna"]
    best = service.design({**base, "strategy": "max_cai"})
    assert best["metrics"]["strategy"] == "max_cai"
    assert best["metrics"]["cai"] > 0.95


@pytest.mark.parametrize(
    "patch, message",
    [
        ({"vector": "nope"}, "Unknown vector"),
        ({"site": "NotAnEnzyme"}, "no insertion site"),
        ({"host": "yeast"}, "Unknown host"),
        ({"strategy": "harmonized"}, "Unknown strategy"),
        ({"seed": "abc"}, "whole number"),
        ({"seed": -1}, "between"),
        ({"protein": "MK1B"}, "Unsupported"),
    ],
)
def test_bad_input_is_a_clean_input_error(patch, message):
    req = {"protein": PEPTIDE, "vector": "pTEST1", "site": "NdeI-XhoI", **patch}
    with pytest.raises(service.InputError, match=message):
        service.design(req)


# ------------------------------------------------------------------- server


@pytest.fixture(scope="module")
def server():
    srv = make_server(port=0)
    thread = threading.Thread(target=srv.serve_forever, daemon=True)
    thread.start()
    yield srv
    srv.shutdown()
    srv.server_close()


def call(server, method, path, body=None, headers=None):
    port = server.server_port
    conn = http.client.HTTPConnection("127.0.0.1", port, timeout=10)
    h = {"Host": f"127.0.0.1:{port}", **(headers or {})}
    payload = None
    if body is not None:
        payload = body if isinstance(body, bytes) else json.dumps(body).encode()
        h.setdefault("Content-Type", "application/json")
    conn.request(method, path, body=payload, headers=h)
    resp = conn.getresponse()
    data = resp.read()
    conn.close()
    return resp, data


def test_serves_the_three_static_files_with_security_headers(server):
    for path, ctype in [
        ("/", "text/html"),
        ("/app.css", "text/css"),
        ("/app.js", "text/javascript"),
    ]:
        resp, data = call(server, "GET", path)
        assert resp.status == 200 and resp.getheader("Content-Type").startswith(ctype)
        assert data
        assert "default-src 'self'" in resp.getheader("Content-Security-Policy")
        assert resp.getheader("X-Content-Type-Options") == "nosniff"


def test_api_meta_and_design_round_trip(server):
    resp, data = call(server, "GET", "/api/meta")
    assert resp.status == 200 and json.loads(data)["vectors"]

    resp, data = call(server, "POST", "/api/design", {"protein": PEPTIDE, "site": "NdeI-XhoI"})
    assert resp.status == 200 and json.loads(data)["ok"] is True


def test_input_errors_are_422_with_a_message(server):
    resp, data = call(server, "POST", "/api/design", {"protein": "MKB", "site": "NdeI-XhoI"})
    assert resp.status == 422 and "Unsupported" in json.loads(data)["error"]


def test_malformed_requests_are_rejected(server):
    assert call(server, "POST", "/api/design", b"not json")[0].status == 400
    assert call(server, "POST", "/api/design", b"[1,2]")[0].status == 400
    resp, _ = call(server, "POST", "/api/design", b"x", {"Content-Type": "text/plain"})
    assert resp.status == 415
    # declare an oversized body without sending it; the server must refuse on the header alone
    resp, _ = call(server, "POST", "/api/design", b"x", {"Content-Length": "1000001"})
    assert resp.status == 413


def test_unknown_paths_and_traversal_are_404(server):
    for path in [
        "/nope",
        "/../pyproject.toml",
        "/static/../service.py",
        "/%2e%2e/service.py",
        "/app.py",
    ]:
        assert call(server, "GET", path)[0].status == 404, path
    assert call(server, "POST", "/api/meta", {})[0].status == 404


def test_foreign_host_header_is_refused(server):
    """DNS-rebinding guard: a hostile page resolving to 127.0.0.1 sends its own Host."""
    resp, _ = call(server, "GET", "/api/meta", headers={"Host": "evil.example"})
    assert resp.status == 403
    resp, _ = call(server, "POST", "/api/design", {"protein": PEPTIDE}, {"Host": "evil.example"})
    assert resp.status == 403


# -------------------------------------------------------------- static files


def static_text(name):
    return resources.files("p2d_ui").joinpath("static", name).read_text(encoding="utf-8")


@pytest.mark.parametrize("name", ["index.html", "app.css", "app.js"])
def test_no_em_or_en_dashes_in_ui_files(name):
    text = static_text(name)
    assert "—" not in text and "–" not in text


def test_ui_never_inserts_server_text_as_html():
    js = static_text("app.js")
    for banned in ["innerHTML", "outerHTML", "insertAdjacentHTML", "document.write", "eval("]:
        assert banned not in js


def test_example_sequence_is_valid_and_designs_cleanly():
    import re

    js = static_text("app.js")
    seq = re.search(r'GFP_EXAMPLE = "([A-Z]+)"', js).group(1)
    assert len(seq) == 237
    r = service.design({"protein": seq, "site": "NdeI-XhoI"})
    assert r["ok"] and not [i for i in r["issues"] if i["severity"] != "info"]


# --------------------------------------------------------------- enzyme picker


def test_meta_says_which_vectors_support_the_picker():
    assert service.meta()["vectors"][0]["picker"] is True


def test_enzymes_lists_single_cutters_in_order_with_positions():
    r = service.enzymes({"vector": "pTEST1"})
    names = [e["name"] for e in r["enzymes"]]
    assert names == ["NcoI", "NdeI", "BamHI", "Acc65I", "KpnI", "EcoRI", "XhoI"]
    assert r["window"] == {"start": 87, "stop_end": 201}
    assert [e["start"] for e in r["enzymes"]] == sorted(e["start"] for e in r["enzymes"])
    xho = r["enzymes"][-1]
    assert xho["aliases"] and xho["site"] == "CTCGAG" and r["protein_length"] is None


def test_upstream_verdict_says_where_the_start_codon_comes_from():
    rows = {e["name"]: e for e in service.enzymes({"protein": PEPTIDE})["enzymes"]}
    assert rows["NcoI"]["as_upstream"]["starts_at"] == "site"
    assert rows["NdeI"]["as_upstream"]["starts_at"] == "vector"
    assert rows["NdeI"]["as_upstream"]["ok"] is True
    # nothing usable lies downstream of the last enzymes
    assert rows["XhoI"]["as_upstream"]["ok"] is False
    assert "pair with it" in rows["XhoI"]["as_upstream"]["reason"]


def test_downstream_verdicts_follow_direction_and_end_compatibility():
    r = service.enzymes({"protein": PEPTIDE, "upstream": "NdeI"})
    down = {e["name"]: e["as_downstream"] for e in r["enzymes"]}
    assert (
        down["NcoI"]["hidden"] and down["NdeI"]["hidden"]
    )  # upstream of the choice: not candidates
    assert down["BamHI"]["ok"] is False and "touch" in down["BamHI"]["reason"]
    assert down["XhoI"]["ok"] is True and down["EcoRI"]["ok"] is True


def test_unusable_enzymes_come_with_a_reason_for_the_protein():
    assert service._protein_conflict("AMA", "ATG", "ecoli_bl21") == (
        "Your protein forces the ATG site near residue 2; no synonymous codons avoid it."
    )
    assert service._protein_conflict(PEPTIDE, "CTCGAG", "ecoli_bl21") is None


def test_unusable_protein_input_does_not_break_the_enzyme_list():
    for bad in ["", "MKB1Z", None]:
        r = service.enzymes({"protein": bad})
        assert r["protein_length"] is None and len(r["enzymes"]) == 7


def test_design_with_an_enzyme_pair_reports_the_real_protein():
    nde = service.design({"protein": PEPTIDE, "upstream": "NdeI", "downstream": "XhoI"})
    assert nde["ok"] and nde["site"] == "NdeI-XhoI"
    assert nde["protein"].startswith("MGSSHHHHHHSSGLVPRGSHM" + PEPTIDE)  # tag kept
    nco = service.design({"protein": PEPTIDE, "upstream": "NcoI", "downstream": "XhoI"})
    assert nco["ok"] and nco["protein"].startswith("MG" + PEPTIDE)  # tag removed


@pytest.mark.parametrize(
    "patch, message",
    [
        ({"upstream": "XhoI", "downstream": "NdeI"}, "must lie upstream"),
        ({"upstream": "NdeI", "downstream": "BamHI"}, "touch or overlap"),
        ({"upstream": "SmaI", "downstream": "XhoI"}, "not a single-cutter"),
        ({"upstream": "NdeI", "downstream": ""}, "both a 5' and a 3' enzyme"),
        ({"upstream": "", "downstream": "", "site": ""}, "Choose a 5' and a 3' enzyme"),
    ],
)
def test_bad_enzyme_choices_are_clean_input_errors(patch, message):
    with pytest.raises(service.InputError, match=message):
        service.design({"protein": PEPTIDE, "vector": "pTEST1", **patch})


def test_enzymes_route_works_and_shares_the_same_guards(server):
    resp, data = call(server, "POST", "/api/enzymes", {"protein": PEPTIDE, "upstream": "NdeI"})
    assert resp.status == 200 and len(json.loads(data)["enzymes"]) == 7
    resp, data = call(server, "POST", "/api/enzymes", {"vector": "nope"})
    assert resp.status == 422 and "Unknown vector" in json.loads(data)["error"]
    resp, _ = call(server, "POST", "/api/enzymes", {}, {"Host": "evil.example"})
    assert resp.status == 403
