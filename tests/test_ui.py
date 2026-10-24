"""The web UI layer: input handling, the JSON service, and the local server."""

import base64
import http.client
import json
import os
import re
import threading
from importlib import resources
from io import StringIO
from pathlib import Path

import pytest
from Bio import SeqIO
from Bio.Seq import Seq

import p2d
from p2d_ui import service
from p2d_ui.server import make_server
from synth import mirror, ptest1_features, snapgene_bytes

SAMPLES = Path(os.environ["P2D_VECTOR_SAMPLES"]) if os.environ.get("P2D_VECTOR_SAMPLES") else None
PEPTIDE = "SKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLKFICTTGKLPVPWPTL"


def b64(data: bytes) -> str:
    return base64.b64encode(data).decode()


@pytest.fixture
def pt1():
    """A freshly loaded bundled pTEST1 (id and start codon)."""
    r = service.load_vector({"bundled": "pTEST1"})
    return r["id"], r["starts"][0]["position"]


def req(pt1, **kw):
    vid, start = pt1
    return {"vector_id": vid, "start": start, "protein": PEPTIDE, **kw}


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


def test_meta_lists_bundled_library_hosts_and_strategies(tmp_path, monkeypatch):
    monkeypatch.setenv("P2D_VECTOR_DIR", str(tmp_path))
    (tmp_path / "mine.gb").write_text("x")
    (tmp_path / "notes.pdf").write_text("x")
    m = service.meta()
    assert m["bundled"] == [{"id": "pTEST1", "name": "pTEST1", "length": 1149}]
    assert m["library"] == ["mine.gb"]
    assert any(h["slug"] == "ecoli_bl21" and h["verified"] for h in m["hosts"])
    assert [s["value"] for s in m["strategies"]] == ["weighted_sample", "max_cai"]


def test_bundled_vector_loads_with_its_analysis():
    r = service.load_vector({"bundled": "pTEST1"})
    assert (r["name"], r["length"], r["circular"], r["flipped"]) == ("pTEST1", 1149, True, False)
    best = r["starts"][0]
    assert best["position"] == 87 and best["confidence"] == "high" and best["evidence"]
    assert best["n_terminal"].startswith("MGSSHHHHHHSSGLVPRGS")
    assert r["site_count"] == 9 and r["single_cutters"] == 9
    labels = {f["label"] for f in r["features"]}
    assert {"start codon", "KanR (placeholder)"} <= labels
    assert all(f["type"] != "source" for f in r["features"])


def test_uploaded_snapgene_file_and_its_reverse_give_the_same_plasmid():
    v = p2d.load_bundled()
    fwd = snapgene_bytes(v.seq, ptest1_features(v))
    rc = snapgene_bytes(
        str(Seq(v.seq).reverse_complement()), mirror(ptest1_features(v), len(v.seq))
    )
    a = service.load_vector({"filename": "a.dna", "data_b64": b64(fwd)})
    b = service.load_vector({"filename": "b.dna", "data_b64": b64(rc)})
    assert (a["flipped"], b["flipped"]) == (False, True)
    assert a["starts"][0]["position"] == b["starts"][0]["position"] == 87
    assert a["site_count"] == b["site_count"] and a["mcs"] == b["mcs"] == {"start": 144, "end": 180}
    assert any("Reverse-complemented" in n for n in b["notes"])
    assert a["source"] == "uploaded file"


def test_pasted_text_loads_and_is_labelled():
    seq = p2d.load_bundled().seq
    r = service.load_vector({"text": f">my plasmid\n{seq}\n", "filename": "mine"})
    assert r["source"] == "pasted text" and r["length"] == 1149
    assert r["starts"][0]["position"] == 87 and r["starts"][0]["confidence"] == "high"
    assert not any("not certain" in w for w in r["warnings"])


def test_the_users_own_vector_folder_is_offered_but_only_by_listed_name(tmp_path, monkeypatch):
    monkeypatch.setenv("P2D_VECTOR_DIR", str(tmp_path))
    text = resources.files("p2d.data").joinpath("pTEST1.gb").read_text()
    (tmp_path / "lab_plasmid.gb").write_text(text)
    assert service.load_vector({"library": "lab_plasmid.gb"})["source"] == "your vector folder"
    for bad in ["../lab_plasmid.gb", "/etc/passwd", "nope.gb", "lab_plasmid.gb/../x"]:
        with pytest.raises(service.InputError, match="not in your vector folder"):
            service.load_vector({"library": bad})


@pytest.mark.parametrize(
    "request_, message",
    [
        ({}, "Choose a plasmid"),
        ({"bundled": "pXYZ"}, "Unknown bundled"),
        ({"data_b64": "!!!not base64!!!"}, "could not be decoded"),
        ({"data_b64": b64(b"\x00\x01\x02" * 400), "filename": "x.bin"}, "binary"),
        (
            {"text": ">p\nMKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALPDAQ\n"},
            "protein",
        ),
        ({"text": "ACGT"}, "bp"),
    ],
)
def test_unusable_vectors_are_clean_input_errors(request_, message):
    with pytest.raises(service.InputError, match=message):
        service.load_vector(request_)


def test_a_plasmid_that_is_no_longer_loaded_says_so():
    with pytest.raises(service.InputError, match="no longer loaded"):
        service.sites({"vector_id": "gone"})


def test_old_plasmids_are_forgotten_but_recent_ones_survive():
    first = service.load_vector({"bundled": "pTEST1"})["id"]
    for _ in range(service.MAX_STORED + 2):
        service.load_vector({"bundled": "pTEST1"})
    with pytest.raises(service.InputError, match="no longer loaded"):
        service.sites({"vector_id": first})


def test_sites_are_listed_by_position_with_a_default_window(pt1):
    r = service.sites(req(pt1))
    assert [s["id"] for s in r["sites"]] == [
        "NcoI:85",
        "NdeI:144",
        "BamHI:150",
        "Acc65I:156",
        "KpnI:156",
        "EcoRI:168",
        "XhoI:174",
    ], "NcoI holds the start codon, so the window starts a little before it"
    assert r["start"] == 87 and r["window"]["start"] == 77
    ncoi = r["sites"][0]
    assert (ncoi["site"], ncoi["overhang"], ncoi["count"], ncoi["number"]) == (
        "CCATGG",
        "CATG",
        1,
        1,
    )


def test_five_prime_verdicts_say_why_not(pt1):
    rows = {s["id"]: s["as_upstream"] for s in service.sites(req(pt1))["sites"]}
    assert rows["NdeI:144"]["ok"] and rows["BamHI:150"]["ok"] and rows["NcoI:85"]["ok"]
    assert not rows["XhoI:174"]["ok"] and "downstream" in rows["XhoI:174"]["reason"]


def test_without_a_confirmed_start_nothing_is_usable_yet(pt1):
    r = service.sites({"vector_id": pt1[0]})
    assert r["start"] is None and r["window"] == {"start": 0, "end": 1149}
    assert all(not s["as_upstream"]["ok"] for s in r["sites"])
    assert "Confirm the start codon" in r["sites"][0]["as_upstream"]["reason"]


def test_three_prime_verdicts_follow_direction_spacing_and_tag(pt1):
    r = service.sites(req(pt1, upstream="BamHI:150"))
    down = {s["id"]: s["as_downstream"] for s in r["sites"]}
    assert down["NdeI:144"]["hidden"] and down["NcoI:85"]["hidden"], "wrong side: not candidates"
    assert down["KpnI:156"]["ok"], "touching sites are allowed (with a warning)"
    assert down["EcoRI:168"]["ok"] and down["XhoI:174"]["ok"]
    assert down["XhoI:174"]["has_tag"] is True


def test_a_site_the_protein_forces_is_unusable_with_the_reason(pt1, monkeypatch):
    monkeypatch.setattr(
        service, "_protein_conflict", lambda p, site, h: "forced" if site == "GAATTC" else None
    )
    r = service.sites(req(pt1, upstream="BamHI:150"))
    eco = next(s for s in r["sites"] if s["id"] == "EcoRI:168")
    assert eco["as_downstream"] == {
        "ok": False,
        "reason": "forced",
        "has_tag": True,
        "tag_name": "His6",
    }


def test_protein_conflict_message_names_the_residue():
    assert service._protein_conflict("AMA", "ATG", "ecoli_bl21") == (
        "Your protein forces the ATG site near residue 2; no synonymous codons avoid it."
    )
    assert service._protein_conflict(PEPTIDE, "CTCGAG", "ecoli_bl21") is None


def test_multi_cutters_can_be_hidden_or_shown():
    seq = p2d.load_bundled().seq
    twice = seq[:400] + "CTCGAG" + seq[406:]  # a second XhoI in the backbone
    r = service.load_vector({"text": twice, "filename": "twice"})
    args = {"vector_id": r["id"], "start": 87, "window": [0, 600]}
    xhos = [s for s in service.sites(args)["sites"] if s["enzyme"] == "XhoI"]
    assert [(s["number"], s["count"]) for s in xhos] == [(1, 2), (2, 2)]
    once = service.sites({**args, "include_multi": False})["sites"]
    assert all(s["count"] == 1 for s in once) and not any(s["enzyme"] == "XhoI" for s in once)


def test_the_window_can_be_chosen_and_is_validated(pt1):
    r = service.sites(req(pt1, window=[149, 157]))
    assert [s["id"] for s in r["sites"]] == ["BamHI:150"]
    for bad, message in [
        ([5, 5], "end after"),
        (["a", 9], "whole numbers"),
        ([1], "whole numbers"),
    ]:
        with pytest.raises(service.InputError, match=message):
            service.sites(req(pt1, window=bad))


@pytest.mark.parametrize(
    "start, message", [("x", "whole number"), (5, "not an ATG"), (10**6, "not an ATG")]
)
def test_a_start_that_is_not_an_atg_is_refused(pt1, start, message):
    with pytest.raises(service.InputError, match=message):
        service.sites({"vector_id": pt1[0], "start": start})


def test_design_returns_everything_the_ui_renders(pt1):
    r = service.design(req(pt1, upstream="NdeI:144", downstream="XhoI:174"))
    assert r["ok"] and r["upstream"] == "NdeI at 145" and r["downstream"] == "XhoI at 175"
    payload_seg = next(s for s in r["segments"] if s["kind"] == "payload")
    assert payload_seg["aa"] == PEPTIDE and r["segments"][-1]["end"] == len(r["protein"])
    a, b = r["coding_span"]
    assert b - a == 3 * len(PEPTIDE) and r["insert_dna"][:6] == "CATATG"
    assert 0 < r["metrics"]["cai"] <= 1
    sizes = r["report"]["sizes"]
    assert sizes["backbone_bp"] + sizes["insert_after_digest_bp"] == r["plasmid_bp"]
    assert sizes["plasmid_bp"] == r["plasmid_bp"]
    assert r["report"]["protein"]["residues"] == len(r["protein"])


def test_the_genbank_output_keeps_the_vectors_annotations(pt1):
    r = service.design(req(pt1, upstream="NdeI:144", downstream="XhoI:174"))
    rec = SeqIO.read(StringIO(r["genbank"]), "genbank")
    assert len(rec.seq) == r["plasmid_bp"]
    labels = {f.qualifiers.get("label", [""])[0]: f for f in rec.features}
    assert {"T7 promoter / lac operator / RBS", "start codon", "KanR (placeholder)"} <= set(labels)
    assert {"insert (p2d)", "fusion ORF"} <= set(labels)
    assert "NdeI" not in labels, "features inside the replaced stretch are dropped"
    kan = labels["KanR (placeholder)"]
    shift = r["plasmid_bp"] - 1149
    assert (int(kan.location.start), int(kan.location.end)) == (669 + shift, 1149 + shift)
    assert r["fasta"].startswith(">")


def test_c_terminus_choice_changes_the_protein_and_the_insert(pt1):
    stop = service.design(req(pt1, upstream="NdeI:144", downstream="XhoI:174", cterm="stop"))
    tag = service.design(req(pt1, upstream="NdeI:144", downstream="XhoI:174", cterm="tag"))
    assert stop["protein"].endswith(PEPTIDE) and stop["insert_dna"].endswith("TAACTCGAG")
    assert (
        tag["protein"].endswith(PEPTIDE + "LEHHHHHH") and "TAA" + "CTCGAG" not in tag["insert_dna"]
    )
    assert stop["cterm"] == "stop" and tag["cterm"] == "tag"


def test_seed_strategy_and_reproducibility(pt1):
    base = req(pt1, upstream="NdeI:144", downstream="XhoI:174", protein=PEPTIDE * 2)
    a = service.design({**base, "seed": 1})["insert_dna"]
    assert a == service.design({**base, "seed": "1"})["insert_dna"]
    assert a != service.design({**base, "seed": 2})["insert_dna"]
    best = service.design({**base, "strategy": "max_cai"})
    assert best["metrics"]["strategy"] == "max_cai" and best["metrics"]["cai"] > 0.95


@pytest.mark.parametrize(
    "patch, message",
    [
        ({"upstream": "XhoI:174", "downstream": "NdeI:144"}, "must lie upstream"),
        ({"upstream": "SmaI:100", "downstream": "XhoI:174"}, "No SmaI site"),
        ({"upstream": "NdeI:144", "downstream": ""}, "both a 5' and a 3' site"),
        ({"upstream": "", "downstream": ""}, "both a 5' and a 3' site"),
        (
            {"upstream": "NdeI:144", "downstream": "XhoI:174", "cterm": "weird"},
            "Unknown C-terminus",
        ),
        ({"upstream": "NdeI:144", "downstream": "XhoI:174", "strategy": "x"}, "Unknown strategy"),
        ({"upstream": "NdeI:144", "downstream": "XhoI:174", "host": "yeast"}, "Unknown host"),
        ({"upstream": "NdeI:144", "downstream": "XhoI:174", "seed": "abc"}, "whole number"),
        ({"upstream": "NdeI:144", "downstream": "XhoI:174", "seed": -1}, "between"),
        ({"upstream": "NdeI:144", "downstream": "XhoI:174", "protein": "MK1B"}, "Unsupported"),
    ],
)
def test_bad_choices_are_clean_input_errors(pt1, patch, message):
    with pytest.raises(service.InputError, match=message):
        service.design(req(pt1, **patch))


def test_a_pair_that_would_cut_the_backbone_is_refused_with_the_position():
    seq = p2d.load_bundled().seq
    r = service.load_vector({"text": seq[:700] + "GGATCC" + seq[706:], "filename": "dup"})
    args = {"vector_id": r["id"], "start": 87, "protein": PEPTIDE, "upstream": "BamHI:150"}
    with pytest.raises(service.InputError, match=r"BamHI at 701"):
        service.design({**args, "downstream": "XhoI:174"})


def test_a_reverse_stored_upload_designs_the_same_construct():
    v = p2d.load_bundled()
    fwd = service.load_vector(
        {"filename": "f.dna", "data_b64": b64(snapgene_bytes(v.seq, ptest1_features(v)))}
    )
    rc = snapgene_bytes(str(Seq(v.seq).reverse_complement()), mirror(ptest1_features(v), 1149))
    rev = service.load_vector({"filename": "r.dna", "data_b64": b64(rc)})
    out = [
        service.design(
            {
                "vector_id": r["id"],
                "start": r["starts"][0]["position"],
                "protein": PEPTIDE,
                "upstream": "BamHI:150",
                "downstream": "XhoI:174",
                "seed": 4,
                "cterm": "tag",
            }
        )
        for r in (fwd, rev)
    ]
    assert out[0]["insert_dna"] == out[1]["insert_dna"] and out[0]["protein"] == out[1]["protein"]
    assert out[0]["report"]["sizes"] == out[1]["report"]["sizes"]


def test_real_pet28a_through_the_service():
    if SAMPLES is None or not (SAMPLES / "pET-28a_plus.dna").is_file():
        pytest.skip("pET-28a_plus.dna not available (set P2D_VECTOR_SAMPLES)")
    data = (SAMPLES / "pET-28a_plus.dna").read_bytes()
    r = service.load_vector({"filename": "pET-28a(+).dna", "data_b64": b64(data)})
    assert r["flipped"] and r["starts"][0]["confidence"] == "high"
    assert r["site_count"] == 60 and r["single_cutters"] == 25
    start = r["starts"][0]["position"]
    listing = service.sites(
        {"vector_id": r["id"], "start": start, "protein": PEPTIDE, "upstream": "BamHI:5166"}
    )
    down = {s["id"]: s["as_downstream"] for s in listing["sites"]}
    assert down["XhoI:5206"]["ok"] and down["XhoI:5206"]["has_tag"] is True
    out = {}
    for mode in ("tag", "stop"):
        out[mode] = service.design(
            {
                "vector_id": r["id"],
                "start": start,
                "protein": PEPTIDE,
                "upstream": "BamHI:5166",
                "downstream": "XhoI:5206",
                "cterm": mode,
                "seed": 2,
            }
        )
        assert out[mode]["ok"], out[mode]["issues"]
    assert out["tag"]["protein"].endswith(PEPTIDE + "LEHHHHHH")
    assert out["stop"]["protein"].endswith(PEPTIDE)
    rec = SeqIO.read(StringIO(out["tag"]["genbank"]), "genbank")
    wanted = {"KanR", "lacI", "T7 promoter", "f1 ori"}
    assert wanted <= {f.qualifiers.get("label", [""])[0] for f in rec.features}


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


def test_the_whole_flow_over_http(server):
    resp, data = call(server, "GET", "/api/meta")
    assert resp.status == 200 and json.loads(data)["bundled"]
    resp, data = call(server, "POST", "/api/vector", {"bundled": "pTEST1"})
    loaded = json.loads(data)
    assert resp.status == 200 and loaded["starts"][0]["position"] == 87
    body = {"vector_id": loaded["id"], "start": 87, "protein": PEPTIDE, "upstream": "BamHI:150"}
    resp, data = call(server, "POST", "/api/sites", body)
    assert resp.status == 200 and len(json.loads(data)["sites"]) == 7
    resp, data = call(server, "POST", "/api/design", {**body, "downstream": "XhoI:174"})
    assert resp.status == 200 and json.loads(data)["ok"] is True


def test_the_vector_route_takes_a_body_larger_than_the_default_json_limit(server):
    seq = "ACGT" * 700_000  # 2.8 MB of sequence, ~3.7 MB once base64-encoded
    big = b64(f">p\n{seq}\n".encode())
    resp, data = call(server, "POST", "/api/vector", {"filename": "big.fa", "data_b64": big})
    # past the body-size gate (a 413 would mean it was cut off), then refused as too big a plasmid
    assert resp.status == 422 and "250000" in json.loads(data)["error"]


def test_input_errors_are_422_with_a_message(server):
    resp, data = call(server, "POST", "/api/vector", {"bundled": "nope"})
    assert resp.status == 422 and "Unknown bundled" in json.loads(data)["error"]


def test_malformed_requests_are_rejected(server):
    assert call(server, "POST", "/api/design", b"not json")[0].status == 400
    assert call(server, "POST", "/api/design", b"[1,2]")[0].status == 400
    resp, _ = call(server, "POST", "/api/design", b"x", {"Content-Type": "text/plain"})
    assert resp.status == 415
    # declare an oversized body without sending it; the server must refuse on the header alone
    resp, _ = call(server, "POST", "/api/design", b"x", {"Content-Length": "1000001"})
    assert resp.status == 413
    resp, _ = call(server, "POST", "/api/vector", b"x", {"Content-Length": "7000001"})
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
    assert call(server, "POST", "/api/enzymes", {})[0].status == 404, (
        "the old name-based route is gone"
    )


def test_foreign_host_header_is_refused(server):
    """DNS-rebinding guard: a hostile page resolving to 127.0.0.1 sends its own Host."""
    resp, _ = call(server, "GET", "/api/meta", headers={"Host": "evil.example"})
    assert resp.status == 403
    for path in ("/api/vector", "/api/sites", "/api/design"):
        resp, _ = call(server, "POST", path, {"bundled": "pTEST1"}, {"Host": "evil.example"})
        assert resp.status == 403


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


def test_example_sequence_is_valid_and_designs_cleanly(pt1):
    import re

    js = static_text("app.js")
    seq = re.search(r'GFP_EXAMPLE = "([A-Z]+)"', js).group(1)
    assert len(seq) == 237
    r = service.design(req(pt1, protein=seq, upstream="NdeI:144", downstream="XhoI:174"))
    assert r["ok"] and not [i for i in r["issues"] if i["severity"] != "info"]


def test_the_result_carries_the_dna_to_order(pt1):
    r = service.design(req(pt1, upstream="NdeI:144", downstream="XhoI:174"))
    o = r["order"]
    assert (len(o["clamp_5"]), len(o["clamp_3"])) == (6, 6), "six extra bases each side by default"
    assert o["dna"] == o["clamp_5"] + r["insert_dna"] + o["clamp_3"] and o["length"] == len(
        o["dna"]
    )
    a, b = o["coding_span"]
    assert b - a == 3 * len(PEPTIDE) and o["dna"][a : a + 3] != ""
    assert o["insert_span"] == [6, 6 + len(r["insert_dna"])]
    assert o["ligation"]["ok"] and "exactly" in o["ligation"]["message"]
    assert 0 < o["gc"] < 1 and o["longest_homopolymer"] >= 1
    assert r["fasta"] == f">{r['filename']}_order\n{o['dna']}\n"
    assert r["insert_fasta"] == f">{r['filename']}_insert\n{r['insert_dna']}\n"


@pytest.mark.parametrize("c5, c3", [(0, 0), (3, 12), ("4", "0"), (30, 30)])
def test_extra_bases_follow_the_request(pt1, c5, c3):
    r = service.design(
        req(pt1, upstream="BamHI:150", downstream="XhoI:174", clamp_5=c5, clamp_3=c3)
    )
    o = r["order"]
    assert (len(o["clamp_5"]), len(o["clamp_3"])) == (int(c5), int(c3))
    assert o["ligation"]["ok"] and r["ok"]


def test_blank_extra_bases_mean_the_default(pt1):
    r = service.design(
        req(pt1, upstream="NdeI:144", downstream="XhoI:174", clamp_5="", clamp_3=None)
    )
    assert (len(r["order"]["clamp_5"]), len(r["order"]["clamp_3"])) == (6, 6)


@pytest.mark.parametrize(
    "patch, message",
    [
        ({"clamp_5": -1}, "between 0 and 30"),
        ({"clamp_3": 31}, "between 0 and 30"),
        ({"clamp_5": "abc"}, "whole numbers"),
        ({"clamp_3": 2.5}, "whole numbers"),
        ({"clamp_5": True}, "whole numbers"),
    ],
)
def test_bad_extra_bases_are_clean_input_errors(pt1, patch, message):
    with pytest.raises(service.InputError, match=message):
        service.design(req(pt1, upstream="NdeI:144", downstream="XhoI:174", **patch))


def test_the_page_asks_for_the_extra_bases():
    html = static_text("index.html")
    for needle in ('id="clamp-5"', 'id="clamp-3"', 'value="6"', "Extra bases for ordering"):
        assert needle in html
    js = static_text("app.js")
    assert "clamp_5" in js and "DNA to order" in js


def test_pairs_group_workable_pairs_by_the_protein_they_give(pt1):
    r = service.pairs(req(pt1))
    assert r["workable"] == r["considered"] == 20  # touching pairs count too
    assert r["total_groups"] == len(r["groups"]) == 5
    first = r["groups"][0]
    assert (first["before"], first["after"], first["n_tag"]) == ("MG", "", False)
    assert first["n_pairs"] == 6 and len(first["pairs"]) == 6
    row = first["pairs"][0]
    assert (row["up"], row["up_label"]) == ("NcoI:85", "NcoI at 86")
    assert set(row) >= {"down", "down_label", "single_cutters", "replaced_bp", "plasmid_bp"}
    counted = sum(g["n_pairs"] for g in r["groups"])
    assert counted == r["workable"], "every pair is in exactly one group"


def test_a_suggested_pair_designs_to_exactly_what_was_promised(pt1):
    r = service.pairs(req(pt1))
    for g in r["groups"]:
        p = g["pairs"][0]
        d = service.design(req(pt1, upstream=p["up"], downstream=p["down"], cterm="stop"))
        assert d["ok"] and d["plasmid_bp"] == p["plasmid_bp"]
        assert d["protein"].startswith(g["before"] + PEPTIDE) and d["protein"].endswith(PEPTIDE)


def test_the_tag_filter_and_the_c_terminus_reach_the_suggestions(pt1):
    keep = service.pairs(req(pt1, n_tag="keep"))["groups"]
    drop = service.pairs(req(pt1, n_tag="remove"))["groups"]
    assert keep and all(g["n_tag"] for g in keep) and drop and not any(g["n_tag"] for g in drop)
    tag = service.pairs(req(pt1, cterm="tag"))["groups"]
    assert all(g["after"].endswith("HHHHHH") for g in tag)


def test_pairs_input_errors_are_clean(pt1):
    for patch, message in [
        ({"n_tag": "maybe"}, "keep or remove"),
        ({"cterm": "x"}, "Unknown C-terminus"),
        ({"host": "yeast"}, "Unknown host"),
        ({"protein": ""}, "Paste a protein"),
        ({"window": [5, 5]}, "end after it starts"),
    ]:
        with pytest.raises(service.InputError, match=message):
            service.pairs(req(pt1, **patch))
    with pytest.raises(service.InputError, match="Choose the start codon"):
        service.pairs({"vector_id": pt1[0], "protein": PEPTIDE})


def test_too_wide_a_window_asks_to_be_narrowed(pt1):
    seq = p2d.load_bundled().seq
    many = seq[:300] + ("GAATTCGGATCCCTCGAGAAGCTTGCTAGC" * 25) + seq[300:]
    r = service.load_vector({"text": many, "filename": "many"})
    with pytest.raises(service.InputError, match="narrow the window"):
        service.pairs(
            {"vector_id": r["id"], "start": 87, "protein": PEPTIDE, "window": [0, len(many)]}
        )


def test_pairs_over_http_and_the_page_has_the_button(server):
    resp, data = call(server, "POST", "/api/vector", {"bundled": "pTEST1"})
    vid = json.loads(data)["id"]
    resp, data = call(
        server, "POST", "/api/pairs", {"vector_id": vid, "start": 87, "protein": PEPTIDE}
    )
    assert resp.status == 200 and json.loads(data)["groups"]
    resp, data = call(server, "POST", "/api/pairs", {"vector_id": vid, "start": 87})
    assert resp.status == 422
    resp, _ = call(server, "POST", "/api/pairs", {}, {"Host": "evil.example"})
    assert resp.status == 403
    assert 'id="suggest-btn"' in static_text("index.html") and "/api/pairs" in static_text("app.js")


def test_meta_lists_the_tags_and_linkers():
    parts = service.meta()["parts"]
    assert {"his6", "flag", "v5"} <= {p["key"] for p in parts["tags"]}
    assert {"gs", "g4s"} <= {p["key"] for p in parts["linkers"]}
    assert all(p["aa"] for rows in parts.values() for p in rows)


def test_design_builds_the_chosen_tags_and_linkers_into_the_protein(pt1):
    d = service.design(
        req(
            pt1,
            upstream="NcoI:85",
            downstream="XhoI:174",
            n_parts=["flag", "gs"],
            c_parts=["strep2"],
        )
    )
    assert d["ok"] and "DYKDDDDKGS" + PEPTIDE + "WSHPQFEK" in d["protein"]


def test_the_suggestions_describe_only_what_the_vector_adds(pt1):
    plain = service.pairs(req(pt1))["groups"]
    built = service.pairs(req(pt1, n_parts=["flag"], c_parts=["gs", "ha"]))["groups"]
    assert [(g["before"], g["after"]) for g in plain] == [(g["before"], g["after"]) for g in built]


def test_unknown_or_malformed_parts_are_refused(pt1):
    for bad in (["nope"], "his6", ["his6"] * 5, [1]):
        for fn in (service.design, service.pairs):
            with pytest.raises(service.InputError):
                fn(req(pt1, upstream="NcoI:85", downstream="XhoI:174", n_parts=bad))


def test_the_page_has_the_tag_and_linker_selects():
    html = static_text("index.html")
    assert all(f'id="{i}"' in html for i in ("n-tag", "n-linker", "c-linker", "c-tag"))


def test_meta_gives_each_host_its_limits_in_percent():
    meta = service.meta()
    assert meta["default_host"] == "ecoli_bl21"
    assert {"hsapiens", "scerevisiae", "bsubtilis"} <= {h["slug"] for h in meta["hosts"]}
    host = next(h for h in meta["hosts"] if h["slug"] == "ecoli_bl21")
    assert host["limits"] == {
        "gc_min": 32,
        "gc_max": 60,
        "window_size": 50,
        "window_min": 28,
        "window_max": 70,
        "repeat_max": 22,
    }


def test_the_limits_chosen_reach_the_design(pt1):
    base = dict(upstream="NcoI:85", downstream="XhoI:174")
    tight = service.design(req(pt1, **base, limits={"gc_min": 55, "gc_max": 60}))
    assert tight["ok"] and 0.55 <= tight["metrics"]["gc"] <= 0.60
    assert tight["metrics"]["window_size"] == 50 and tight["metrics"]["gc_window_low"] is not None
    open_ = service.design(req(pt1, **base, limits={"gc_min": 0, "gc_max": 100}))
    assert open_["metrics"]["gc"] != tight["metrics"]["gc"]


def test_bad_limits_are_refused(pt1):
    base = dict(upstream="NcoI:85", downstream="XhoI:174")
    for limits, message in [
        ({"gc_min": 70, "gc_max": 30}, "above the upper"),
        ({"window_size": 3}, "between 10 and 500"),
        ({"nope": 1}, "Unknown limit"),
        ({"gc_min": "lots"}, "must be a number"),
        ({"gc_min": True}, "must be a number"),
        ("gc", "named numbers"),
    ]:
        for call, extra in ((service.design, base), (service.pairs, {})):
            with pytest.raises(service.InputError, match=message):
                call(req(pt1, **extra, limits=limits))


def test_the_page_has_the_limit_inputs():
    html = static_text("index.html")
    for i in ("gc-min", "gc-max", "window-size", "window-min", "window-max", "repeat-max"):
        assert f'id="lim-{i}"' in html


def test_every_host_can_design_through_the_service(pt1):
    names = {}
    for slug in service.meta()["hosts"]:
        d = service.design(req(pt1, upstream="NcoI:85", downstream="XhoI:174", host=slug["slug"]))
        assert d["ok"], slug["slug"]
        names[slug["slug"]] = d["insert_fasta"]
    assert len(set(names.values())) == len(names), "each host gives its own codons"


def test_each_host_note_reaches_the_page_data():
    hosts = {h["slug"]: h for h in service.meta()["hosts"]}
    assert "K-12" in hosts["ecoli_bl21"]["note"]
    assert "not for a cell line" in hosts["hsapiens"]["note"]


BATCH_FASTA = (
    ">first one\nKVFLDWINEAYQRGTRVLAEMAKRGDEFVKRLIAEGHDPFEVLKELGYSE\n"
    ">second\nMKVLAAGIVGLLL\n>oops\nMK#\n"
)


def batch_req(pt1, **extra):
    return req(pt1, upstream="NdeI:144", downstream="XhoI:174", fasta=BATCH_FASTA, **extra)


def test_batch_designs_every_protein_and_reports_the_one_that_fails(pt1):
    r = service.batch(batch_req(pt1))
    assert (r["n"], r["n_ok"]) == (3, 2)
    assert [(x["name"], x["status"]) for x in r["rows"]] == [
        ("first", "ok"),
        ("second", "ok"),
        ("oops", "failed"),
    ]
    assert all("order_dna" not in x for x in r["rows"]), "the DNA is in the downloads"
    assert r["rows"][2]["error"] and r["filename"] == "pTEST1_NdeI-XhoI_batch"
    assert r["orders_fasta"].count(">") == 2 and "first_order" in r["orders_fasta"]


def test_batch_zip_holds_the_summary_orders_and_constructs(pt1):
    import base64
    import io
    import zipfile

    r = service.batch(batch_req(pt1))
    z = zipfile.ZipFile(io.BytesIO(base64.b64decode(r["zip_b64"])))
    assert z.testzip() is None
    assert sorted(z.namelist()) == [
        "constructs/first.gb",
        "constructs/second.gb",
        "failed.txt",
        "inserts.fasta",
        "orders.fasta",
        "summary.csv",
    ]
    assert z.read("summary.csv").decode() == r["csv"]
    assert service.batch(batch_req(pt1))["zip_b64"] == r["zip_b64"], "same bytes every time"
    assert {i.date_time for i in z.infolist()} == {(2020, 1, 1, 0, 0, 0)}


def test_batch_uses_the_same_options_as_a_single_design(pt1):
    first_protein = BATCH_FASTA.split("\n")[1]
    one = service.design(
        req(
            pt1,
            upstream="NdeI:144",
            downstream="XhoI:174",
            n_parts=["flag"],
            clamp_5=9,
            protein=first_protein,
        )
    )
    many = service.batch(batch_req(pt1, n_parts=["flag"], clamp_5=9, seed=5))
    first = next(x for x in many["rows"] if x["name"] == "first")
    assert first["fusion_residues"] == len(one["protein"]) and first["status"] == "ok"
    for bad in ({"clamp_5": 99}, {"host": "nowhere"}, {"cterm": "x"}, {"n_parts": ["nope"]}):
        with pytest.raises(service.InputError):
            service.batch(batch_req(pt1, **bad))


def test_batch_input_errors_are_clean(pt1):
    for fasta, message in [
        ("", "at least one protein"),
        (None, "at least one protein"),
        (5, "at least one protein"),
        (">\n", "No protein records|no sequence"),
    ]:
        try:
            r = service.batch(req(pt1, upstream="NdeI:144", downstream="XhoI:174", fasta=fasta))
        except service.InputError as exc:
            assert re.search(message, str(exc))
        else:
            assert r["n_ok"] == 0 and r["rows"][0]["error"]
    many = "".join(f">p{i}\nMKV\n" for i in range(service.MAX_BATCH + 1))
    with pytest.raises(service.InputError, match="too many"):
        service.batch(req(pt1, upstream="NdeI:144", downstream="XhoI:174", fasta=many))
    with pytest.raises(service.InputError, match="5' and a 3' site"):
        service.batch(req(pt1, fasta=BATCH_FASTA))


def test_batch_over_http_and_the_page_has_the_panel(server):
    resp, data = call(server, "POST", "/api/vector", {"bundled": "pTEST1"})
    vid = json.loads(data)["id"]
    body = {
        "vector_id": vid,
        "start": 87,
        "upstream": "NdeI:144",
        "downstream": "XhoI:174",
        "fasta": BATCH_FASTA,
    }
    resp, data = call(server, "POST", "/api/batch", body)
    assert resp.status == 200 and json.loads(data)["n_ok"] == 2
    resp, _ = call(server, "POST", "/api/batch", {}, {"Host": "evil.example"})
    assert resp.status == 403
    html = static_text("index.html")
    assert all(
        f'id="{i}"' in html for i in ("batch-box", "batch-input", "batch-btn", "batch-result")
    )
    assert "/api/batch" in static_text("app.js")


def test_released_hosts_come_first_and_the_others_are_marked():
    meta = service.meta()
    status = {h["slug"]: h["status"] for h in meta["hosts"]}
    assert status["ecoli_bl21"] == status["bsubtilis"] == "released"
    assert all(status[s] == "development" for s in ("hsapiens", "mmusculus", "scerevisiae"))
    order = [h["status"] for h in meta["hosts"]]
    assert order == sorted(order, key=lambda s: s != "released"), "released hosts are listed first"
    assert meta["default_host"] == "ecoli_bl21"


def test_designing_for_a_host_in_development_carries_a_warning(pt1):
    for slug, expected in (("ecoli_bl21", False), ("hsapiens", True)):
        d = service.design(req(pt1, upstream="NdeI:144", downstream="XhoI:174", host=slug))
        codes = [i["code"] for i in d["issues"]]
        assert ("host_in_development" in codes) is expected, slug
        assert d["ok"], "a warning, not a refusal"


def test_the_page_grid_has_the_form_and_one_stage_so_the_result_stays_in_its_column():
    """A third direct child of <main> once pushed the result under the form."""
    from html.parser import HTMLParser

    class Tree(HTMLParser):
        def __init__(self):
            super().__init__()
            self.depth, self.in_main, self.children = 0, False, []

        def handle_starttag(self, tag, attrs):
            if tag == "main":
                self.in_main, self.depth = True, 0
            elif self.in_main:
                if self.depth == 0 and tag not in ("input", "br"):
                    self.children.append((tag, dict(attrs).get("id") or dict(attrs).get("class")))
                if tag not in ("input", "br", "meta", "link", "img", "hr"):
                    self.depth += 1

        def handle_endtag(self, tag):
            if tag == "main":
                self.in_main = False
            elif self.in_main and tag not in ("input", "br", "meta", "link", "img", "hr"):
                self.depth -= 1

    tree = Tree()
    tree.feed(static_text("index.html"))
    assert tree.children == [("form", "design-form"), ("div", "stage")], tree.children
    html = static_text("index.html")
    stage = html[html.index('class="stage"') :]
    assert stage.index('id="result"') < stage.index('id="batch-box"'), "the result comes first"


def test_the_fonts_are_served_from_the_package_and_licensed(server):
    css = static_text("app.css")
    for name in ("Geist-Variable.woff2", "GeistMono-Variable.woff2"):
        assert f"/fonts/{name}" in css
        resp, data = call(server, "GET", f"/fonts/{name}")
        assert resp.status == 200 and resp.getheader("Content-Type") == "font/woff2"
        assert data[:4] == b"wOF2" and 20_000 < len(data) < 500_000
    assert "Open Font License" in (
        resources.files("p2d_ui").joinpath("static", "fonts", "OFL.txt").read_text()
    )
    resp, _ = call(server, "GET", "/fonts/../app.py")
    assert resp.status == 404


def test_no_stylesheet_or_script_reaches_outside_this_computer():
    for name in ("app.css", "app.js", "index.html"):
        text = static_text(name)
        assert "https://" not in text.replace("https://github.com", "") or name == "index.html"
        assert "fonts.googleapis" not in text and "cdn." not in text


def test_two_stop_codons_can_be_asked_for_and_other_counts_are_refused(pt1):
    base = dict(upstream="NdeI:144", downstream="XhoI:174", cterm="stop")
    one = service.design(req(pt1, **base))
    two = service.design(req(pt1, **base, stop_codons=2))
    assert one["insert_bp"] + 3 == two["insert_bp"] and "TAATGA" in two["insert_dna"]
    assert "TAATGA" not in one["insert_dna"]
    for bad in (0, 3, "x", True):
        with pytest.raises(service.InputError, match="must be 1 or 2"):
            service.design(req(pt1, **base, stop_codons=bad))
    assert 'id="stop-codons"' in static_text("index.html")
