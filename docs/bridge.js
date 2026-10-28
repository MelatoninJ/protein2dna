/* p2d in the browser.
 *
 * app.js was written against the local server: it POSTs to /api/vector,
 * /api/sites, /api/pairs, /api/design, /api/batch and GETs /api/meta. Here
 * there is no server. Pyodide runs the same Python package inside the page,
 * and this file answers those same requests from it, so app.js runs unchanged.
 *
 * The error contract mirrors p2d_ui/server.py exactly:
 *   service.InputError -> 422 with the message
 *   anything else      -> 500 with a generic message
 */
(function () {
  "use strict";

  var PYODIDE = "https://cdn.jsdelivr.net/pyodide/v0.29.5/full/pyodide.mjs";
  var WHEEL = "protein2dna-0.3.0.dev1-py3-none-any.whl";

  var ROUTES = {
    "/api/vector": "load_vector",
    "/api/sites": "sites",
    "/api/pairs": "pairs",
    "/api/design": "design",
    "/api/batch": "batch"
  };

  /* ---------------------------------------------------------------- curtain */

  var style = document.createElement("style");
  style.textContent = [
    "#p2d-boot{position:fixed;inset:0;z-index:999;display:grid;place-items:center;",
    "  background:var(--bg,#0c0d10);font-family:var(--font-sans,system-ui,sans-serif);",
    "  transition:opacity .35s ease}",
    "#p2d-boot.gone{opacity:0;pointer-events:none}",
    ".p2d-boot-card{width:min(90vw,440px);color:var(--text,#ececf0)}",
    ".p2d-boot-h{margin:0 0 6px;font-size:19px;font-weight:600;letter-spacing:-.01em}",
    ".p2d-boot-s{margin:0 0 18px;font-size:14px;color:var(--muted,#9a9ca8)}",
    ".p2d-boot-bar{height:3px;border-radius:2px;overflow:hidden;",
    "  background:var(--sunken,#1b1c22)}",
    ".p2d-boot-bar i{display:block;height:100%;width:0;border-radius:2px;",
    "  background:var(--accent,#8da0ff);transition:width .5s ease}",
    ".p2d-boot-n{margin:18px 0 0;font-size:12.5px;line-height:1.6;color:var(--muted,#9a9ca8)}",
    ".p2d-boot-e{color:var(--error,#ff8a80)}"
  ].join("");
  document.head.appendChild(style);

  var boot = document.createElement("div");
  boot.id = "p2d-boot";
  boot.innerHTML =
    '<div class="p2d-boot-card">' +
    '<p class="p2d-boot-h">Starting p2d</p>' +
    '<p class="p2d-boot-s" id="p2d-boot-step">loading Python…</p>' +
    '<div class="p2d-boot-bar"><i id="p2d-boot-fill"></i></div>' +
    '<p class="p2d-boot-n">Python and Biopython are fetched once (about 12 MB, then cached) and ' +
    'run inside this tab. There is no server, so nothing you paste is sent anywhere.</p>' +
    "</div>";

  function mount() {
    if (document.body) document.body.appendChild(boot);
    else document.addEventListener("DOMContentLoaded", mount);
  }
  mount();

  function step(text, pct) {
    var s = document.getElementById("p2d-boot-step");
    var f = document.getElementById("p2d-boot-fill");
    if (s) s.textContent = text;
    if (f) f.style.width = pct + "%";
  }

  function fail(message) {
    var s = document.getElementById("p2d-boot-step");
    if (s) {
      s.className = "p2d-boot-s p2d-boot-e";
      s.textContent = message;
    }
  }

  /* ------------------------------------------------------------------ boot */

  // Runs in Python. Returns a JSON string so no proxy objects cross the border.
  // ?debug=1 adds the Python traceback to a failed answer. Off by default so
  // the page behaves exactly like the local app.
  var DEBUG = /[?&]debug=1\b/.test(location.search);

  var SHIM = [
    "import json, traceback",
    "import p2d_ui.service as _svc",
    "",
    "_p2d_debug = False",
    "",
    "def _p2d_call(name, payload):",
    "    fn = getattr(_svc, name)",
    "    try:",
    "        # A JS null arrives as a JsNull object rather than None, so this",
    "        # bridge sends the empty string to mean 'no body'.",
    "        out = fn() if not payload else fn(json.loads(payload))",
    "    except _svc.InputError as exc:",
    "        return json.dumps({'status': 422, 'body': {'error': str(exc)}})",
    "    except Exception:",
    "        body = {'error': 'The design failed unexpectedly.'}",
    "        if _p2d_debug:",
    "            body['detail'] = traceback.format_exc()",
    "        return json.dumps({'status': 500, 'body': body})",
    "    return json.dumps({'status': 200, 'body': out})",
    ""
  ].join("\n");

  var call = null;

  var ready = (async function () {
    try {
      step("loading Python…", 10);
      var mod = await import(PYODIDE);
      var py = await mod.loadPyodide();

      step("loading Biopython…", 40);
      await py.loadPackage(["micropip", "biopython", "pyyaml"]);

      step("loading codon tables…", 65);
      var micropip = py.pyimport("micropip");
      await micropip.install("python_codon_tables");

      step("loading p2d…", 85);
      // deps:false — biopython, pyyaml and the codon tables are already in, and
      // typer is only used by the command line, which does not exist here.
      await micropip.install(new URL(WHEEL, document.baseURI).href, false, false);

      step("ready", 100);
      py.runPython(SHIM);
      if (DEBUG) py.runPython("_p2d_debug = True");
      call = py.globals.get("_p2d_call");

      boot.classList.add("gone");
      setTimeout(function () {
        if (boot.parentNode) boot.parentNode.removeChild(boot);
      }, 400);
      return true;
    } catch (err) {
      fail("Could not start p2d: " + (err && err.message ? err.message : err));
      throw err;
    }
  })();

  /* ------------------------------------------------------------- the shim */

  var realFetch = window.fetch.bind(window);

  function reply(raw) {
    var out = JSON.parse(raw);
    return new Response(JSON.stringify(out.body), {
      status: out.status,
      headers: { "Content-Type": "application/json" }
    });
  }

  window.fetch = async function (input, init) {
    var url = typeof input === "string" ? input : (input && input.url) || "";
    var path = url.replace(/^[a-z]+:\/\/[^/]+/i, "").split("?")[0];

    if (path.endsWith("/api/meta")) {
      await ready;
      return reply(call("meta", ""));
    }

    var hit = Object.keys(ROUTES).filter(function (r) {
      return path.endsWith(r);
    })[0];

    if (hit) {
      await ready;
      var body = (init && init.body) || "{}";
      if (typeof body !== "string") body = "{}";
      return reply(call(ROUTES[hit], body));
    }

    return realFetch(input, init);
  };
})();
