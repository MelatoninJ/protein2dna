"""Local web server for the p2d UI.

Standard library only, so ``pip install p2d`` is enough to run ``p2d ui``.  It
is a single-user tool bound to loopback: it serves three fixed static files and
a two-endpoint JSON API, and refuses requests whose Host header is not the
loopback address it was bound to (which also blocks DNS-rebinding attacks).
"""

from __future__ import annotations

import json
import threading
import webbrowser
from http import HTTPStatus
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from importlib import resources

from . import service

MAX_BODY = 1_000_000  # bytes; 5000 residues is ~5 KB, so this is generous
STATIC = {
    "/": ("index.html", "text/html; charset=utf-8"),
    "/app.css": ("app.css", "text/css; charset=utf-8"),
    "/app.js": ("app.js", "text/javascript; charset=utf-8"),
}
SECURITY_HEADERS = {
    "Cache-Control": "no-store",
    "X-Content-Type-Options": "nosniff",
    "Referrer-Policy": "no-referrer",
    "Content-Security-Policy": (
        "default-src 'self'; img-src 'self' data:; style-src 'self'; "
        "script-src 'self'; connect-src 'self'; base-uri 'none'; form-action 'none'"
    ),
}


def _static_bytes(name: str) -> bytes:
    return resources.files("p2d_ui").joinpath("static", name).read_bytes()


POST_ROUTES = {
    "/api/vector": service.load_vector,
    "/api/sites": service.sites,
    "/api/design": service.design,
}
# an uploaded plasmid arrives base64-encoded, so its route may carry up to ~5 MB of file
ROUTE_LIMITS = {"/api/vector": 7_000_000}


class Handler(BaseHTTPRequestHandler):
    server_version = "p2d-ui"

    # ------------------------------------------------------------ plumbing
    def log_message(self, format: str, *args) -> None:  # noqa: A002 - stdlib signature
        pass

    def _send(self, status: int, body: bytes, content_type: str) -> None:
        self.send_response(status)
        self.send_header("Content-Type", content_type)
        self.send_header("Content-Length", str(len(body)))
        for k, v in SECURITY_HEADERS.items():
            self.send_header(k, v)
        self.end_headers()
        self.wfile.write(body)

    def _json(self, status: int, payload: dict) -> None:
        self._send(status, json.dumps(payload).encode(), "application/json; charset=utf-8")

    def _host_ok(self) -> bool:
        allowed = {f"127.0.0.1:{self.server.server_port}", f"localhost:{self.server.server_port}"}
        return self.headers.get("Host", "") in allowed

    # -------------------------------------------------------------- routes
    def do_GET(self) -> None:
        if not self._host_ok():
            return self._json(HTTPStatus.FORBIDDEN, {"error": "Unexpected Host header."})
        path = self.path.split("?", 1)[0]
        if path in STATIC:
            name, ctype = STATIC[path]
            return self._send(HTTPStatus.OK, _static_bytes(name), ctype)
        if path == "/api/meta":
            return self._json(HTTPStatus.OK, service.meta())
        self._json(HTTPStatus.NOT_FOUND, {"error": "Not found."})

    def do_POST(self) -> None:
        if not self._host_ok():
            return self._json(HTTPStatus.FORBIDDEN, {"error": "Unexpected Host header."})
        route = self.path.split("?", 1)[0]
        handler = POST_ROUTES.get(route)
        if handler is None:
            return self._json(HTTPStatus.NOT_FOUND, {"error": "Not found."})
        if not self.headers.get("Content-Type", "").startswith("application/json"):
            return self._json(HTTPStatus.UNSUPPORTED_MEDIA_TYPE, {"error": "Send JSON."})
        try:
            length = int(self.headers.get("Content-Length", "0"))
        except ValueError:
            length = -1
        if not 0 < length <= ROUTE_LIMITS.get(route, MAX_BODY):
            return self._json(HTTPStatus.REQUEST_ENTITY_TOO_LARGE, {"error": "Bad request size."})
        try:
            request = json.loads(self.rfile.read(length))
            if not isinstance(request, dict):
                raise ValueError
        except ValueError:
            return self._json(HTTPStatus.BAD_REQUEST, {"error": "Body must be a JSON object."})
        try:
            result = handler(request)
        except service.InputError as exc:
            return self._json(HTTPStatus.UNPROCESSABLE_ENTITY, {"error": str(exc)})
        except Exception:  # noqa: BLE001 - never leak internals to the browser
            return self._json(
                HTTPStatus.INTERNAL_SERVER_ERROR, {"error": "The design failed unexpectedly."}
            )
        self._json(HTTPStatus.OK, result)


def make_server(host: str = "127.0.0.1", port: int = 8765) -> ThreadingHTTPServer:
    return ThreadingHTTPServer((host, port), Handler)


def serve(port: int = 8765, open_browser: bool = True) -> None:
    server = make_server(port=port)
    url = f"http://127.0.0.1:{server.server_port}/"
    print(f"p2d UI running at {url}  (Ctrl+C to stop)")
    if open_browser:
        threading.Timer(0.4, webbrowser.open, args=(url,)).start()
    try:
        server.serve_forever()
    except KeyboardInterrupt:
        pass
    finally:
        server.server_close()
