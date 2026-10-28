#!/usr/bin/env bash
# Rebuild docs/ from src/ after changing the p2d UI.
#
# docs/ is a snapshot of the package: app.js, app.css, index.html and a wheel.
# Editing src/p2d_ui/static/ alone does nothing to the published page — this
# script copies the new files across and rebuilds the wheel.
#
# Run it from the root of the protein2dna repository.

set -euo pipefail

SRC="src/p2d_ui/static"
DOCS="docs"

[ -d "$SRC" ] || { echo "run this from the protein2dna repo root"; exit 1; }
[ -d "$DOCS" ] || { echo "no docs/ folder — create it first"; exit 1; }

echo "==> wheel"
rm -f "$DOCS"/*.whl
python -m build --wheel -o "$DOCS" >/dev/null
WHEEL=$(basename "$DOCS"/*.whl)
echo "    $WHEEL"

echo "==> static files"
cp "$SRC/app.js" "$DOCS/app.js"
cp "$SRC/app.css" "$DOCS/app.css"
mkdir -p "$DOCS/fonts" && cp "$SRC"/fonts/* "$DOCS/fonts/"

python - "$SRC/index.html" "$DOCS/index.html" "$WHEEL" <<'PY'
import pathlib, re, sys

src, dst, wheel = (pathlib.Path(sys.argv[1]), pathlib.Path(sys.argv[2]), sys.argv[3])
s = src.read_text(encoding="utf-8")

# 1. absolute paths -> relative, because the page lives in a subfolder
s = s.replace('href="/app.css"', 'href="app.css"')
s = s.replace('src="/app.js"', 'src="app.js"')
s = s.replace('<a class="brand" href="/"', '<a class="brand" href="./"')

# 2. bridge.js must run before app.js
assert '<script src="app.js" defer></script>' in s
s = s.replace('<script src="app.js" defer></script>',
              '<!-- bridge.js must run first: it answers the /api/ calls app.js makes. -->\n'
              '  <script src="bridge.js"></script>\n'
              '  <script src="app.js" defer></script>')

# 3. bump the cache-busting number so browsers pick the new files up
old = dst.read_text(encoding="utf-8") if dst.exists() else ""
ver = max((int(n) for n in re.findall(r'\?v=(\d+)', old)), default=1) + 1
for name in ("app.css", "bridge.js", "app.js"):
    s = s.replace(f'"{name}"', f'"{name}?v={ver}"')

# 4. the local-only promise is true here for a different reason
s = s.replace("<title>p2d | Protein to DNA</title>",
              "<title>p2d | Protein to DNA, in your browser</title>")
old_line = "<p>Runs entirely on this computer. Nothing you paste or upload leaves it.</p>"
assert old_line in s, "the footer sentence moved — update sync-docs.sh"
s = s.replace(old_line,
    "<p>Runs entirely in this tab. Python itself is downloaded once and executed by your "
    "browser, so nothing you paste or upload is sent anywhere — there is no server to "
    "send it to.</p>")

dst.write_text(s, encoding="utf-8")
print(f"    index.html rewritten, cache version v={ver}")
PY

echo "==> wheel name in bridge.js"
python - "$DOCS/bridge.js" "$WHEEL" <<'PY'
import pathlib, re, sys
p = pathlib.Path(sys.argv[1]); s = p.read_text(encoding="utf-8")
s = re.sub(r'var WHEEL = "[^"]+";', f'var WHEEL = "{sys.argv[2]}";', s)
p.write_text(s, encoding="utf-8")
print("    pointed at", sys.argv[2])
PY

echo
echo "done. test it before pushing:"
echo "    cd $DOCS && python3 -m http.server 8080"
