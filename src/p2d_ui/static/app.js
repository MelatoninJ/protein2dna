"use strict";
/* p2d UI. Vanilla JS, no build step. Server text is only ever inserted with
   textContent, never as HTML. */

const GFP_EXAMPLE = "SKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLKFICTTGKLPVPWPTLVTTFSYGVQCFSRYPDHMKQHDFFKSAMPEGYVQERTIFFKDDGNYKTRAEVKFEGDTLVNRIELKGIDFKEDGNILGHKLEYNYNSHNVYIMADKQKNGIKVNFKIRHNIEDGSVQLADHYQQNTPIGDGPVLLPDNHYLSTQSALSKDPNEKRDHMVLLEFVTAAGITHGMDELYK"; // UniProt P42212 without its first Met (NdeI supplies it)
const LABELS = {
  "vector (N-terminal)": "N-terminal vector",
  "scar: 5' junction": "5' scar",
  "insert": "Your protein",
  "scar: 3' junction": "3' scar",
  "vector (C-terminal)": "C-terminal vector",
};
const SCAR_NOTE = "Residues added by the restriction sites.";
const RESIDUES = /^[ACDEFGHIKLMNPQRSTVWY]*$/;

const $ = (id) => document.getElementById(id);
let meta = null;
let busy = false;

function h(tag, attrs, ...kids) {
  const node = document.createElement(tag);
  for (const [k, v] of Object.entries(attrs || {})) {
    if (v === false || v == null) continue;
    if (k === "class") node.className = v;
    else if (k === "text") node.textContent = v;
    else if (k.startsWith("on")) node.addEventListener(k.slice(2), v);
    else node.setAttribute(k, v === true ? "" : v);
  }
  for (const kid of kids.flat()) {
    if (kid == null || kid === false) continue;
    node.append(kid.nodeType ? kid : document.createTextNode(String(kid)));
  }
  return node;
}

/* ------------------------------------------------------------------ theme */
function effectiveTheme() {
  const set = document.documentElement.dataset.theme;
  if (set) return set;
  return matchMedia("(prefers-color-scheme: dark)").matches ? "dark" : "light";
}
function initTheme() {
  try {
    const saved = localStorage.getItem("p2d-theme");
    if (saved === "light" || saved === "dark") document.documentElement.dataset.theme = saved;
  } catch (_) { /* storage can be unavailable; the system theme still works */ }
  $("theme-toggle").addEventListener("click", () => {
    const next = effectiveTheme() === "dark" ? "light" : "dark";
    document.documentElement.dataset.theme = next;
    try { localStorage.setItem("p2d-theme", next); } catch (_) {}
  });
}

/* ------------------------------------------------------------------ input */
function cleanProtein(text) {
  const lines = text.replace(/\r/g, "").split("\n");
  const body = [];
  let header = false;
  for (const line of lines) {
    if (line.startsWith(">")) {
      if (header || body.length) break;
      header = true;
      continue;
    }
    body.push(line);
  }
  return body.join("").replace(/[\s\d]+/g, "").toUpperCase().replace(/\*/g, "");
}

function checkProtein(showEmpty) {
  const seq = cleanProtein($("protein").value);
  const bad = [...new Set([...seq].filter((c) => !RESIDUES.test(c)))].sort();
  $("protein-count").textContent = `${seq.length} aa`;
  let msg = "";
  if (bad.length) msg = `Unsupported character(s): ${bad.join(", ")}. Use the 20 standard one-letter codes.`;
  else if (seq.length > meta.max_residues) msg = `That is ${seq.length} residues; the limit is ${meta.max_residues}.`;
  else if (showEmpty && !seq.length) msg = "Paste a protein sequence to design.";
  const err = $("protein-error");
  err.hidden = !msg;
  err.textContent = msg;
  if (msg) $("protein").setAttribute("aria-invalid", "true");
  else $("protein").removeAttribute("aria-invalid");
  return !msg;
}

/* --------------------------------------------------------------- form init */
async function post(path, body) {
  const res = await fetch(path, { method: "POST", headers: { "Content-Type": "application/json" }, body: JSON.stringify(body) });
  const data = await res.json();
  if (!res.ok) throw new Error(data.error || "Request failed.");
  return data;
}

function showHostNote() {
  const host = meta.hosts.find((x) => x.slug === $("host").value);
  $("host-note").textContent = host
    ? host.verified ? "Codon table from a cited source." : "Codon table is an unverified placeholder."
    : "";
}

/* ---------------------------------------------------------- enzyme picker */
let enz = null;      // last /api/enzymes answer
let enzSeq = 0;      // drops answers that arrive out of order
let enzTimer = null;

const pos1 = (e) => e.start + 1; // show 1-based positions like a plasmid map

function optionFor(e) {
  const alias = e.aliases.length ? ` Also sold as ${e.aliases.join(", ")}.` : "";
  return h("option", { value: e.name, title: `${e.site}, ${e.overhang_kind} ${e.overhang} end.${alias}`, text: `${e.name} at ${pos1(e)}` });
}

function fillSelect(select, rows, verdict, keep, prefer) {
  const usable = rows.filter((e) => verdict(e).ok);
  const wanted = usable.find((e) => e.name === keep) || usable.find((e) => e.name === prefer) || usable[0];
  select.replaceChildren(...usable.map(optionFor));
  if (!usable.length) select.append(h("option", { value: "", text: "Nothing usable" }));
  if (wanted) select.value = wanted.name;
  select.disabled = !usable.length;
  return wanted ? wanted.name : "";
}

async function refreshEnzymes() {
  const my = ++enzSeq;
  const body = () => ({ vector: $("vector").value, host: $("host").value, protein: $("protein").value, upstream: $("upstream").value });
  try {
    let data = await post("/api/enzymes", body());
    if (my !== enzSeq) return;
    const upBefore = $("upstream").value;
    const up = fillSelect($("upstream"), data.enzymes, (e) => e.as_upstream, upBefore, "NdeI");
    if (up !== upBefore) {            // the 5' choice changed, so the 3' verdicts must be recomputed
      data = await post("/api/enzymes", { ...body(), upstream: up });
      if (my !== enzSeq) return;
    }
    enz = data;
    const downRows = data.enzymes.filter((e) => e.as_downstream && !e.as_downstream.hidden);
    fillSelect($("downstream"), downRows, (e) => e.as_downstream, $("downstream").value, "XhoI");
    renderMcs();
    renderWhy();
    renderPairNote();
  } catch (err) {
    if (my !== enzSeq) return;
    $("upstream").replaceChildren();
    $("downstream").replaceChildren();
    $("mcs").replaceChildren();
    $("pair-note").textContent = err.message;
    $("why").hidden = true;
  }
}

function scheduleEnzymes() {
  clearTimeout(enzTimer);
  enzTimer = setTimeout(refreshEnzymes, 300);
}

function renderMcs() {
  const box = $("mcs");
  const { start, stop_end } = enz.window;
  const span = stop_end - start;
  const pct = (p) => `${Math.min(100, Math.max(0, ((p - start) / span) * 100))}%`;
  const up = $("upstream").value, down = $("downstream").value;
  // enzymes that share a position share one tick and one label
  const byPos = new Map();
  for (const e of enz.enzymes) byPos.set(e.start, [...(byPos.get(e.start) || []), e]);
  const nodes = [h("div", { class: "mcs-bar" }), h("span", { class: "mcs-end l", text: "start" }), h("span", { class: "mcs-end r", text: "stop" })];
  [...byPos.entries()].sort((a, b) => a[0] - b[0]).forEach(([p, group], i) => {
    const role = group.some((e) => e.name === up) ? "up" : group.some((e) => e.name === down) ? "down" : "";
    const off = !role && group.every((e) => !e.as_upstream.ok && !(e.as_downstream && e.as_downstream.ok));
    const tick = h("i", { class: `mcs-tick ${role} ${off ? "off" : ""}` });
    const label = h("span", { class: `mcs-label r${i % 3} ${role} ${off ? "off" : ""}`, text: group.map((e) => e.name).join("/") });
    tick.style.left = pct(p);
    label.style.left = pct(p);
    nodes.push(tick, label);
  });
  box.replaceChildren(...nodes);
  const named = (n) => (n ? n : "none chosen");
  box.setAttribute("aria-label", `Enzyme sites along the vector expression region. 5' enzyme: ${named(up)}. 3' enzyme: ${named(down)}.`);
}

function renderWhy() {
  const items = [];
  for (const e of enz.enzymes) {
    if (!e.as_upstream.ok) items.push([e.name, `as the 5' enzyme: ${e.as_upstream.reason}`]);
    const d = e.as_downstream;
    if (d && !d.ok && !d.hidden && d.reason) items.push([e.name, `as the 3' partner of ${$("upstream").value}: ${d.reason}`]);
  }
  $("why-list").replaceChildren(...items.map(([n, r]) => h("li", null, h("b", { text: n }), " ", r)));
  $("why").hidden = !items.length;
}

function renderPairNote() {
  const up = enz.enzymes.find((e) => e.name === $("upstream").value);
  const down = $("downstream").value;
  if (!up) { $("pair-note").textContent = "No enzyme pair can clone this protein into this vector."; return; }
  const start = up.as_upstream.starts_at === "site"
    ? `${up.name} holds the vector's start codon. The vector's N-terminal tag is replaced and the ATG in the site becomes your Met, so leave the first M out of your sequence.`
    : `Translation starts at the vector's own start codon, so the vector sequence between it and ${up.name} stays on as an N-terminal tag.` +
      (up.site.endsWith("ATG") ? ` The ATG in ${up.name} becomes the first Met of your protein, so leave your own first M out.` : "");
  const downRow = enz.enzymes.find((e) => e.name === down);
  const cut = (e) => `${e.name} cuts ${e.site} and leaves a ${e.overhang_kind} ${e.overhang} end`;
  const ends = downRow ? ` ${cut(up)}; ${cut(downRow)}. The insert goes in one way round only.` : "";
  $("pair-note").textContent = `${start}${ends}`;
}

function initForm() {
  $("vector").replaceChildren(
    ...meta.vectors.map((v) => h("option", { value: v.name, text: `${v.name} (${v.length} bp)` })));
  $("host").replaceChildren(...meta.hosts.map((x) => h("option", { value: x.slug, text: x.name })));
  $("strategies").replaceChildren(
    ...meta.strategies.map((s, i) =>
      h("label", { class: "choice" },
        h("input", { type: "radio", name: "strategy", value: s.value, checked: i === 0 }),
        h("b", { text: s.label }),
        h("span", { text: s.note }))));
  showHostNote();

  $("vector").addEventListener("change", () => { $("upstream").replaceChildren(); $("downstream").replaceChildren(); refreshEnzymes(); });
  $("host").addEventListener("change", () => { showHostNote(); refreshEnzymes(); });
  $("upstream").addEventListener("change", refreshEnzymes);
  $("downstream").addEventListener("change", () => { if (enz) { renderMcs(); renderPairNote(); } });
  $("protein").addEventListener("input", () => { checkProtein(false); scheduleEnzymes(); });
  $("protein").addEventListener("keydown", (e) => {
    if (e.key === "Enter" && (e.metaKey || e.ctrlKey)) { e.preventDefault(); $("design-form").requestSubmit(); }
  });
  $("load-example").addEventListener("click", loadExample);
  $("design-form").addEventListener("submit", onSubmit);
  refreshEnzymes();
}

function loadExample() {
  $("protein").value = GFP_EXAMPLE;
  checkProtein(false);
  refreshEnzymes();
  $("protein").focus();
}

/* ---------------------------------------------------------------- results */
function setResult(...nodes) {
  const box = $("result");
  box.classList.remove("fresh");
  box.replaceChildren(...nodes);
  void box.offsetWidth; // restart the entry animation
  box.classList.add("fresh");
}

function renderEmpty() {
  setResult(h("div", { class: "empty" },
    h("h2", { text: "Paste a protein, get DNA you can clone" }),
    h("ol", null,
      h("li", null, h("b", { text: "Paste a sequence." }), " Plain amino acids or FASTA."),
      h("li", null, h("b", { text: "Choose the vector and enzymes." }), " The reading frame is worked out for both junctions."),
      h("li", null, h("b", { text: "Read the fusion map." }), " It shows the exact protein your construct will make, scar residues included."),
      h("li", null, h("b", { text: "Take the DNA." }), " Copy the insert, or download FASTA or a full GenBank plasmid.")),
    h("button", { type: "button", class: "secondary try", onclick: loadExample, text: "Load GFP example" })));
}

function renderLoading() {
  setResult(h("div", { class: "skel", "aria-busy": "true" },
    h("div", { class: "skel-head" }, h("i", { class: "s1" }), h("i", { class: "s2" })),
    h("i", { class: "s3" }),
    h("i", { class: "s4" })));
}

function renderError(message) {
  setResult(h("div", { class: "error-box", role: "alert" },
    h("h2", { text: "Could not design that" }),
    h("p", { text: message })));
}

function stat(label, value) {
  return [h("dt", { text: label }), h("dd", { text: value })];
}

function mapSegment(seg, total) {
  // short segments show their residues; long vector runs get a plain label that always fits
  const text = seg.kind === "payload" ? `${seg.aa.length} aa` : seg.aa.length <= 6 ? seg.aa : "vector";
  const node = h("div", { class: `seg ${seg.kind === "payload" ? "payload" : seg.kind === "scar" ? "scar" : "vector"}`,
    title: `${LABELS[seg.label] || seg.label}: residues ${seg.start} to ${seg.end}`, text });
  node.style.setProperty("--w", String(Math.max(seg.aa.length, total * 0.04)));
  return node;
}

function legendRow(seg) {
  const long = seg.aa.length > 90;
  const aa = h("p", { class: `aa${long ? " clip" : ""}`, text: seg.aa });
  const toggle = long
    ? h("button", { type: "button", class: "link", text: "Show all",
        onclick: (e) => { const open = aa.classList.toggle("clip"); e.target.textContent = open ? "Show all" : "Show less"; } })
    : null;
  return h("div", { class: "legend-row" },
    h("div", null,
      h("span", { class: "what", text: LABELS[seg.label] || seg.label }),
      h("span", { class: "range", text: seg.start === seg.end ? `residue ${seg.start}` : `residues ${seg.start} to ${seg.end}` })),
    h("div", null, aa, seg.kind === "scar" ? h("p", { class: "help", text: SCAR_NOTE }) : null, toggle));
}

function dnaView(r) {
  const [a, b] = r.coding_span;
  const seq = r.insert_dna;
  const rows = [];
  for (let i = 0; i < seq.length; i += 60) {
    const end = Math.min(i + 60, seq.length);
    const parts = [];
    const cut = (from, to, cls) => { if (to > from) parts.push(h("span", { class: cls, text: seq.slice(from, to) })); };
    cut(i, Math.min(end, a), "flank");
    cut(Math.max(i, a), Math.min(end, b), "cds");
    cut(Math.max(i, b), end, "flank");
    rows.push(h("div", { class: "dna-line" }, h("span", { class: "n", text: String(i + 1) }), h("span", null, parts)));
  }
  return h("div", { class: "dna", tabindex: "0", role: "group", "aria-label": "Insert DNA sequence" }, rows);
}

function download(name, text) {
  const url = URL.createObjectURL(new Blob([text], { type: "text/plain" }));
  const a = h("a", { href: url, download: name });
  document.body.append(a);
  a.click();
  a.remove();
  setTimeout(() => URL.revokeObjectURL(url), 1000);
}

async function copyText(btn, text) {
  const original = btn.textContent;
  try {
    await navigator.clipboard.writeText(text);
    btn.textContent = "Copied";
  } catch (_) {
    btn.textContent = "Copy failed";
  }
  setTimeout(() => { btn.textContent = original; }, 1600);
}

function renderResult(r) {
  const total = r.segments.reduce((n, s) => n + s.aa.length, 0);
  const verdict = h("div", { class: `verdict ${r.ok ? "pass" : "fail"}` },
    h("h2", { text: r.ok ? "Passed the round-trip check" : "Failed validation" }),
    h("p", { class: "where", text: `${r.vector}, ${r.site.replace("-", " to ")}, ${r.payload_length} residues` }));

  const stats = h("dl", { class: "stats" },
    stat("Insert", `${r.insert_bp} bp`),
    stat("Plasmid", `${r.plasmid_bp} bp`),
    stat("CAI", r.metrics.cai == null ? "n/a" : r.metrics.cai.toFixed(2)),
    stat("GC", r.metrics.gc == null ? "n/a" : `${(r.metrics.gc * 100).toFixed(1)}%`),
    stat("Rare codons", String(r.metrics.n_rare ?? 0)));

  const map = h("section", { "aria-labelledby": "map-h" },
    h("div", { class: "map-title" }, h("h2", { id: "map-h", text: "Fusion protein" }), h("p", { text: `${total} residues from the start codon` })),
    h("div", { class: "map", role: "img", "aria-label": "Fusion protein map" }, r.segments.map((s) => mapSegment(s, total))),
    h("div", { class: "legend" }, r.segments.map(legendRow)));

  const issues = h("section", { "aria-labelledby": "issues-h" },
    h("h2", { id: "issues-h", text: "Checks" }),
    r.issues.length
      ? h("ul", { class: "issues" }, r.issues.map((i) =>
          h("li", { class: `issue ${i.severity}` },
            h("div", null, h("span", { class: "sev", text: i.severity[0].toUpperCase() + i.severity.slice(1) })),
            h("div", null, h("p", { text: i.message }), h("span", { class: "code", text: i.code })))))
      : h("p", { class: "none", text: "No warnings. The insert is in frame at both junctions and carries no extra cutting sites." }));

  const copyBtn = h("button", { type: "button", class: "secondary", text: "Copy DNA" });
  copyBtn.addEventListener("click", () => copyText(copyBtn, r.insert_dna));
  const dna = h("section", { "aria-labelledby": "dna-h" },
    h("div", { class: "dna-head" },
      h("h2", { id: "dna-h", text: "Insert DNA" }),
      h("div", { class: "actions" }, copyBtn,
        h("button", { type: "button", class: "secondary", text: "Download FASTA", onclick: () => download(`${r.filename}_insert.fasta`, r.fasta) }),
        h("button", { type: "button", class: "secondary", text: "Download GenBank", onclick: () => download(`${r.filename}.gb`, r.genbank) }))),
    dnaView(r),
    h("p", { class: "key", text: "Dark bases are your coding region. Lighter bases are the cloning sites and linker bases around it." }));

  setResult(h("div", null, verdict, stats), map, issues, dna);
}

/* ----------------------------------------------------------------- submit */
async function onSubmit(e) {
  e.preventDefault();
  if (busy || !checkProtein(true)) { if (!busy) $("protein").focus(); return; }
  if (!$("upstream").value || !$("downstream").value) { renderError("Choose a usable 5' and 3' enzyme first."); return; }
  busy = true;
  const btn = $("submit");
  btn.disabled = true;
  btn.textContent = "Designing...";
  renderLoading();
  const form = new FormData($("design-form"));
  const body = {
    protein: form.get("protein"), vector: form.get("vector"),
    upstream: form.get("upstream"), downstream: form.get("downstream"),
    host: form.get("host"), strategy: form.get("strategy"), seed: form.get("seed"),
  };
  try {
    const res = await fetch("/api/design", { method: "POST", headers: { "Content-Type": "application/json" }, body: JSON.stringify(body) });
    const data = await res.json();
    if (!res.ok) renderError(data.error || "The design failed.");
    else { renderResult(data); $("result").focus({ preventScroll: false }); }
  } catch (_) {
    renderError("Could not reach the p2d server. Check that `p2d ui` is still running in your terminal.");
  } finally {
    busy = false;
    btn.disabled = false;
    btn.textContent = "Design insert";
  }
}

/* ------------------------------------------------------------------- boot */
async function boot() {
  initTheme();
  try {
    meta = await (await fetch("/api/meta")).json();
  } catch (_) {
    renderError("Could not load the p2d server data. Check that `p2d ui` is running.");
    return;
  }
  initForm();
  renderEmpty();
}
document.addEventListener("DOMContentLoaded", boot);
