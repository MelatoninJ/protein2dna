"use strict";
/* p2d UI. Vanilla JS, no build step. Server text is only ever inserted with
   textContent, never as HTML. */

const GFP_EXAMPLE = "SKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLKFICTTGKLPVPWPTLVTTFSYGVQCFSRYPDHMKQHDFFKSAMPEGYVQERTIFFKDDGNYKTRAEVKFEGDTLVNRIELKGIDFKEDGNILGHKLEYNYNSHNVYIMADKQKNGIKVNFKIRHNIEDGSVQLADHYQQNTPIGDGPVLLPDNHYLSTQSALSKDPNEKRDHMVLLEFVTAAGITHGMDELYK"; // UniProt P42212 without its first Met (NdeI and NcoI sites supply one)
const LABELS = {
  "vector (N-terminal)": "N-terminal vector",
  "scar: 5' junction": "5' scar",
  "insert": "Your protein",
  "scar: 3' junction": "3' scar",
  "vector (C-terminal)": "C-terminal vector",
};
const SCAR_NOTE = "Residues added by the restriction sites.";
const RESIDUES = /^[ACDEFGHIKLMNPQRSTVWY]*$/;
const MAX_UPLOAD = 5 * 1024 * 1024;
const MAX_LABEL_ROWS = 6;
const ROLE_TEXT = { "backbone (kept)": "backbone, kept", "replaced piece": "replaced piece", "insert": "insert", "fragment": "fragment" };

const $ = (id) => document.getElementById(id);
const bp = (n) => `${n.toLocaleString("en-US")} bp`;
const pos1 = (n) => n + 1; // plasmid maps number from 1

let meta = null;
let vec = null;          // the loaded plasmid (answer of /api/vector)
let enz = null;          // the last /api/sites answer
let winTouched = false;  // true once the user typed their own window
let busy = false;
let siteSeq = 0;         // drops answers that arrive out of order
let siteTimer = null;

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

async function post(path, body) {
  const res = await fetch(path, { method: "POST", headers: { "Content-Type": "application/json" }, body: JSON.stringify(body) });
  let data = {};
  try { data = await res.json(); } catch (_) { /* a non-JSON failure is reported below */ }
  if (!res.ok) throw new Error(data.error || "Request failed.");
  return data;
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

/* ---------------------------------------------------------------- protein */
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

/* ------------------------------------------------------------ plasmid input */
const loaded = new Map(); // id -> { label, vec } for plasmids loaded this session

function vectorOptions(selected) {
  const groups = [];
  groups.push(h("optgroup", { label: "Bundled" }, meta.bundled.map((b) => h("option", { value: `bundled:${b.id}`, text: `${b.name} (${bp(b.length)})` }))));
  if (meta.library.length) groups.push(h("optgroup", { label: "Your vector folder" }, meta.library.map((f) => h("option", { value: `library:${f}`, text: f }))));
  if (loaded.size) groups.push(h("optgroup", { label: "Loaded this session" }, [...loaded].map(([id, e]) => h("option", { value: `loaded:${id}`, text: e.label }))));
  $("vector").replaceChildren(...groups);
  if (selected) $("vector").value = selected;
}

function setStatus(text, isError) {
  const el = $("vector-status");
  el.textContent = text || "";
  el.classList.toggle("field-error", !!isError);
}

async function loadVector(request, label) {
  setStatus("Reading the plasmid...");
  try {
    const v = await post("/api/vector", request);
    loaded.set(v.id, { label: `${v.name} (${bp(v.length)})`, vec: v });
    vectorOptions(`loaded:${v.id}`);
    setStatus("");
    showVector(v);
  } catch (err) {
    setStatus(err.message, true);
  }
}

function onVectorSelect() {
  const [kind, ...rest] = $("vector").value.split(":");
  const value = rest.join(":");
  if (kind === "loaded") return showVector(loaded.get(value).vec);
  return loadVector(kind === "bundled" ? { bundled: value } : { library: value });
}

function toBase64(buffer) {
  const bytes = new Uint8Array(buffer);
  let binary = "";
  for (let i = 0; i < bytes.length; i += 0x8000) binary += String.fromCharCode(...bytes.subarray(i, i + 0x8000));
  return btoa(binary);
}

async function onFile(e) {
  const file = e.target.files[0];
  e.target.value = "";
  if (!file) return;
  if (file.size > MAX_UPLOAD) return setStatus(`${file.name} is larger than 5 MB.`, true);
  const data = toBase64(await file.arrayBuffer());
  loadVector({ filename: file.name, data_b64: data });
}

function showVector(v) {
  vec = v;
  winTouched = false;
  enz = null;
  $("vector-card").hidden = false;

  const topology = v.circular ? "circular" : "linear";
  const extra = [`Source: ${v.source}`, v.flipped ? "Stored on the opposite strand, flipped to read forward" : null, ...v.notes.filter((n) => !/Reverse-complemented/.test(n)).map((n) => n.replace(/\.$/, ""))].filter(Boolean);
  $("vector-summary").replaceChildren(`${v.name}, ${bp(v.length)}, ${topology}`, h("small", { text: `${v.site_count} restriction sites, ${v.single_cutters} of them cut once. ${extra.join(". ")}.` }));
  $("vector-warnings").replaceChildren(...v.warnings.map((w) => h("li", { text: w })));

  const startSelect = $("start");
  startSelect.replaceChildren(
    ...v.starts.map((s, i) => h("option", { value: String(s.position), text: `ATG at ${pos1(s.position)} (${s.confidence})${i === 0 ? "" : ""}` })),
    h("option", { value: "manual", text: "Other position..." }));
  startSelect.value = v.starts.length ? String(v.starts[0].position) : "manual";
  onStartChange(false);

  const feats = v.features;
  $("features-summary").textContent = `Annotated features (${feats.length})`;
  $("features-list").replaceChildren(...feats.map((f) => h("li", null, h("b", { text: f.label }), " ", h("span", { class: "pos", text: `${f.type}, ${pos1(f.start)} to ${f.end}, ${f.strand === -1 ? "minus" : "plus"} strand` }))));
  $("features-box").hidden = !feats.length;

  $("win-mcs").hidden = !v.mcs;
  $("sites-field").disabled = false;
  refreshSites();
}

function currentStart() {
  if ($("start").value === "manual") {
    const n = parseInt($("start-manual").value, 10);
    return Number.isFinite(n) && n >= 1 ? n - 1 : null;
  }
  const n = parseInt($("start").value, 10);
  return Number.isFinite(n) ? n : null;
}

function onStartChange(refresh = true) {
  const manual = $("start").value === "manual";
  $("start-manual-box").hidden = !manual;
  const cand = vec.starts.find((s) => String(s.position) === $("start").value);
  const note = $("start-note");
  if (manual) note.textContent = "You are choosing the ATG yourself. It must be the start codon the host really uses.";
  else if (cand) {
    const nterm = cand.n_terminal ? ` The vector encodes ${cand.n_terminal.length > 34 ? cand.n_terminal.slice(0, 34) + "..." : cand.n_terminal} before your site.` : "";
    note.textContent = `${cand.evidence.join("; ")}.${nterm}`;
  }
  if (refresh) { winTouched = false; refreshSites(); }
}

/* ----------------------------------------------------------- restriction sites */
function siteLabel(s) {
  const which = s.count > 1 ? ` (${s.number} of ${s.count})` : "";
  return `${s.enzyme} at ${pos1(s.start)}${which}`;
}

function optionFor(s) {
  const alias = s.aliases.length ? ` Also sold as ${s.aliases.join(", ")}.` : "";
  return h("option", { value: s.id, title: `${s.site}, leaves a ${s.overhang_kind} ${s.overhang} end.${alias}`, text: siteLabel(s) });
}

function fillSelect(select, rows, verdict, keep, prefer) {
  const usable = rows.filter((s) => verdict(s).ok);
  const wanted = usable.find((s) => s.id === keep) || prefer(usable);
  select.replaceChildren(...usable.map(optionFor));
  if (!usable.length) select.append(h("option", { value: "", text: "Nothing usable" }));
  if (wanted) select.value = wanted.id;
  select.disabled = !usable.length;
  return wanted ? wanted.id : "";
}

function windowParam() {
  if (!winTouched) return undefined;
  const a = parseInt($("win-from").value, 10), b = parseInt($("win-to").value, 10);
  return Number.isFinite(a) && Number.isFinite(b) ? [a - 1, b] : undefined;
}

function siteBody(extra) {
  return {
    vector_id: vec.id, start: currentStart(), protein: $("protein").value, host: $("host").value,
    window: windowParam(), include_multi: $("multi").checked, cterm: document.querySelector('input[name="cterm"]:checked').value,
    ...extra,
  };
}

async function refreshSites() {
  if (!vec) return;
  const my = ++siteSeq;
  try {
    let data = await post("/api/sites", siteBody({ upstream: $("upstream").value }));
    if (my !== siteSeq) return;
    const before = $("upstream").value;
    const up = fillSelect($("upstream"), data.sites, (s) => s.as_upstream, before, (u) => u.find((s) => s.enzyme === "NdeI") || u.find((s) => s.enzyme === "NcoI") || u[0]);
    if (up !== before) {
      data = await post("/api/sites", siteBody({ upstream: up }));
      if (my !== siteSeq) return;
    }
    enz = data;
    if (!winTouched) { $("win-from").value = pos1(data.window.start); $("win-to").value = data.window.end; }
    const downRows = data.sites.filter((s) => s.as_downstream && !s.as_downstream.hidden);
    fillSelect($("downstream"), downRows, (s) => s.as_downstream, $("downstream").value, (u) => u.find((s) => s.enzyme === "XhoI") || u[u.length - 1]);
    syncCterm();
    renderMap();
    renderWhy();
    renderPairNote();
  } catch (err) {
    if (my !== siteSeq) return;
    $("upstream").replaceChildren();
    $("downstream").replaceChildren();
    $("mcs").replaceChildren();
    $("pair-note").textContent = err.message;
    $("why").hidden = true;
  }
}

function scheduleSites() {
  clearTimeout(siteTimer);
  siteTimer = setTimeout(refreshSites, 300);
}

function syncCterm() {
  const row = enz && enz.sites.find((s) => s.id === $("downstream").value);
  const hasTag = !!(row && row.as_downstream && row.as_downstream.has_tag);
  $("cterm-tag").disabled = !hasTag;
  $("cterm-tag-choice").classList.toggle("off", !hasTag);
  $("cterm-tag-note").textContent = hasTag
    ? "Adds only the bases that put the tag in frame with your protein."
    : "No His tag follows the chosen 3' site.";
  if (!hasTag && $("cterm-tag").checked) document.querySelector('input[name="cterm"][value="stop"]').checked = true;
}

function renderMap() {
  const box = $("mcs");
  const { start, end } = enz.window;
  const span = Math.max(1, end - start);
  const pct = (p) => `${Math.min(100, Math.max(0, ((p - start) / span) * 100))}%`;
  const up = $("upstream").value, down = $("downstream").value;
  const nodes = [h("div", { class: "mcs-bar" }), h("span", { class: "mcs-end l", text: pos1(start) }), h("span", { class: "mcs-end r", text: String(end) })];

  // annotated features that fall in the window, stacked in up to three lanes
  const lanes = [-1, -1, -1];
  const inWindow = vec.features.filter((f) => f.end > start && f.start < end && f.type !== "misc_feature" || /start codon|ATG/.test(f.label) && f.end > start && f.start < end);
  for (const f of inWindow.slice(0, 40)) {
    const lane = lanes.findIndex((last) => last <= f.start);
    if (lane < 0) continue;
    lanes[lane] = f.end;
    const bar = h("i", { class: `mcs-feat l${lane} ${f.type === "CDS" ? "cds" : ""} ${/start codon|^ATG$/.test(f.label) ? "start" : ""}`, title: `${f.label} (${pos1(f.start)} to ${f.end})` });
    bar.style.left = pct(f.start);
    bar.style.width = `max(3px, calc(${pct(f.end)} - ${pct(f.start)}))`;
    nodes.push(bar);
  }

  // place each label on the first row where it clears the previous one, using real text widths
  const width = box.clientWidth || 360;
  const rowEnd = [];
  const place = (px, text) => {
    const w = text.length * 6.6 + 8;
    for (let r = 0; r < MAX_LABEL_ROWS; r++) {
      if (rowEnd[r] === undefined || rowEnd[r] <= px - w / 2) { rowEnd[r] = px + w / 2; return r; }
    }
    return -1;
  };
  enz.sites.forEach((s) => {
    const role = s.id === up ? "up" : s.id === down ? "down" : "";
    const usable = s.as_upstream.ok || (s.as_downstream && s.as_downstream.ok);
    const off = !role && !usable;
    const tick = h("i", { class: `mcs-tick ${s.count > 1 ? "multi" : ""} ${role} ${off ? "off" : ""}`, title: siteLabel(s) });
    tick.style.left = pct(s.start);
    nodes.push(tick);
    const row = place(((s.start - start) / span) * width, s.enzyme);
    if (row >= 0) {
      const label = h("span", { class: `mcs-label ${role} ${off ? "off" : ""}`, text: s.enzyme, title: siteLabel(s) });
      label.style.left = pct(s.start);
      label.style.top = `${30 + row * 13}px`;
      nodes.push(label);
    }
  });
  box.style.height = `${40 + Math.max(1, rowEnd.length) * 13}px`;
  box.replaceChildren(...nodes);
  box.setAttribute("aria-label", `Restriction sites from ${pos1(start)} to ${end}. 5' site: ${up || "none"}. 3' site: ${down || "none"}.`);

  const names = inWindow.slice(0, 12).map((f) => f.label);
  $("mcs-key").textContent = names.length
    ? `In this window: ${names.join(", ")}. Dashed ticks are enzymes that cut more than once.`
    : "Dashed ticks are enzymes that cut more than once.";
}

function renderWhy() {
  const items = [];
  for (const s of enz.sites) {
    if (!s.as_upstream.ok) items.push([siteLabel(s), `as the 5' site: ${s.as_upstream.reason}`]);
    const d = s.as_downstream;
    if (d && !d.ok && !d.hidden && d.reason) items.push([siteLabel(s), `as the 3' partner of ${$("upstream").selectedOptions[0]?.text || "the 5' site"}: ${d.reason}`]);
  }
  $("why-list").replaceChildren(...items.map(([n, r]) => h("li", null, h("b", { text: n }), " ", r)));
  $("why").hidden = !items.length;
}

function renderPairNote() {
  const up = enz.sites.find((s) => s.id === $("upstream").value);
  const down = enz.sites.find((s) => s.id === $("downstream").value);
  if (!up) { $("pair-note").textContent = "No pair of sites in this window can clone this protein. Widen the window or change the start codon."; return; }
  const start = enz.start;
  const holdsStart = start != null && up.start <= start && start < up.end;
  let text = holdsStart
    ? `${up.enzyme} holds the start codon, so the vector's N-terminal sequence is replaced and the ATG in the site becomes your Met. Leave the first M out of your sequence.`
    : `Translation starts at the vector's start codon, so the vector sequence between it and ${up.enzyme} stays on as an N-terminal tag.`;
  if (!holdsStart && up.site.endsWith("ATG")) text += ` The ATG in ${up.enzyme} becomes the first Met of your protein, so leave your own first M out.`;
  if (down) {
    const cut = (s) => `${s.enzyme} cuts ${s.site} and leaves a ${s.overhang_kind} ${s.overhang} end`;
    text += ` ${cut(up)}; ${cut(down)}. The insert goes in one way round only.`;
  }
  $("pair-note").textContent = text;
}

/* ------------------------------------------------------------------ results */
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
      h("li", null, h("b", { text: "Choose a plasmid." }), " Use a bundled one, or upload your own GenBank, SnapGene or FASTA file."),
      h("li", null, h("b", { text: "Pick two restriction sites." }), " Every site is listed by position, and sites the plasmid has several of can be told apart."),
      h("li", null, h("b", { text: "Read the construct summary." }), " Fusion protein, plasmid size, a check digest, and the DNA to order.")),
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
  return h("div", null, h("dt", { text: label }), h("dd", { text: value }));
}

function sizeItem(label, value, note) {
  return h("div", null, h("dt", { text: label }), h("dd", null, bp(value), note ? h("small", { text: note }) : null));
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

function bandLine(title, fragments) {
  return h("p", null, `${title} `,
    fragments.flatMap((f, i) => [i ? " + " : "", h("b", { text: bp(f.length) }), " ", h("span", { class: "role", text: `(${ROLE_TEXT[f.role] || f.role})` })]));
}

function constructSummary(r) {
  const rep = r.report, s = rep.sizes;
  const sign = s.net_change_bp >= 0 ? "+" : "";
  const sizes = h("dl", { class: "sizes" },
    sizeItem("Plasmid before", s.vector_bp),
    sizeItem("Piece replaced", s.replaced_bp, `${r.upstream} to ${r.downstream}`),
    sizeItem("Backbone kept", s.backbone_bp),
    sizeItem("Insert, with its sites", s.insert_bp, `${bp(s.insert_after_digest_bp)} once cut`),
    sizeItem("Plasmid after", s.plasmid_bp, `${sign}${s.net_change_bp.toLocaleString("en-US")} bp`));

  const parts = h("ul", { class: "parts" }, rep.insert_parts.map((p) =>
    h("li", null, h("span", { class: "what", text: p.name }), h("span", { class: "seq", text: p.dna || bp(p.bp) }))));

  const pads = [];
  if (rep.frame.pad_5prime) pads.push(`5' ${rep.frame.pad_5prime} (${rep.frame.pad_5prime.length} base${rep.frame.pad_5prime.length > 1 ? "s" : ""})`);
  if (rep.frame.pad_3prime) pads.push(`3' ${rep.frame.pad_3prime} (${rep.frame.pad_3prime.length} base${rep.frame.pad_3prime.length > 1 ? "s" : ""})`);
  const padText = pads.length
    ? `Bases added so the reading frame lines up: ${pads.join(", ")}. A different choice of sites or C-terminus changes these.`
    : "No extra bases were needed to keep the reading frame.";

  const dg = rep.digest_check;
  const bands = h("div", { class: "bands" },
    bandLine(`Cut the plasmid with ${dg.enzymes.join(" + ")}:`, dg.vector),
    bandLine("Cut the finished plasmid the same way:", dg.final));

  const gone = rep.replaced_features.filter((f) => f.fully_removed);
  const essential = gone.filter((f) => f.essential);
  const flags = [];
  if (essential.length) flags.push(h("p", { class: "flag bad", text: `This choice removes ${essential.map((f) => f.label).join(", ")} from the plasmid. Check that you meant to replace this much.` }));
  const others = gone.filter((f) => !f.essential);
  if (others.length) flags.push(h("p", { class: "flag", text: `Removed with the replaced piece: ${others.map((f) => f.label).join(", ")}.` }));

  const pr = rep.protein;
  return h("section", { "aria-labelledby": "sum-h" },
    h("h2", { id: "sum-h", text: "Construct summary" }),
    sizes,
    h("h3", { class: "sub", text: "What the insert is made of" }), parts,
    h("p", { class: "help", text: padText }),
    h("h3", { class: "sub", text: "Check digest" }), bands,
    h("h3", { class: "sub", text: "Fusion protein" }),
    h("p", { text: `${pr.residues} residues, ${pr.mass_kda} kDa, pI ${pr.pi}. Your protein is ${pr.payload_residues} of them.` }),
    ...flags);
}

function renderResult(r) {
  const total = r.segments.reduce((n, s) => n + s.aa.length, 0);
  const verdict = h("div", { class: `verdict ${r.ok ? "pass" : "fail"}` },
    h("h2", { text: r.ok ? "Passed the round-trip check" : "Failed validation" }),
    h("p", { class: "where", text: `${r.vector}, ${r.upstream} to ${r.downstream}, ${r.payload_length} residues` }));

  const stats = h("dl", { class: "stats" },
    stat("Insert", bp(r.insert_bp)),
    stat("Plasmid", bp(r.plasmid_bp)),
    stat("CAI", r.metrics.cai == null ? "n/a" : r.metrics.cai.toFixed(2)),
    stat("GC", r.metrics.gc == null ? "n/a" : `${(r.metrics.gc * 100).toFixed(1)}%`),
    stat("Rare codons", String(r.metrics.n_rare ?? 0)));

  const map = h("section", { "aria-labelledby": "map-h" },
    h("div", { class: "map-title" }, h("h2", { id: "map-h", text: "Fusion protein map" }), h("p", { text: `${total} residues from the start codon` })),
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
    h("p", { class: "key", text: "Dark bases are your coding region. Lighter bases are the cloning sites and linker bases around it. The GenBank file is the whole finished plasmid with the vector's annotations." }));

  setResult(h("div", null, verdict, stats), constructSummary(r), map, issues, dna);
}

/* ------------------------------------------------------------------- submit */
function loadExample() {
  $("protein").value = GFP_EXAMPLE;
  checkProtein(false);
  scheduleSites();
  $("protein").focus();
}

async function onSubmit(e) {
  e.preventDefault();
  if (busy || !checkProtein(true)) { if (!busy) $("protein").focus(); return; }
  if (!vec) { renderError("Choose or upload a plasmid first."); return; }
  if (currentStart() == null) { renderError("Choose where translation starts first."); return; }
  if (!$("upstream").value || !$("downstream").value) { renderError("Choose a usable 5' and 3' site first."); return; }
  busy = true;
  const btn = $("submit");
  btn.disabled = true;
  btn.textContent = "Designing...";
  renderLoading();
  const form = new FormData($("design-form"));
  const body = {
    vector_id: vec.id, start: currentStart(), protein: form.get("protein"),
    upstream: form.get("upstream"), downstream: form.get("downstream"), cterm: form.get("cterm"),
    host: form.get("host"), strategy: form.get("strategy"), seed: form.get("seed"),
  };
  try {
    const data = await post("/api/design", body);
    renderResult(data);
    $("result").focus({ preventScroll: false });
  } catch (err) {
    renderError(err.name === "TypeError" ? "Could not reach the p2d server. Check that `p2d ui` is still running in your terminal." : err.message);
  } finally {
    busy = false;
    btn.disabled = false;
    btn.textContent = "Design insert";
  }
}

/* --------------------------------------------------------------------- boot */
function showHostNote() {
  const host = meta.hosts.find((x) => x.slug === $("host").value);
  $("host-note").textContent = host ? (host.verified ? "Codon table from a cited source." : "Codon table is an unverified placeholder.") : "";
}

function initForm() {
  vectorOptions();
  $("host").replaceChildren(...meta.hosts.map((x) => h("option", { value: x.slug, text: x.name })));
  $("strategies").replaceChildren(
    ...meta.strategies.map((s, i) =>
      h("label", { class: "choice" },
        h("input", { type: "radio", name: "strategy", value: s.value, checked: i === 0 }),
        h("b", { text: s.label }),
        h("span", { text: s.note }))));
  showHostNote();

  $("vector").addEventListener("change", onVectorSelect);
  $("upload-btn").addEventListener("click", () => $("vector-file").click());
  $("vector-file").addEventListener("change", onFile);
  $("paste-btn").addEventListener("click", () => {
    const open = $("paste-box").hidden;
    $("paste-box").hidden = !open;
    $("paste-btn").setAttribute("aria-expanded", String(open));
    if (open) $("paste-text").focus();
  });
  $("paste-use").addEventListener("click", () => {
    const text = $("paste-text").value.trim();
    if (!text) return setStatus("Paste a plasmid sequence first.", true);
    loadVector({ text, filename: "pasted" });
  });
  $("start").addEventListener("change", () => onStartChange(true));
  $("start-manual").addEventListener("input", () => { winTouched = false; scheduleSites(); });
  $("host").addEventListener("change", () => { showHostNote(); refreshSites(); });
  $("multi").addEventListener("change", refreshSites);
  $("upstream").addEventListener("change", refreshSites);
  $("downstream").addEventListener("change", () => { if (enz) { syncCterm(); renderMap(); renderPairNote(); } });
  document.querySelectorAll('input[name="cterm"]').forEach((r) => r.addEventListener("change", refreshSites));
  for (const id of ["win-from", "win-to"]) $(id).addEventListener("input", () => { winTouched = true; scheduleSites(); });
  $("win-all").addEventListener("click", () => { winTouched = true; $("win-from").value = 1; $("win-to").value = vec.length; refreshSites(); });
  $("win-orf").addEventListener("click", () => { winTouched = false; refreshSites(); });
  $("win-mcs").addEventListener("click", () => {
    if (!vec || !vec.mcs) return;
    winTouched = true;
    $("win-from").value = Math.max(1, pos1(vec.mcs.start) - 20);
    $("win-to").value = Math.min(vec.length, vec.mcs.end + 20);
    refreshSites();
  });

  $("protein").addEventListener("input", () => { checkProtein(false); scheduleSites(); });
  $("protein").addEventListener("keydown", (e) => {
    if (e.key === "Enter" && (e.metaKey || e.ctrlKey)) { e.preventDefault(); $("design-form").requestSubmit(); }
  });
  $("load-example").addEventListener("click", loadExample);
  $("design-form").addEventListener("submit", onSubmit);
}

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
  if (meta.bundled.length) loadVector({ bundled: meta.bundled[0].id });
}
document.addEventListener("DOMContentLoaded", boot);
