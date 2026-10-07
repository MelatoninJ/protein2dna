# Contributing

```bash
git clone https://github.com/MelatoninJ/protein2dna && cd protein2dna
uv venv && source .venv/bin/activate
uv pip install -e ".[dev]"
pytest -q && ruff check .
```

## The easiest useful contribution: a host profile

Host-specific behaviour is data, not code. Adding an expression host means
adding **one YAML file** to `src/p2d/hosts/` and nothing else.

1. Copy `src/p2d/hosts/ecoli_bl21.yaml`.
2. Replace the codon usage table with a cited one — Kazusa CUTG or HIVE-CUT —
   and fill in `codon_usage.source` with the exact table and download date.
   Set `verified: true` only when you have done this.
3. Set the host-specific `forbidden_motifs`. For eukaryotic hosts this
   normally means polyadenylation signals (`AATAAA`), AU-rich elements
   (`ATTTA`), and splice donor/acceptor consensus sequences. For bacteria it
   means internal Shine–Dalgarno sequences and Chi sites.
4. Add the slug to `test_structure.py::test_every_residue_has_a_codon` —
   `test_codon_table_is_self_consistent` and `test_codon_table_is_complete`
   will catch most transcription errors for you.

Wanted: *S. cerevisiae*, *P. pastoris*, human/CHO, *S. frugiperda*,
*N. benthamiana*.

## Adding a vector

Bundled vectors live in `src/p2d/vectors.py`. A GenBank file alone is not
enough — you also need `InsertionSite` definitions, because no GenBank record
reliably tells you where the reading frame starts for a given enzyme pair.
Each site needs the enzyme pair, the start-codon source, and (for
`StartSource.VECTOR`) the coordinate of the vector's ATG.

Every bundled vector must come with a golden test in
`tests/test_roundtrip.py::GOLDEN` giving the exact expected fusion protein.
A vector without a golden test will not be merged.

## Scope

Open an issue before starting anything on the v0.3+ roadmap (GC windows and mRNA
folding, automatic enzyme-pair ranking, Addgene integration, multi-chain linkage).
These have design decisions attached that are worth settling in an issue first.
The rules a change must respect are listed under "Architecture and rules" below.

## Testing against real plasmids

Real vectors are not in the repository because their licences differ. Tests that need them
read the folder named by `P2D_VECTOR_SAMPLES` and are skipped when it is unset:

```bash
P2D_VECTOR_SAMPLES=~/vectors pytest -q
```

The suite looks for `pET-28a_plus.dna` (SnapGene) and `U13853.gb`, `U13852.gb`, `U13850.gb`
(pGEX-4T-1, 3X, 2T from NCBI). Please do not commit plasmid files.

## Architecture and rules

Read this before changing anything.

### What this is

`p2d` reverse-translates a designed protein (RFdiffusion / ProteinMPNN output,
nanobodies, binders) into DNA that is actually cloneable into a chosen vector.
The differentiator is not codon optimisation — it is that the vector is a
first-class input, so the tool can tell the user **what protein actually comes
out** after ligation, scar residues and all.

**Current state: v0.2.** Data model, frame engine, assembly, round-trip
validation and the constrained codon optimiser (`p2d.optimize`) are done and
tested. `reverse_translate` in `frame.py` is kept only as the naive baseline.

### Commands

```bash
uv pip install -e ".[dev]"   # or: pip install -e ".[dev]"
pytest -q                    # all tests must pass
ruff check . && ruff format --check .
python tools/build_test_vector.py   # regenerates src/p2d/data/pTEST1.gb
python tools/build_codon_table.py --write   # regenerates the ecoli_bl21 codon_usage block
```

### Architecture

```
src/p2d/
  model.py      Part / Chain / ConstructPlan, tag+linker+signal-peptide library
  vector.py     Vector, InsertionSite, enzyme search, ORF re-detection
  host.py       HostProfile loader (YAML -> dataclass)
  hosts/*.yaml  one file per expression host; this is where host behaviour lives
  frame.py      THE FRAME ENGINE. Junction arithmetic, pad selection,
                fusion preview; calls optimize_coding() for the coding region
  enzymes.py    every cloning-grade site of a vector by position (list_cut_sites),
                the rules for a chosen pair of *sites* (plan_cut), and the older
                name-based picker for enzymes that cut once (candidate_enzymes)
  report.py     sizes, check-digest fragments, frame pads and protein stats of a design
  vectorio.py   read a plasmid from GenBank, SnapGene .dna, FASTA or pasted text, and
                reverse-complement it if the expression cassette is on the minus strand
  analysis.py   start-codon candidates with evidence and confidence; reports an annotated
                MCS as a landmark only (the MCS is never guessed)
  issues.py     Issue / Severity (shared so optimize/ need not import frame)
  optimize/     automaton.py  Aho-Corasick DFA (IUPAC, both strands)
                core.py       DP over (residue, DFA state); strategies; metrics
  assemble.py   splice insert into vector, run the hard gates
  vectors.py    bundled vector registry + InsertionSite definitions
  cli.py        thin typer wrapper, no logic (`p2d ui` starts the web UI)
  data/         pTEST1.gb (synthetic test vector)
src/p2d_ui/     Local web UI. A sibling package, never imported by p2d except
                lazily from `p2d ui`. stdlib server, no build step.
  service.py    validation + JSON shaping; the only place that calls the library
  server.py     loopback-only http.server (Host-header check, CSP, size limits)
  static/       index.html, app.css, app.js (vanilla; server text via textContent only)
```

Data flow:

```
ConstructPlan + Vector + InsertionSite + HostProfile
        -> design_insert()  -> InsertDesign   (frame.py)
        -> assemble()       -> Assembly       (assemble.py)
```

#### The core arithmetic

```
final  = vector[:up.start] + insert + vector[down.end:]
insert = up_site + pad_up + coding + stop? + pad_down + down_site

pad_up   = (3 - (up.start + len(up_site) - orf_start) % 3) % 3
pad_down = (-len(down_site)) % 3
```

`orf_start` comes from one of two sources (`StartSource`): the vector's own ATG
(`VECTOR`, e.g. an N-terminal His-tag upstream of the MCS) or the ATG inside the
upstream recognition site (`UPSTREAM_SITE`, e.g. NdeI `CAT|ATG`,
`start_offset_in_site=3`).

### Invariants — do not break these

1. **Nothing under `src/p2d/` imports a UI framework.** Enforced by
   `tests/test_structure.py::test_library_imports_no_ui_framework`. The web UI
   (v0.5) will be a shell over the library, never the place logic lives.
2. **`ConstructPlan.chains` is always a list.** Single-chain is `len == 1`. Do
   not add a single-chain fast path. The multi-chain strategies in
   `LinkageStrategy` raise `NotImplementedError` on purpose.
3. **Never trust GenBank annotation labels.** ORFs and restriction sites are
   re-derived from the sequence. Features are for reporting only.
4. **The round-trip gate is sacred.** `assemble()` translates the assembled
   plasmid from the vector's own start codon and checks the payload survives
   intact, in frame, with a terminating stop. Never weaken, skip or
   special-case this.
5. **Host behaviour is data, not branches.** No `if host == "ecoli"`. New host
   = one YAML file in `src/p2d/hosts/`.

### Public API — breaking changes need a reason

```python
p2d.run(payload, *, vector=None, site="NdeI-XhoI", host="ecoli_bl21")
p2d.design_insert(plan, vector, site, host) -> InsertDesign
p2d.assemble(design) -> Assembly
p2d.load_bundled(name="pTEST1") -> Vector
p2d.load_host(slug_or_path) -> HostProfile
Assembly.ok / .protein / .seq / .insert_dna / .fusion_report() / .to_record()
```

### Known placeholders — do not treat as real data

- **Pairs are chosen as two individual sites** (`CutSite`, id like `XhoI:5206`). A second
  site of either enzyme inside the replaced piece is harmless; one in the kept backbone
  is refused (`backbone_conflicts`) because a complete digest would cut there too.
  `Vector.cloning_region` only serves the older name-based `candidate_enzymes`.
- **The C-terminus is a choice** (`CTerm`): `STOP`, `TAG` (reads into the vector's His
  tag, adding whatever bases align it) or `VECTOR_FRAME` (keep the vector's frame, the
  library default). In pET-28a(+) XhoI and the tag are one base off the start codon's
  frame, so `VECTOR_FRAME` reads vector junk; the UI defaults to `STOP`.
- **The user's own vectors** live in `P2D_VECTOR_DIR` (default `~/.p2d/vectors`) or are
  uploaded; they are held in memory only and never written to the project.
- **Real vectors are not in the repo** (SnapGene/Addgene/vendor licences differ). Tests
  that need them read the folder named by `P2D_VECTOR_SAMPLES` and skip when it is unset.
- **`Vector.expression_start` is curated, not read from annotations.** It is the ATG
  the ribosome really uses; `enzymes.py` and the 3' pad arithmetic hang on it. A vector
  without it cannot use the enzyme picker.
- **Start codon vs site ATG.** A site like NdeI (CAT|ATG) has an ATG, but the vector's
  own start codon upstream still initiates, so the tag stays. Only a site that *holds*
  the vector's start codon (NcoI, CC|ATG|G) removes the tag. `assemble()` enforces this
  with the `wrong_start` gate: it translates from `Vector.expression_start` and refuses a
  design whose declared start gives a different protein.
- **The `ecoli_bl21` codon table is K-12, not BL21.** It is generated by
  `tools/build_codon_table.py` from `python_codon_tables` (table
  `e_coli_316407`); never edit the block by hand and **never fabricate codon
  frequencies**. The K-12-for-BL21 substitution is recorded in
  `codon_usage.note`; keep it recorded.
- **Host `forbidden_motifs` (internal Shine-Dalgarno, Chi) are not yet passed to
  the optimiser** — only the two enzymes' sites are. Do not describe them as
  enforced.
- **`pTEST1` is a synthetic 1149 bp test vector**, not a real plasmid. It mimics
  pET-28a topology so the frame tests are meaningful. Do not describe it as
  pET-28a anywhere.

### Testing conventions

- `tests/test_roundtrip.py` holds the **golden set**: exact expected fusion
  proteins per (vector, site, payload). These assert at the *protein* level, so
  changing codon choice must not break them. If a change breaks a golden test,
  the change is wrong until proven otherwise.
- `tests/test_structure.py` holds vector parsing, host loading and the
  architecture-boundary tests.
- Any new bundled vector needs a golden test. No exceptions.

### Style

Ruff, line length 100, `select = ["E","F","I","UP","B","SIM"]`. Type hints on
public functions. Docstrings explain *why*, not *what*. Comments are for
non-obvious biology or non-obvious algorithmic choices — not narration.

### Roadmap

| version | adds |
|---|---|
| 0.1 | data model, frame engine, round-trip validation, CLI — **done** |
| 0.2 | Aho–Corasick DFA + DP codon optimiser; verified codon tables — **done** |
| 0.3 | GC windows, 5′ ΔG (ViennaRNA, optional dep), repeats, synthesis limits |
| 0.4 | enzyme-pair search, Dam/Dcm methylation, buffer compatibility |
| 0.5 | Addgene API, Gibson assembly, Streamlit UI |
| 1.0 | multi-chain linkage, mammalian hosts, Golden Gate |

UI rules (`tests/test_ui.py` enforces some): no logic in `p2d_ui` beyond input
validation and shaping; no `innerHTML`; no em or en dashes in UI text; one accent
colour, one radius; both light and dark themes; honour `prefers-reduced-motion`.

Next: v0.3 (non-local constraints as a local-search stage after the DP, keeping the DFA as a hard constraint).
