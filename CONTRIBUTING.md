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

## Architecture rules

Three rules the tests enforce or the reviews will:

1. **Nothing under `src/p2d/` imports a UI framework.** Enforced by
   `test_library_imports_no_ui_framework`. The eventual web UI is a shell over
   the library, never the place logic lives.
2. **Chains are always a list.** Do not add a single-chain fast path that
   bypasses the model.
3. **Never trust annotation labels.** Re-derive ORFs and sites from the
   sequence; use the record's features for reporting only.

## Scope

Open an issue before starting anything on the v0.2+ roadmap (the codon
optimiser, enzyme-pair search, Addgene integration, multi-chain linkage) —
these have design decisions attached that are worth settling in an issue
first.


## Testing against real plasmids

Real vectors are not in the repository because their licences differ. Tests that need them
read the folder named by `P2D_VECTOR_SAMPLES` and are skipped when it is unset:

```bash
P2D_VECTOR_SAMPLES=~/vectors pytest -q
```

The suite looks for `pET-28a_plus.dna` (SnapGene) and `U13853.gb`, `U13852.gb`, `U13850.gb`
(pGEX-4T-1, 3X, 2T from NCBI). Please do not commit plasmid files.
