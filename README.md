# p2d — protein-to-DNA reverse translation that knows about your vector

Reverse-translate a designed protein (RFdiffusion / ProteinMPNN output, a
nanobody, a binder) into DNA that is **actually cloneable into the vector you
picked** — correct reading frame through both junctions, restriction sites at
the ends, no surprise stop codons, and a preview of the fusion protein you will
really get, scar residues included.

> **Status: v0.2 — frame engine + constrained codon optimiser.**
> Codons are chosen by a DFA-constrained dynamic programme, so the insert is
> provably free of both enzymes' recognition sites on either strand, including
> sites that would straddle the vector junctions. The default strategy samples
> codons above a usage floor (seeded, reproducible); `MAX_CAI` is available.
> **Not done yet:** GC windows, 5′ mRNA folding, repeat and synthesis checks
> (v0.3), and the host's `forbidden_motifs` are not yet fed to the optimiser.
> The bundled `ecoli_bl21` table is *E. coli* K-12 usage from
> `python_codon_tables`, used for BL21 as is standard — see the profile's
> `codon_usage.note`. Review a design before you order it.

## Install

p2d is not on PyPI yet; install it from source (Python 3.10 or newer):

```bash
git clone https://github.com/MelatoninJ/protein2dna
cd protein2dna
python -m venv .venv && source .venv/bin/activate
pip install -e .
p2d --help
```

For development use `pip install -e ".[dev]"`, then `pytest -q`.

## Web UI

```bash
pip install -e .
p2d ui            # opens http://127.0.0.1:8765 in your browser
```

1. **Pick a plasmid.** Use the bundled test vector, one of your own files, an upload, or
   pasted text. GenBank (`.gb`), SnapGene (`.dna`) and FASTA are read; a plasmid stored on
   the opposite strand (SnapGene's own pET-28a(+) is) is flipped automatically. Put files
   you do not want to re-upload in `~/.p2d/vectors` (or set `P2D_VECTOR_DIR`) and they show
   up in the list; nothing is copied or sent anywhere.
2. **Confirm where translation starts.** p2d offers the likely start codons with their
   evidence (an annotated start codon, a ribosome binding site, a T7 promoter) and you
   choose.
3. **Choose two restriction sites by position.** Every site in the window is listed, and
   an enzyme the plasmid cuts several times is shown as "2 of 3" and so on. Two sites work
   when no other site of either enzyme lies in the part you keep; each unavailable site
   comes with the reason (wrong way round, touching, same ends, forced by your protein).
4. **Choose the C-terminus.** A stop codon, or read on into the vector's His tag (p2d adds
   the bases that put it in frame; some vectors, pET-28a(+) among them, have the tag one
   base off the start codon's frame).
5. **Read the construct summary.** Fusion protein, plasmid size before and after, the
   piece replaced, the bases added for the reading frame, a check digest (fragment sizes
   before and after) and the DNA to order, plus the whole finished plasmid as GenBank with
   the vector's annotations carried over.

It runs only on your computer, binds to loopback, and needs no dependencies beyond the
library itself.

## Why this exists

Codon-optimisation tools give you a coding sequence. They do not tell you what
protein comes out the other end after you ligate it into pET-28a, because they
never look at the vector. That last step is where constructs actually fail:
an off-by-one at the junction, a C-terminal His-tag that is out of frame, a
second His-tag you did not realise the vector already supplied, a stop codon
you left in that silently deletes the tag.

p2d makes the vector a first-class input.

## Example

```python
import p2d

a = p2d.run("MKVFLDWINEAYQRGTRVLAEMAKRGDEFVKRLIAEGHDPFEVLKELGYSE")

print(a.ok)  # round-trip validation passed
print(a.fusion_report())  # what you actually get, residue by residue
print(a.insert_dna)  # order this
a.to_record()  # full annotated plasmid, writable as GenBank
```

```
      1-1     scar: 5' junction        M
      2-52    insert                   MKVFLDWINEAYQRGTRVLAEMAKRG...
     53-54    scar: 3' junction        LE
     55-60    vector (C-terminal)      HHHHHH
  total: 60 aa
```

The `LE` is the XhoI site (CTCGAG) being translated. It is not noise — it is
two real residues that will be in your purified protein, and most tools never
show it to you.

## CLI

```bash
p2d hosts                                  # list host profiles
p2d vectors                                # list bundled vectors
p2d design MKVFLDWINEAYQ... --site NdeI-XhoI --gb out.gb
```

## What it checks

- **Frame** through both junctions, from the vector's own start codon
- **Pad bases** inserted automatically when an enzyme pair would otherwise
  shift the frame, chosen to avoid stop codons and new restriction sites
- **Internal sites** of the cloning enzymes, removed by synonymous swap
- **Stop codon policy** — omitted when the vector supplies a C-terminal tag,
  added (as a double stop) when it does not
- **Duplicate tags** between insert and vector
- **Round-trip**: the assembled plasmid is translated from the vector's start
  codon and the result must contain the payload intact, with no frameshift,
  no internal stop, and a terminating stop

## Design rules

Three, and they are load-bearing:

1. **Chains are a list, always.** Single-chain is `len(chains) == 1`. The
   multi-chain strategies (Duet vectors, polycistronic, 2A polyproteins) are
   registered in `LinkageStrategy` and raise `NotImplementedError` rather than
   requiring a model rewrite later.
2. **Hosts are data, not branches.** A new expression host is one YAML file in
   `src/p2d/hosts/` — see `CONTRIBUTING.md`.
3. **Never trust annotation labels.** ORFs and sites are re-derived from the
   sequence; the record's own features are used only for reporting.

And one more, enforced by a test: nothing under `src/p2d/` may import a UI
framework. The library has to stay usable from a batch script.

## Roadmap

| version | adds |
|---|---|
| **0.1** | data model, frame engine, round-trip validation, CLI ✅ |
| **0.2** | Aho-Corasick DFA + DP codon optimiser (both strands, IUPAC, junction-aware, unavoidable sites reported); cited codon table; `MAX_CAI` and seeded `WEIGHTED_SAMPLE` strategies (`HARMONIZED` not implemented) ✅ |
| 0.3 | GC windows, 5′ ΔG (ViennaRNA), repeats, synthesis manufacturability |
| 0.4 | enzyme-pair search, Dam/Dcm methylation, buffer compatibility |
| 0.5 | Addgene API, Gibson assembly, Streamlit web UI |
| 1.0 | multi-chain linkage strategies, mammalian host profiles, Golden Gate |

## Bundled vector

`pTEST1` is a **synthetic** 1149 bp pET-like vector used by the test suite. It
is not a real plasmid; do not order it. Real vectors come from
`Vector.from_genbank()`.

## Licence

MIT, see `LICENSE`. ViennaRNA is an optional extra (`pip install protein2dna[fold]`) so that
the core package stays permissively licensed and installs without a compiler.

### Third-party data

* **Codon usage** comes from [`python_codon_tables`](https://pypi.org/project/python_codon_tables/)
  (CC0, derived from Kazusa CUTG). The bundled `ecoli_bl21` profile uses the *E. coli* K-12
  table for BL21, which is standard practice and is recorded in the profile.
* **Plasmid files are not redistributed.** `pTEST1` is synthetic. SnapGene, Addgene and vendor
  maps carry their own terms (SnapGene's, for example, require attribution and a licence for
  commercial use), so keep them on your own computer: put them in `~/.p2d/vectors` or upload
  them in the UI, and they never leave it.
* Depends on Biopython (Biopython License), PyYAML and Typer (MIT).
