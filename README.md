# p2d — protein-to-DNA reverse translation that knows about your vector

Reverse-translate a designed protein into DNA that is **actually cloneable into the vector you
entered** — correct reading frame through both junctions, restriction sites at
the ends, no unexpected stop codons, and a preview of the fusion protein you will
really get, scar residues included.

> **Status: v0.2 — frame engine + constrained codon optimiser.**
> Codons are chosen by a DFA-constrained dynamic programme, so the insert is
> provably free of both enzymes' recognition sites on either strand, including
> sites that would straddle the vector junctions.
> The bundled `ecoli_bl21` table is *E. coli* K-12 usage from
> `python_codon_tables`, used for BL21 as is standard — see the profile's
> `codon_usage.note`. Review a design before you order it.

## Two ways to use it

**A local website.** Run `p2d ui` and use p2d in your browser: pick or upload a plasmid,
choose two restriction sites on a map, and read the construct summary. For who will not code, there is nothing to
code in this way. It runs on your own computer (there is no hosted version).

**A Python module.** `import p2d` to script it. The website is a thin layer over this
module, so both give the same answer for the same input.

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

## The website

```bash
p2d ui
```

1. **Pick a plasmid.** Use the bundled test vector, one of your own files, an upload, or
   pasted text. GenBank (`.gb`), SnapGene (`.dna`) and FASTA are read; a plasmid stored on
   the opposite strand (ex. pET-28a(+)) is flipped automatically. Put files
   you do not want to re-upload in `~/.p2d/vectors` (or set `P2D_VECTOR_DIR`) and they show
   up in the list; nothing is copied or sent anywhere.
2. **Confirm where translation starts.** p2d offers the likely start codons with their
   evidence (an annotated start codon, a ribosome binding site, a T7 promoter) and you
   choose.
3. **Choose two restriction sites by position.** Every restriction site is listed, and
   an enzyme the plasmid cuts several times is shown as "2 of 3" and so on. Two sites work
   when no other site of either enzyme lies in the part you keep; each unavailable site
   comes with the reason (wrong way round, touching, same ends, forced by your protein).
   Not sure which pair? **Suggest site pairs** designs every pair in the window for your
   protein and groups them by the protein each gives (for example NcoI drops the vector's
   N-terminal tag, NdeI keeps it), fewest added residues first, with a **Use** button that
   fills in the two sites. A window with hundreds of pairs must be narrowed first.
4. **Add tags and linkers if you want them.** An N-terminal tag and linker, and a
   C-terminal linker and tag (His6, His8, FLAG, Myc, Strep-II, HA, V5; GS, G4S, (G4S)3,
   EAAAK) are built onto your protein and shown as their own parts in the summary.
5. **Choose the C-terminus.** A stop codon, or read on into the vector's own tag (His, FLAG,
   Myc, Strep-II, HA or V5; p2d adds
   the bases that put it in frame; some vectors, pET-28a(+) among them, have the tag one
   base off the start codon's frame).
6. **Set the extra bases for ordering.** A restriction site at the very end of a synthetic
   fragment often fails to cut, so p2d adds bases outside each site (six per side by default,
   NEB's [general guideline](https://www.neb.com/tools-and-resources/usage-guidelines/cleavage-close-to-the-end-of-dna-fragments);
   some enzymes need fewer). They are cut away and never reach the plasmid.
7. **Read the summary and take the DNA to order.** Fusion protein, plasmid size before and
   after, the piece replaced, the bases added for the reading frame, and a check digest
   (fragment sizes before and after). The order is checked in silico: it is cut with
   Biopython's cut positions and ligated into the vector, and must give the planned plasmid
   exactly. The whole finished plasmid is available as GenBank with the vector's annotations.

It runs only on your computer, binds to loopback, and needs no dependencies beyond the
library itself. The same steps are available from code, below.

## The Python module

### Quick start

```python
import p2d

# the bundled synthetic test plasmid; NdeI supplies the start Met, so leave a leading M out
a = p2d.run("KVFLDWINEAYQRGTRVLAEMAKRGDEFVKRLIAEGHDPFEVLKELGYSE", site="NdeI-XhoI")

print(a.ok)  # the round-trip validation passed
print(a.fusion_report())  # the protein you actually get, residue by residue
print(a.insert_dna)  # the DNA to order
a.to_record()  # the whole annotated plasmid, writable as GenBank
```

```
      1-19    vector (N-terminal)      MGSSHHHHHHSSGLVPRGS
     20-21    scar: 5' junction        HM
     22-71    insert                   KVFLDWINEAYQRGTRVLAEMAKRGDEFVKRLIAEGHDPFEVLKELGYSE
     72-73    scar: 3' junction        LE
     74-79    vector (C-terminal)      HHHHHH
  total: 79 aa
```

The `HM` and `LE` are the NdeI and XhoI sites being translated. They are residues that will be in your
expressed protein, and most tools never show them to you.
The N-terminal tag stays because NdeI lies downstream of this vector's start codon, which is
exactly the kind of thing that goes wrong when the vector is ignored.

### With your own plasmid

```python
from Bio import SeqIO

import p2d
from p2d.analysis import analyze_vector, confirm
from p2d.enzymes import find_cut_site, list_cut_sites, plan_cut
from p2d.report import construct_report
from p2d.vectorio import read_vector_file

# 1. read your plasmid (GenBank, SnapGene .dna or FASTA); it is flipped if stored backwards
analysis = analyze_vector(read_vector_file("my_plasmid.dna"))
print(analysis.best_start.position + 1, analysis.best_start.confidence)  # where translation starts
vector = confirm(analysis)  # accepts the start codon (pass start=... to choose your own)

# 2. every restriction site, by position (an enzyme with several sites is numbered)
for site in list_cut_sites(vector, window=(5160, 5215)):
    print(site.id, site.label)

# 3. choose two sites and design
site = plan_cut(vector, find_cut_site(vector, "BamHI:5166"), find_cut_site(vector, "XhoI:5206"))
plan = p2d.ConstructPlan.single_chain("KVFLDWINEAYQRGTRVLAEMAKRGDEFVKRLIAEGHDPFEVLKELGYSE")
design = p2d.design_insert(plan, vector, site, p2d.load_host("ecoli_bl21"), cterm=p2d.CTerm.TAG)
assembly = p2d.assemble(design)

# 4. read the result
print(assembly.ok, assembly.protein)
print(construct_report(assembly)["sizes"])  # plasmid before/after, replaced piece, insert, ...
SeqIO.write(assembly.to_record(), "finished_plasmid.gb", "genbank")

# 5. the sequence to order: the insert plus extra bases outside each site (6 per side by default)
order = p2d.make_order(assembly)
print(order.dna)
print(order.ligation.message)  # the order was cut and ligated into the vector in silico
```

`plan_cut` refuses a pair that would not work and says why (wrong way round, touching, same
ends, or another site of the enzyme in the part of the plasmid you keep).

### The main entry points

| Task | Function |
|---|---|
| One call, bundled plasmid | `p2d.run(protein, site=...)` |
| Read a plasmid | `p2d.vectorio.read_vector_file / read_vector_text / read_vector_bytes` |
| Find the start codon | `p2d.analysis.analyze_vector`, then `confirm` |
| List and choose sites | `p2d.enzymes.list_cut_sites`, `find_cut_site`, `plan_cut` |
| Design and validate | `p2d.design_insert` (`CTerm`, `CodonStrategy`), `p2d.assemble` |
| Sizes, check digest, protein mass | `p2d.report.construct_report` |
| Just the codon optimiser | `p2d.optimize.optimize_coding` |

The result of `assemble` has `.ok`, `.protein`, `.seq`, `.insert_dna`, `.fusion_report()`
and `.to_record()`.

## Why this exists

Codon-optimisation tools give you a coding sequence. They do not tell you what
protein comes out the other end after you ligate it into a vector, because they
never look at the vector. That last step is where considerations are needed:
an off-by-one at the junction, a C-terminal His-tag that is out of frame, a
second His-tag you did not realise the vector already supplied, a stop codon
you left in that silently deletes the tag.

## CLI

```bash
p2d hosts                                  # list host profiles
p2d vectors                                # list bundled vectors
p2d design KVFLDWINEAYQ... --site NdeI-XhoI --gb out.gb
p2d ui                                     # the local website
```

## What it checks

- **Frame**: The finished plasmid is translated from the vector's own
  start codon, and a design whose assumed start differs from the real one is rejected
- **Pad bases** added automatically when a pair of sites would shift the frame, chosen to
  avoid stop codons and new restriction sites, and reported to you
- **Cloning sites kept out of the coding DNA**, on both strands and across the junctions
  with the vector, by a search; a site your protein forces is reported, not hidden
- **Host preferences**: internal Shine-Dalgarno, Chi sites and long homopolymer runs are
  avoided when a path exists, and any that the protein forces are named in a warning
- **The order itself**: extra bases outside each site that create no new cloning site or host
  motif, and an in-silico digest and ligation of the whole order that must reproduce the plasmid
- **Which sites you can use**: order, spacing, compatible ends, and whether a second site of
  the same enzyme would cut the backbone you keep
- **The C-terminus**: a stop codon, or a read into the vector's His tag with the bases that
  put it in frame
- **Duplicate tags** between insert and vector
- **GC content and repeats**: after the codons are chosen, p2d swaps synonymous codons to bring
  overall GC, GC in sliding windows, and direct and inverted repeats inside limits you can
  change (website: *GC and repeat limits*). A swap is refused if it would add a cloning site
  or a host motif; whatever the protein itself forces stays and is reported
- **Round-trip**: the payload must survive intact and in frame, with no internal stop and a
  terminating stop

## What it does not do yet

p2d makes promises about the *sequence*: the reading frame, which sites are in the DNA, the
fusion protein you get, and the sizes. It does not yet promise anything about how much protein
you will make.

- **No experimental validation.** Codons follow a published usage table (the *E. coli* K-12
  table stands in for BL21). That they raise expression in your hands has not been tested.
- **Not implemented:** 5′ mRNA folding, and Dam/Dcm
  methylation and double-digest buffers. p2d also does not model any synthesis vendor's rules
  (length, GC or complexity limits): those differ between vendors and change, so they are out of
  scope. Check a design against the rules of whoever makes it.
- **Host limits are defaults, not citations.** The homopolymer lengths, motifs, GC range, GC
  window and repeat length in the host profile come from the project template, not from a
  reference. They never block a design (what cannot be met is reported), but do not read them
  as vetted thresholds. Swapping codons to meet them lowers CAI; the page shows both.
- **Scope.** One protein chain at a time, and restriction-site cloning only (no Gibson or
  Golden Gate). Tested on bacterial vectors: pET-28a(+), pGEX-4T-1/3X/2T and a synthetic one.
  Eukaryotic vectors are not supported.

Review a design before you order it.

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
   

## Roadmap

| version | adds |
|---|---|
| **0.1** | data model, frame engine, round-trip validation, CLI ✅ |
| **0.2** | Aho-Corasick DFA + DP codon optimiser (both strands, IUPAC, junction-aware, unavoidable sites reported); cited codon table; `MAX_CAI` and seeded `WEIGHTED_SAMPLE` strategies ✅ |
| **0.2.x** | the local website; reading your own plasmid (GenBank, SnapGene, FASTA); start-codon detection; picking two sites by position, with suggested pairs ranked by what they add to the protein; C-terminus choice; tags and linkers; construct report; GC and repeat refinement with adjustable limits; host motifs and homopolymer limits as soft preferences; the DNA to order with extra bases outside the sites and an in-silico digest and ligation check ✅ |
| 0.3 | first PyPI release |
| 0.4 | 5′ mRNA structure (ViennaRNA, optional), batch design from a FASTA, the `HARMONIZED` strategy, more host profiles |
| 0.5 | Dam/Dcm methylation and buffers, validation on more real vectors, a hosted version of the website if it can run in the browser, a Korean interface |
| later | Gibson and Golden Gate assembly, multi-chain constructs, eukaryotic vectors, a vector-free mode, Addgene |

## Bundled vector

`pTEST1` is a **synthetic** 1149 bp pET-like vector used by the test suite.
Do not order it. Use your own plasmids: `read_vector_file()` from
Python, or the website's upload and `~/.p2d/vectors` folder.

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

---

Developed with AI assistance (Claude).
