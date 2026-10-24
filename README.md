# p2d — protein-to-DNA reverse translation that knows about your vector

Reverse-translate a designed protein into DNA that is **actually cloneable into the vector you
entered** — correct reading frame through both junctions, restriction sites at
the ends, no unexpected stop codons, and a preview of the fusion protein you will
really get, scar residues included.

> **Status: v0.2 — frame engine + constrained codon optimiser.**
> Codons are chosen by a DFA-constrained dynamic programme, so the insert is
> provably free of both enzymes' recognition sites on either strand, including
> sites that would straddle the vector junctions.
> The bundled codon tables come from `python_codon_tables` (see *Hosts* below). The
> `ecoli_bl21` table is *E. coli* K-12 usage, used for BL21 as is standard; each profile's
> `codon_usage.note` says what its table is. Review a design before you order it.

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
   choose. A plasmid with no annotation at all is read from its sequence: an ATG that follows a
   ribosome-binding-site-like motif and is part of an NcoI or NdeI site, or lies behind a T7
   promoter or lac operator, is offered (pTrc99A, pBAD24, pET-3a/11a/30a and pColdI as published
   at NCBI). A vector whose own sequence holds no start codon (pKK223-3, pT7-5, pBAD18/33, where
   the insert brings the ATG) is reported as having none, and you give the position yourself.
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

8. **Many proteins at once.** Open **Design many proteins at once**, paste FASTA text or load a
   file, and press **Design all**. Every protein uses the plasmid, sites, C-terminus, tags,
   host, limits and extra bases chosen above. A protein that cannot be cloned is named in the
   table with the reason and the others carry on. Download the summary (CSV), the DNA to order
   (FASTA), or everything, with one annotated GenBank file per construct, as a zip.

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

### Many proteins at once

```python
import p2d
from p2d.batch import design_batch, read_proteins, write_batch
from p2d.enzymes import resolve_site

vector = p2d.load_bundled()
proteins = read_proteins(open("binders.fasta").read())  # one FASTA record per protein
up, down = resolve_site(vector, "NdeI"), resolve_site(vector, "XhoI")
result = design_batch(proteins, vector, up, down, p2d.load_host("ecoli_bl21"))

print(result.n_ok, "of", len(result.entries), "designed")
for entry in result.entries:
    print(entry.name, entry.ok, entry.message)  # message is None when it worked
write_batch(result, "results")  # summary.csv, orders.fasta, inserts.fasta, constructs/*.gb
```

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
p2d batch binders.fasta --out results/     # every protein in a FASTA file, into one vector
p2d ui                                     # the local website
```

`p2d batch` takes the same choices as the website: `--vector` (a GenBank, SnapGene or FASTA
plasmid) with `--start` unless the file shows the start codon outright, `--upstream` and
`--downstream` (an enzyme name when it cuts once, or an id such as `XhoI:5206`), `--host`,
`--cterm`, `--n-part` and `--c-part` for tags and linkers, `--clamp-5` and `--clamp-3`. It
writes `summary.csv`, `orders.fasta`, `inserts.fasta`, `constructs/<name>.gb` and, if anything
failed, `failed.txt`. It exits 1 when some protein failed (the rest are still written) and 2
for a mistake in the options.

## What it checks

- **Frame**: The finished plasmid is translated from the vector's own
  start codon, and a design whose assumed start differs from the real one is rejected
- **Pad bases** added automatically when a pair of sites would shift the frame, chosen to
  avoid stop codons and new restriction sites, and reported to you
- **Cloning sites kept out of the coding DNA**, on both strands and across the junctions
  with the vector, by a search; a site your protein forces is reported, not hidden
- **Host preferences**: a host's own motifs (internal Shine-Dalgarno sequences in *E. coli* and
  *B. subtilis*, Chi sites in *E. coli*; polyadenylation signals and AU-rich elements in human
  and mouse) and runs of one base longer than the host's own genes contain are avoided when a
  path exists, and any that the protein forces are named in a warning
- **The order itself**: extra bases outside each site that create no new cloning site or host
  motif, and an in-silico digest and ligation of the whole order that must reproduce the plasmid
- **Which sites you can use**: order, spacing, compatible ends, and whether a second site of
  the same enzyme would cut the backbone you keep
- **The C-terminus**: one stop codon (the host's most used; two in a row on request), or a read
  into the vector's His tag with the bases that put it in frame. In *E. coli* the base after
  the stop codon changes how well translation ends, so p2d reports it
- **Sites that touch**: two sites with no base between them are offered with a warning; sites
  that share bases are refused
- **Duplicate tags** between insert and vector
- **GC content and repeats**: after the codons are chosen, p2d swaps synonymous codons to bring
  overall GC, GC in sliding windows, and direct and inverted repeats inside limits you can
  change (website: *GC and repeat limits*). A swap is refused if it would add a cloning site
  or a host motif; whatever the protein itself forces stays and is reported
- **Round-trip**: the payload must survive intact and in frame, with no internal stop and a
  terminating stop

## Checked against published constructs

p2d was checked against plasmids other people built. Given the empty pET-28a(+), the two sites and
the inserted protein, it rebuilt the same fusion protein for all 9 comparable constructs on NCBI (8
genes, 6 pairs of sites), three of them down to the exact plasmid size. 15 records cannot test it
(their reason is given). The test covers the junctions and the reading frame, not expression,
and only one vector family. See [VALIDATION.md](VALIDATION.md).

## What it does not do yet

p2d makes promises about the *sequence*: the reading frame, which sites are in the DNA, the
fusion protein you get, and the sizes. It does not yet promise anything about how much protein
you will make.

- **No experimental validation.** Codons follow a published usage table for the species (the
  *E. coli* K-12 table stands in for BL21; the human, mouse, yeast, *B. subtilis*, fly and worm
  tables are whole-genome usage, not tuned to a cell line, strain or expression level). That
  they raise expression in your hands has not been tested.
- **mRNA folding and its energy are deliberately not included.** The evidence that 5′ structure
  limits translation in *E. coli* is real (see [Kudla 2009](https://pmc.ncbi.nlm.nih.gov/articles/PMC3902468) and
  [Goodman 2013](https://pubmed.ncbi.nlm.nih.gov/24072823/)), but it is
  correlational and built on a few reporters, the sources use different software and windows,
  and none gives a free-energy cut-off, so p2d would have nothing citable to check a design
  against. Choose the codons near the start with that in mind if it matters for your protein.
- **Not implemented:** Dam/Dcm
  methylation and double-digest buffers. p2d also does not model any synthesis vendor's rules
  (length, GC or complexity limits): those differ between vendors and change, so they are out of
  scope. Check a design against the rules of whoever makes it.
- **Eukaryotic host limits are defaults, not citations.** The bacterial profiles are sourced (see
  *Hosts*). The human, mouse, yeast, fly and worm profiles, which are in development, still use
  homopolymer lengths, motifs, GC range, GC window and repeat length from the project template.
  They never block a design (what cannot be met is reported), but do not read them as vetted.
  Swapping codons to meet limits lowers CAI; the page shows both.
- **Scope.** One protein chain at a time, and restriction-site cloning only (no Gibson or
  Golden Gate). Tested on bacterial vectors: pET-28a(+), pGEX-4T-1/3X/2T, pUC19 and a
  synthetic one, plus a mammalian vector (mOrange2-C1) for the reading frame, the pad bases and
  the fusion protein, which agree with a hand design. Choose the host profile that matches
  where the vector is used; the default is *E. coli*.

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
   

## Hosts

`p2d hosts` lists the profiles. Each one is a YAML file generated from a table in
`python_codon_tables`, never typed by hand.

| profile | species | table | status |
|---|---|---|---|
| `ecoli_bl21` (default) | *Escherichia coli* BL21(DE3) | K-12 usage | released |
| `bsubtilis` | *Bacillus subtilis* | *B. subtilis* | released |
| `hsapiens` | human | *Homo sapiens* | in development |
| `mmusculus` | mouse | *Mus musculus* | in development |
| `scerevisiae` | baker's yeast | *Saccharomyces cerevisiae* | in development |
| `dmelanogaster` | fruit fly | *Drosophila melanogaster* | in development |
| `celegans` | roundworm | *Caenorhabditis elegans* | in development |

**Bacterial hosts are released; every eukaryotic profile is in development.** A profile in
development works, but the website and `p2d hosts` mark it, and every design made with it
carries a warning (`host_in_development`). The reason is that nothing beyond a reading-frame
check against one mammalian vector has been done for eukaryotes: their vectors, their start
and polyadenylation signals and their expression behaviour are not validated here.

**What is cited in the bacterial profiles.** The codon table is *E. coli* K-12 usage (for BL21)
or *B. subtilis* usage from `python_codon_tables`. Everything else is sourced or measured:

* *E. coli* rare codons AGA, AGG, AUA, CUA, GGA, CCC, CGG: the codons the Rosetta 2 strains carry
  extra tRNA genes for (Novagen/Merck, cat. 71397).
* Internal Shine-Dalgarno-like sequence `AGGAGG` (both bacteria): Li, Oh and Weissman 2012, *Nature*
  484:538. Chi site `GCTGGTGG` (*E. coli*): the RecBCD recombination hotspot (Smith laboratory).
* GC and homopolymer limits are measured on the host's own genome (RefSeq NC_012971.2 for BL21(DE3),
  NC_000964.3 for *B. subtilis* 168) by `tools/build_limits.py`: gene GC and GC in 50 bp windows
  take the central 99% of the host's coding sequences, and a base's run limit is the longest run
  found in at least three coding sequences. *E. coli*: GC 32 to 60%, windows 28 to 70%, runs A and
  T 9, G and C 7. The earlier hand-chosen defaults were stricter than natural genes.
* Direct repeats of 23 bp or more are avoided, the minimal length for efficient RecBC-dependent
  recombination in *E. coli* (Shen and Huang 1986, *Genetics* 112:441). No threshold is cited for
  inverted repeats, so that check is off by default.
* The stop-codon note in *E. coli* is Poole et al. 1995 (*EMBO Journal*): the base after the stop
  changes termination from about 7% to about 80%.

The eukaryotic profiles are in development and keep their old defaults: rare codons are the codons
used for under 10% of their amino acid's codons, and the motifs (`AATAAA`, `ATTAAA`, `ATTTA` for
human and mouse) are textbook sequences not checked against a reference here.

## Roadmap

| version | adds |
|---|---|
| **0.1** | data model, frame engine, round-trip validation, CLI ✅ |
| **0.2** | Aho-Corasick DFA + DP codon optimiser (both strands, IUPAC, junction-aware, unavoidable sites reported); cited codon table; `MAX_CAI` and seeded `WEIGHTED_SAMPLE` strategies ✅ |
| **0.2.x** | the local website; reading your own plasmid (GenBank, SnapGene, FASTA); start-codon detection; picking two sites by position, with suggested pairs ranked by what they add to the protein; C-terminus choice; tags and linkers; construct report; GC and repeat refinement with adjustable limits; batch design from a FASTA; host profiles (bacteria released, eukaryotes in development); host motifs and homopolymer limits as soft preferences; the DNA to order with extra bases outside the sites and an in-silico digest and ligation check ✅ |
| 0.3 | first PyPI release |
| 0.4 | the `HARMONIZED` strategy |
| 0.5 | Dam/Dcm methylation and buffers, validation on more real vectors, a hosted version of the website if it can run in the browser, a Korean interface |
| later | Gibson and Golden Gate assembly, multi-chain constructs, eukaryotic vectors, a vector-free mode, Addgene |

## Bundled vector

`pTEST1` is a **synthetic** 1149 bp pET-like vector used by the test suite.
Do not order it. Use your own plasmids: `read_vector_file()` from
Python, or the website's upload and `~/.p2d/vectors` folder.

## Licence

MIT, see `LICENSE`. p2d has no compiled or non-permissive dependencies.

### Third-party data

* **Codon usage** comes from [`python_codon_tables`](https://pypi.org/project/python_codon_tables/)
  (CC0, derived from Kazusa CUTG). The `ecoli_bl21` profile uses the *E. coli* K-12 table for
  BL21, which is standard practice and is recorded in the profile.
* **Plasmid files are not redistributed.** `pTEST1` is synthetic. SnapGene, Addgene and vendor
  maps carry their own terms (SnapGene's, for example, require attribution and a licence for
  commercial use), so keep them on your own computer: put them in `~/.p2d/vectors` or upload
  them in the UI, and they never leave it.
* **Fonts**: the website ships Geist and Geist Mono (Vercel, SIL Open Font License 1.1, licence text
  in `src/p2d_ui/static/fonts/OFL.txt`), so it fetches nothing from the internet.
* Depends on Biopython (Biopython License), PyYAML and Typer (MIT).
