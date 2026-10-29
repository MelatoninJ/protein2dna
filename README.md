# p2d

[![CI](https://github.com/MelatoninJ/protein2dna/actions/workflows/ci.yml/badge.svg)](https://github.com/MelatoninJ/protein2dna/actions/workflows/ci.yml)
[![License: MIT](https://img.shields.io/badge/license-MIT-blue.svg)](LICENSE)
![Python 3.10+](https://img.shields.io/badge/python-3.10%2B-blue.svg)

Reverse-translate a protein into DNA that clones into the vector you actually have.

Most codon-optimisation tools return a coding sequence and stop. p2d also reads your plasmid, so
it can tell you what protein comes out after you ligate the insert in: the reading frame through
both junctions, the residues the restriction sites add, a C-terminal tag that is out of frame, a
His-tag the vector already carries.

**Try it in the browser: <https://melatoninj.github.io/protein2dna/>** (nothing is uploaded; see
[The website](#the-website)).

```python
import p2d

a = p2d.run("KVFLDWINEAYQRGTRVLAEMAKRGDEFVKRLIAEGHDPFEVLKELGYSE", site="NdeI-XhoI")
print(a.fusion_report())
```

```
      1-19    vector (N-terminal)      MGSSHHHHHHSSGLVPRGS
     20-21    scar: 5' junction        HM
     22-71    insert                   KVFLDWINEAYQRGTRVLAEMAKRGDEFVKRLIAEGHDPFEVLKELGYSE
     72-73    scar: 3' junction        LE
     74-79    vector (C-terminal)      HHHHHH
  total: 79 aa
```

`HM` and `LE` are the NdeI and XhoI sites read as codons. They end up in your protein, and a tool
that ignores the vector never shows them. The N-terminal tag stays because NdeI lies downstream
of this vector's start codon.

## Contents

- [Install](#install)
- [The website](#the-website)
- [Python](#python)
- [Command line](#command-line)
- [How it works](#how-it-works)
- [Checked against published constructs](#checked-against-published-constructs)
- [Hosts](#hosts)
- [Limitations](#limitations)
- [Development](#development)
- [References](#references)
- [License](#license)

## Install

p2d is not on PyPI yet. Install it from source (Python 3.10 or newer):

```bash
git clone https://github.com/MelatoninJ/protein2dna
cd protein2dna
python -m venv .venv && source .venv/bin/activate
pip install -e .
p2d --help
```

For development use `pip install -e ".[dev]"`.

## The website

```bash
p2d ui
```

This opens a local page (loopback only, no extra dependencies). The same page is published at
<https://melatoninj.github.io/protein2dna/>, where Python runs inside your browser tab through
Pyodide, so there is no server and nothing you paste or upload leaves your computer. The first
load downloads Python and takes 10 to 20 seconds.

1. **Pick a plasmid.** Use the bundled test vector, upload GenBank (`.gb`), SnapGene (`.dna`) or
   FASTA, or paste text. A plasmid stored on the opposite strand, such as pET-28a(+), is flipped.
   Files kept in `~/.p2d/vectors` (or `P2D_VECTOR_DIR`) show up in the list.
2. **Confirm where translation starts.** p2d lists the likely start codons with their evidence.
   A plasmid with no annotation is read from its sequence (see [How it works](#how-it-works)). If
   none can be found, give the position yourself, or write it in a FASTA header:
   `>my_plasmid start=5307 mcs=5305-5380` (1-based, on the strand as written; `start` must be an
   ATG, `mcs` marks the cloning region).
3. **Choose two restriction sites.** Every site is listed by position, and an enzyme that cuts
   several times is shown as "2 of 3". A site that cannot be used says why (wrong way round,
   touching, same ends, forced by your protein). **Suggest site pairs** designs every pair in the
   window and groups them by the protein each gives, fewest added residues first.
4. **Add tags and linkers (optional).** N-terminal and C-terminal tags (His6, His8, FLAG, Myc,
   Strep-II, HA, V5) and linkers (GS, G4S, (G4S)3, EAAAK) are built onto the protein and shown as
   separate parts.
5. **Choose the C-terminus.** A stop codon, or reading on into the vector's own tag, in which case
   p2d adds the bases that put it in frame.
6. **Set the extra bases for ordering.** An enzyme cuts poorly at the very end of a fragment, so
   six bases (NEB's [general guideline][neb-ends]; some enzymes need fewer) are added outside each
   site. They are cut away and never reach the plasmid.
7. **Read the summary and take the DNA to order.** Fusion protein, plasmid size before and after,
   the piece replaced, and a check digest. The order is cut and ligated into the vector in
   silico and must give the planned plasmid exactly. The finished plasmid can be saved as
   GenBank with the vector's annotations.

**Many proteins at once.** Open *Design many proteins at once*, paste FASTA or load a file, and
press *Design all*. Every protein uses the plasmid, sites, tags and host chosen above; one that
cannot be cloned is named with the reason and the others carry on. Download the summary (CSV),
the DNA to order (FASTA), or everything as a zip with one GenBank file per construct.

## Python

```python
from Bio import SeqIO

import p2d
from p2d.analysis import analyze_vector, confirm
from p2d.enzymes import find_cut_site, list_cut_sites, plan_cut
from p2d.report import construct_report
from p2d.vectorio import read_vector_file

# 1. read the plasmid (flipped if stored backwards) and confirm the start codon
analysis = analyze_vector(read_vector_file("my_plasmid.dna"))
print(analysis.best_start.position + 1, analysis.best_start.confidence)
vector = confirm(analysis)  # pass start=... to choose your own

# 2. every restriction site, by position
for site in list_cut_sites(vector, window=(5160, 5215)):
    print(site.id, site.label)

# 3. choose two sites and design
site = plan_cut(vector, find_cut_site(vector, "BamHI:5166"), find_cut_site(vector, "XhoI:5206"))
plan = p2d.ConstructPlan.single_chain("KVFLDWINEAYQRGTRVLAEMAKRGDEFVKRLIAEGHDPFEVLKELGYSE")
design = p2d.design_insert(plan, vector, site, p2d.load_host("ecoli_bl21"), cterm=p2d.CTerm.TAG)
assembly = p2d.assemble(design)

# 4. read the result
print(assembly.ok, assembly.protein)
print(construct_report(assembly)["sizes"])
SeqIO.write(assembly.to_record(), "finished_plasmid.gb", "genbank")

# 5. the sequence to order: the insert plus extra bases outside each site
order = p2d.make_order(assembly)
print(order.dna)
print(order.ligation.message)
```

`plan_cut` refuses a pair that would not work and says why.

Designing many proteins into one vector:

```python
import p2d
from p2d.batch import design_batch, read_proteins, write_batch
from p2d.enzymes import resolve_site

vector = p2d.load_bundled()
proteins = read_proteins(open("binders.fasta").read())
up, down = resolve_site(vector, "NdeI"), resolve_site(vector, "XhoI")
result = design_batch(proteins, vector, up, down, p2d.load_host("ecoli_bl21"))

print(result.n_ok, "of", len(result.entries), "designed")
write_batch(result, "results")  # summary.csv, orders.fasta, inserts.fasta, constructs/*.gb
```

| Task | Function |
|---|---|
| One call, bundled plasmid | `p2d.run(protein, site=...)` |
| Read a plasmid | `p2d.vectorio.read_vector_file`, `read_vector_text`, `read_vector_bytes` |
| Find the start codon | `p2d.analysis.analyze_vector`, then `confirm` |
| List and choose sites | `p2d.enzymes.list_cut_sites`, `find_cut_site`, `plan_cut` |
| Design and validate | `p2d.design_insert` (`CTerm`, `CodonStrategy`), `p2d.assemble` |
| Sizes, check digest, protein mass | `p2d.report.construct_report` |
| Only the codon optimiser | `p2d.optimize.optimize_coding` |

An `assemble` result has `.ok`, `.protein`, `.seq`, `.insert_dna`, `.fusion_report()` and
`.to_record()`.

## Command line

```bash
p2d hosts                                  # list host profiles
p2d vectors                                # list bundled vectors
p2d design KVFLDWINEAYQ... --site NdeI-XhoI --gb out.gb
p2d batch binders.fasta --out results/     # every protein in a FASTA file, into one vector
p2d ui                                     # the local website
```

`p2d batch` takes the same choices as the website: `--vector`, `--start` (unless the file shows
the start codon), `--upstream` and `--downstream` (an enzyme name when it cuts once, or an id
such as `XhoI:5206`), `--host`, `--cterm`, `--n-part` and `--c-part`, `--clamp-5` and
`--clamp-3`. It writes `summary.csv`, `orders.fasta`, `inserts.fasta`, `constructs/<name>.gb` and,
if anything failed, `failed.txt`. It exits 1 when some protein failed (the rest are still
written) and 2 for a mistake in the options.

## How it works

**Frame engine.** The finished plasmid is `vector[:up.start] + insert + vector[down.end:]`, with
`insert = up_site + pad_up + coding + pad_down + down_site`. The pads are the bases that keep the
reading frame; they are chosen to avoid stop codons and new copies of either site. The finished
plasmid is then translated from the vector's own start codon, and a design that does not give the
intended protein is rejected.

**Codon optimiser.** Forbidden sites (the two enzymes' sites, host motifs) are compiled into an
Aho-Corasick automaton [1] over both strands. Codon choice is a shortest path through
(residue, automaton state), so the result contains no forbidden pattern, junctions with the vector
included, and is optimal for the per-codon score. The score is codon usage from the host's table,
as in the codon adaptation index [2], or a weighted random draw (`WEIGHTED_SAMPLE`, seeded). If
no path exists, the site cannot be avoided by synonymous changes and p2d reports it. A site your
protein forces is reported, not hidden.

**GC and repeat refinement.** A local search then swaps synonymous codons to bring overall GC, GC
in 50 bp windows and direct and inverted repeats inside limits you can change. A swap is refused
if it adds a forbidden pattern. Swapping codons lowers the codon adaptation index; the website
shows both.

**Start codon.** In order of evidence: an annotated start codon, an annotated CDS that reads into
the cloning region, then a ribosome-binding-site-like motif (the Shine-Dalgarno sequence [4]) a few
bases before an ATG. For a plasmid with no annotation, an ATG counts when an `AGGAGG`-like motif
lies 4 to 10 bases before it and the ATG is part of an NcoI or NdeI site, or lies behind a T7
promoter or a lac operator. (The best spacing in *E. coli* is studied in [6]; the window here is
deliberately wide.) A Kozak consensus [5] with a long ORF behind it is offered for vectors without
a bacterial binding site. These rules are heuristics, checked on 37 NCBI vectors: NCBI pTrc99A,
pBAD24, pET-3a/11a/30a and pColdI are found; pKK223-3, pT7-5 and pBAD18/33, whose own sequence has
no start codon, are reported as having none. The cloning sites are never guessed.

**Order check.** Extra bases outside each site are chosen so they create no new cloning site or
host motif. The whole order is then cut with Biopython's cut positions [11] and ligated into the
vector in silico; the product must equal the planned plasmid. Restriction enzyme data come from
Biopython's `Bio.Restriction`, which is built from REBASE [12].

**What else it checks.**

- Pad bases are added when a pair of sites would shift the frame, and reported.
- Sites that touch (no base between them) are offered with a warning; sites that share bases are
  refused.
- Which sites you can use: order, spacing, compatible ends, and whether a second site of the same
  enzyme would cut the backbone you keep.
- The C-terminus: one stop codon (the host's most used; two in a row on request), or a read into
  the vector's tag. In *E. coli* the base after the stop codon changes how well translation ends
  [8], so p2d reports it.
- A tag in the insert that the vector already supplies.
- Host preferences, avoided when a path exists and named in a warning when the protein forces
  them: internal Shine-Dalgarno-like sequences [7] and Chi sites [9] in bacteria, polyadenylation
  signals and AU-rich elements in human and mouse, and runs of one base longer than the host's
  own genes contain.

p2d does not make promises about expression level. See [Limitations](#limitations).

Related tools that design or optimise sequences include DNA Chisel [13], GeneDesign [14] and
JCat [15]. p2d is narrower: it starts from your vector and protein, and treats the reading frame,
the scars and the sites at the junctions as the problem to solve.

## Checked against published constructs

p2d was checked against plasmids other people built. Given the empty pET-28a(+), the two sites and
the inserted protein, it rebuilt the same fusion protein for all 9 comparable constructs on NCBI
(8 genes, 6 pairs of sites), three of them down to the exact plasmid size. 15 more records cannot
test it, and the reason for each is given. The check covers the junctions and the reading frame,
not expression, and only one vector family. Details are in [VALIDATION.md](VALIDATION.md).

## Hosts

`p2d hosts` lists the profiles. Each is a YAML file in `src/p2d/hosts/`; a new host is a new
file. The codon tables come from `python_codon_tables`, derived from the Kazusa codon usage
database [3].

| Profile | Species | Table | Status |
|---|---|---|---|
| `ecoli_bl21` (default) | *Escherichia coli* BL21(DE3) | K-12 usage | released |
| `bsubtilis` | *Bacillus subtilis* | *B. subtilis* | released |
| `hsapiens` | human | *Homo sapiens* | in development |
| `mmusculus` | mouse | *Mus musculus* | in development |
| `scerevisiae` | baker's yeast | *S. cerevisiae* | in development |
| `dmelanogaster` | fruit fly | *D. melanogaster* | in development |
| `celegans` | roundworm | *C. elegans* | in development |

The bacterial profiles are released. The eukaryotic ones work but are in development: the website
and `p2d hosts` mark them, and every design made with one carries a `host_in_development`
warning, because nothing beyond a reading-frame check against one mammalian vector has been done
for them.

**Where the bacterial values come from.**

- *E. coli* rare codons AGA, AGG, AUA, CUA, GGA, CCC and CGG are the codons the Rosetta 2 strains
  carry extra tRNA genes for (Novagen/Merck, cat. 71397).
- The internal Shine-Dalgarno-like sequence `AGGAGG` (both bacteria) is from Li, Oh and Weissman
  [7]. The Chi site `GCTGGTGG` (*E. coli*) is the RecBCD recombination hotspot [9].
- GC and homopolymer limits are measured on the host's own genome (RefSeq NC_012971.2 for
  BL21(DE3), NC_000964.3 for *B. subtilis* 168) by `tools/build_limits.py`: gene GC and GC in
  50 bp windows take the central 99% of the host's coding sequences, and a base's run limit is
  the longest run found in at least three coding sequences. For *E. coli*: GC 32 to 60%, windows
  28 to 70%, runs of A and T up to 9 and of G and C up to 7.
- Direct repeats of 23 bp or more are avoided: the minimal length for efficient RecBC-dependent
  recombination in *E. coli* is 23 to 27 bp [10]. No threshold is cited for inverted repeats, so
  that check is off by default.
- The stop-codon note for *E. coli* follows [8]: the base after the stop moves termination from
  about 7% to about 80%.

The eukaryotic profiles keep their earlier defaults: rare codons are codons used for under 10% of
their amino acid's codons, and the motifs (`AATAAA`, `ATTAAA`, `ATTTA` for human and mouse) are
textbook sequences not checked against a reference here.

## Limitations

p2d makes promises about the sequence: the reading frame, the sites in the DNA, the fusion protein
and the sizes. It makes none about how much protein you will get.

- **No experimental validation.** Codons follow a published usage table for the species. The
  *E. coli* K-12 table stands in for BL21, and the other tables are whole-genome usage, not
  tuned to a cell line, strain or expression level. That they raise expression in your hands has
  not been tested.
- **No mRNA folding or free energy**, on purpose. 5′ structure does limit translation in
  *E. coli* [16, 17], but the evidence is correlational, rests on a few reporters, and no source
  gives a free-energy cut-off that p2d could check a design against.
- **Not implemented:** Dam/Dcm methylation and double-digest buffers. p2d does not model any
  synthesis vendor's rules (length, GC or complexity limits); they differ between vendors and
  change, so check a design against the rules of whoever makes it.
- **Eukaryotic limits are defaults, not citations.** Homopolymer lengths, motifs, GC range, GC
  window and repeat length come from the project template. They never block a design, but do not
  read them as vetted.
- **Scope.** One protein chain at a time, and restriction-site cloning only (no Gibson or Golden
  Gate). Tested on bacterial vectors (pET-28a(+), pGEX-4T-1/3X/2T, pUC19, a synthetic one) and one
  mammalian vector (mOrange2-C1) for the reading frame, the pad bases and the fusion protein.
  Choose the host profile that matches where the vector is used.
- **Size limits:** 5 MB per plasmid file, 250 kb of sequence, 5000 residues per protein and 100
  proteins per batch.

Review a design before you order it.

## Development

```bash
pip install -e ".[dev]"
pytest -q && ruff check .
```

Tests that need real plasmids read them from the folder in `P2D_VECTOR_SAMPLES` and skip when it
is not set; plasmid files are not in the repository because their licences differ. See
[CONTRIBUTING.md](CONTRIBUTING.md) for adding a host or a vector, and `sync-docs.sh` for
rebuilding the browser version in `docs/` after changing the website.

Design rules the code keeps:

1. A construct is a list of chains, even when there is only one. The multi-chain strategies are
   registered in `LinkageStrategy` and raise `NotImplementedError`.
2. Hosts are data, not branches.
3. Annotation labels are not trusted. ORFs and sites are derived from the sequence; the record's
   own features are used only for reporting.

Roadmap: the first PyPI release, a `HARMONIZED` codon strategy, Dam/Dcm methylation and buffers,
validation on more real vectors, a Korean interface. Later: Gibson and Golden Gate assembly,
multi-chain constructs and eukaryotic vectors.

## References

1. Aho AV, Corasick MJ. Efficient string matching: an aid to bibliographic search. *Commun ACM*
   1975;18(6):333-340. [doi:10.1145/360825.360855](https://doi.org/10.1145/360825.360855)
2. Sharp PM, Li WH. The codon adaptation index: a measure of directional synonymous codon usage
   bias, and its potential applications. *Nucleic Acids Res* 1987;15(3):1281-1295.
   [doi:10.1093/nar/15.3.1281](https://doi.org/10.1093/nar/15.3.1281)
3. Nakamura Y, Gojobori T, Ikemura T. Codon usage tabulated from international DNA sequence
   databases: status for the year 2000. *Nucleic Acids Res* 2000;28(1):292.
   [doi:10.1093/nar/28.1.292](https://doi.org/10.1093/nar/28.1.292)
4. Shine J, Dalgarno L. The 3'-terminal sequence of *Escherichia coli* 16S ribosomal RNA:
   complementarity to nonsense triplets and ribosome binding sites. *Proc Natl Acad Sci USA*
   1974;71(4):1342-1346. [doi:10.1073/pnas.71.4.1342](https://doi.org/10.1073/pnas.71.4.1342)
5. Kozak M. An analysis of 5'-noncoding sequences from 699 vertebrate messenger RNAs.
   *Nucleic Acids Res* 1987;15(20):8125-8148.
   [doi:10.1093/nar/15.20.8125](https://doi.org/10.1093/nar/15.20.8125)
6. Chen H, Bjerknes M, Kumar R, Jay E. Determination of the optimal aligned spacing between the
   Shine-Dalgarno sequence and the translation initiation codon of *Escherichia coli* mRNAs.
   *Nucleic Acids Res* 1994;22(23):4953-4957.
   [doi:10.1093/nar/22.23.4953](https://doi.org/10.1093/nar/22.23.4953)
7. Li GW, Oh E, Weissman JS. The anti-Shine-Dalgarno sequence drives translational pausing and
   codon choice in bacteria. *Nature* 2012;484(7395):538-541.
   [doi:10.1038/nature10965](https://doi.org/10.1038/nature10965)
8. Poole ES, Brown CM, Tate WP. The identity of the base following the stop codon determines the
   efficiency of in vivo translational termination in *Escherichia coli*. *EMBO J*
   1995;14(1):151-158. [doi:10.1002/j.1460-2075.1995.tb06985.x](https://doi.org/10.1002/j.1460-2075.1995.tb06985.x)
9. Taylor AF, Amundsen SK, Smith GR. Unexpected DNA context-dependence identifies a new
   determinant of Chi recombination hotspots. *Nucleic Acids Res* 2016;44(17):8216-8228.
   [doi:10.1093/nar/gkw541](https://doi.org/10.1093/nar/gkw541)
10. Shen P, Huang HV. Homologous recombination in *Escherichia coli*: dependence on substrate
    length and homology. *Genetics* 1986;112(3):441-457.
    [doi:10.1093/genetics/112.3.441](https://doi.org/10.1093/genetics/112.3.441)
11. Cock PJA, Antao T, Chang JT, et al. Biopython: freely available Python tools for computational
    molecular biology and bioinformatics. *Bioinformatics* 2009;25(11):1422-1423.
    [doi:10.1093/bioinformatics/btp163](https://doi.org/10.1093/bioinformatics/btp163)
12. Roberts RJ, Vincze T, Posfai J, Macelis D. REBASE: a database for DNA restriction and
    modification: enzymes, genes and genomes. *Nucleic Acids Res* 2015;43(D1):D298-D299.
    [doi:10.1093/nar/gku1046](https://doi.org/10.1093/nar/gku1046)
13. Zulkower V, Rosser S. DNA Chisel, a versatile sequence optimizer. *Bioinformatics*
    2020;36(16):4508-4509. [doi:10.1093/bioinformatics/btaa558](https://doi.org/10.1093/bioinformatics/btaa558)
14. Richardson SM, Wheelan SJ, Yarrington RM, Boeke JD. GeneDesign: rapid, automated design of
    multikilobase synthetic genes. *Genome Res* 2006;16(4):550-556.
    [doi:10.1101/gr.4431306](https://doi.org/10.1101/gr.4431306)
15. Grote A, Hiller K, Scheer M, et al. JCat: a novel tool to adapt codon usage of a target gene
    to its potential expression host. *Nucleic Acids Res* 2005;33(Web Server issue):W526-W531.
    [doi:10.1093/nar/gki376](https://doi.org/10.1093/nar/gki376)
16. Kudla G, Murray AW, Tollervey D, Plotkin JB. Coding-sequence determinants of gene expression
    in *Escherichia coli*. *Science* 2009;324(5924):255-258.
    [doi:10.1126/science.1170160](https://doi.org/10.1126/science.1170160)
17. Goodman DB, Church GM, Kosuri S. Causes and effects of N-terminal codon bias in bacterial
    genes. *Science* 2013;342(6157):475-479.
    [doi:10.1126/science.1241934](https://doi.org/10.1126/science.1241934)

The bibliographic details were checked against Crossref and PubMed. Not every method here is from
a paper: the start-codon rules, the ranking of site pairs and the weights in the refinement are
this project's own choices.

## License

MIT, see [LICENSE](LICENSE). p2d has no compiled or non-permissive dependencies.

- **Codon usage** comes from [`python_codon_tables`](https://pypi.org/project/python_codon_tables/)
  (CC0, derived from Kazusa CUTG).
- **Plasmid files are not redistributed.** `pTEST1`, the bundled vector, is a synthetic 1149 bp
  pET-like plasmid for the tests; do not order it. SnapGene, Addgene and vendor maps carry their
  own terms, so keep them on your computer (`~/.p2d/vectors` or the upload box).
- **Fonts:** the website ships Geist and Geist Mono (Vercel, SIL Open Font License 1.1; text in
  `src/p2d_ui/static/fonts/OFL.txt`), so it fetches nothing from the internet.
- Depends on Biopython (Biopython License), PyYAML and Typer (MIT).

The framework was built and revised together with Claude.

[neb-ends]: https://www.neb.com/tools-and-resources/usage-guidelines/cleavage-close-to-the-end-of-dna-fragments
