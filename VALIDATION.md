# Checking p2d against published constructs

p2d turns a protein into DNA for a chosen vector. A strong test is a plasmid someone else already
built: take their finished construct from GenBank, take the empty vector it came from, give p2d only
the vector, the two restriction sites and the inserted protein, and ask whether the plasmid p2d
makes translates to the **same fusion protein**.

`tools/compare_genbank.py` does this for any vector and any GenBank records:

```bash
python tools/compare_genbank.py pET-28a_plus.dna constructs/*.gb
```

For each record it lines the construct up with the vector, finds which piece of the vector was replaced,
names the two cloning sites that bound it (they stay in the construct), reads the protein the
construct really makes from the vector's own start codon, runs p2d, and compares the two proteins
residue by residue. Codons are expected to differ, because p2d recodes the gene.

## What was compared

Every complete circular gene construct on **pET-28a(+)** that NCBI title searches for `pET28a` and
`pET-28a` returned (24 records), compared with the empty pET-28a(+). The list is not claimed to be
exhaustive, and records that were not gene constructs (modified vectors) were left out. Every one of
the 24 contains the vector's backbone unchanged, so the comparison starts from the same sequence the
authors used.

| result | records | meaning |
|---|---|---|
| match | 9 | p2d rebuilt exactly the same fusion protein |
| not comparable | 15 | the record cannot test p2d; the reason is printed for each |

The 9 that match:

| record | construct | sites | residues | plasmid, real / p2d |
|---|---|---|---|---|
| MK847897.1 | pET28A-blaNMC-A | NcoI, BamHI | 338 | 6,419 / 6,290 bp |
| MK847898.1 | pET28A-blaBIC-1 | BmtI, NotI | 317 | 6,196 / 6,196 bp |
| MK847901.1 | pET28A-blaNDM-1 | BamHI, EcoRI | 304 | 6,320 / 6,182 bp |
| MK847902.1 | pET28A-blaVIM-1 | NcoI, BamHI | 278 | 6,193 / 6,110 bp |
| MK847903.1 | pET28A-blaVIM-2 | BamHI, SacI | 300 | 6,264 / 6,164 bp |
| MK847904.1 | pET28A-blaIMP-1 | BamHI, SacI | 280 | 6,163 / 6,104 bp |
| OK264409.1 | pET28a_ZZ_proinsulin | NcoI, NdeI | 210 | 5,948 / 5,945 bp |
| OL989372.1 | pET-28a-ZaTdT12 | NdeI, XhoI | 533 | 6,835 / 6,835 bp |
| OL989373.1 | pET-28a-ZaTdT12-R335L-K337G | NdeI, XhoI | 533 | 6,835 / 6,835 bp |

That is eight different genes (the last two differ by two residues), six different pairs of sites,
and constructs in which the vector's own tag is kept (NdeI, BamHI) or removed (NcoI). Three plasmids
come out at exactly the same size as the real one. In the ZZ_proinsulin record the authors ended the
gene with two stop codons and p2d with one (3 bp). In the other five the authors' insert carried
bases outside the coding sequence (up to 138 bp), which p2d does not add.

The 15 that are not comparable:

* **11 are not read through from the vector's start codon.** In these records the inserted gene is out
  of frame with the His-tag start codon, so the fusion protein stops after 4 to 61 residues while the
  gene is 280 to 657 codons long. Rebuilding that short stretch would prove almost nothing, so they
  are excluded by a rule (the real protein must cover at least 80% of the insert).
* **4 were cloned without a restriction site at the junction** (the gene sits directly on the start
  codon, so they were made by a seamless method). p2d clones with restriction sites only.

Five records built on **pUC19** were also tried. None can be compared: they are expression
cassettes with their own promoters, not fusions to the lacZ-alpha start.

## A pair p2d used to refuse

One real construct (MK847901.1, pET28A-blaNDM-1) was cloned with BamHI and EcoRI, which sit side
by side in pET-28a(+) with no base between them. p2d used to refuse any two sites that touch. These
authors used that pair, so p2d now offers it, with a `sites_adjacent` warning that the second
enzyme cuts right at a fragment end and the digest should be checked. Sites that share bases are
still refused. The record now matches.

## What this does and does not show

It shows that, for these real constructs, p2d's reading-frame arithmetic, its handling of the
start codon, the residues the restriction sites add, and its choice of sites reproduce what other
people actually built. The inserted protein is read from the real record, so the test is about the
**junctions and the frame**, not about the gene.

It does **not** show that the protein is made, folds, or is made well. It says nothing about codon
choice. The evidence is also narrow: one vector family (pET-28a(+)), constructs mostly from two
laboratories (six of the nine matches come from one series), and bacterial vectors only. It does not
cover mammalian vectors, other tags, or other kinds of cloning.

## Reproducing it

The empty vector and the records are not in the repository (the vector is a SnapGene file with its
own terms). Fetch a record from NCBI with
`efetch.fcgi?db=nuccore&rettype=gb&retmode=text&id=MK847897.1` and keep it as `MK847897.1.gb`. The
tests in `tests/test_genbank_comparison.py` use the files named in `CONTRIBUTING.md` when
`P2D_VECTOR_SAMPLES` is set, and otherwise run on constructs assembled by hand from the bundled test
vector, which also check that the tool reports a wrong junction as a mismatch.
