# Changelog

All notable changes are listed here. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/) and the project uses
[semantic versioning](https://semver.org/).

## [Unreleased]

### Added
- Plasmids with no annotation: the start of lacZ-alpha (MTMITPS, in the pUC and pBluescript
  families) is recognised from the sequence, so an unannotated pUC19 is oriented and its start
  codon offered. It was checked against a hand-designed GFP into pUC19 cloning.

## [0.2.0] - 2026-10-14

First tagged release.

### Added
- Codon optimiser: an Aho-Corasick automaton over forbidden DNA patterns (both strands, IUPAC)
  with a dynamic program over (residue, automaton state), so no forbidden site appears, including
  across the junctions with the vector. `MAX_CAI` and seeded `WEIGHTED_SAMPLE` strategies.
- Frame engine: the 5' and 3' pad bases, the vector's real start codon, and a `wrong_start` gate
  that translates the finished plasmid and checks it gives the protein reported.
- Reading your own plasmid (GenBank, SnapGene `.dna`, FASTA or pasted text), flipping a plasmid
  stored on the opposite strand, and start-codon candidates with the evidence for each.
- Choosing two restriction sites by position, with the reason when a site is unusable, plus
  plasmid sizes, the replaced piece and a check digest.
- **Site-pair suggestions**: every workable pair designed and ranked by what it adds to the
  protein, grouped by the protein each gives, with a filter for keeping or removing the vector's
  N-terminal tag.
- **Tags and linkers**: N- and C-terminal tags (His6, His8, FLAG, Myc, Strep-II, HA, V5) and
  linkers built onto the protein; the vector's own C-terminal tag is detected for the same set.
- C-terminus choice: stop codon, read into the vector's tag, or keep the vector's frame.
- Host motifs and homopolymer limits as soft preferences.
- **GC and repeat refinement**: a local search after the optimiser that swaps synonymous codons
  to bring overall GC, windowed GC and direct and inverted repeats inside limits you can change;
  it never adds a cloning site or host motif, and reports what the protein forces.
- **DNA to order**: extra bases outside each site, chosen to add no site or motif, and an in-silico
  digest and ligation that must give the planned plasmid exactly.
- Local web interface (`p2d ui`) and a command line.

### Known limits
- One protein chain, restriction-site cloning only, bacterial vectors tested.
- The host profile's rare codons, motifs, homopolymer, GC and repeat limits are uncited defaults.
- 5' mRNA folding, Dam/Dcm methylation, batch design and other hosts are not implemented.

[0.2.0]: https://github.com/MelatoninJ/protein2dna/releases/tag/v0.2.0
