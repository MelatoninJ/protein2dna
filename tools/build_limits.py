"""Derive a bacterial host's GC and homopolymer limits from its own genome.

The limits steer the local search that follows the codon optimiser (see p2d/optimize/refine.py).
Instead of choosing them by hand they are measured on the host's genome, so each number can be
traced to a RefSeq record and a stated rule:

* gene GC and GC in a 50 bp window: the central 99% of what the host's coding sequences contain
  (0.5th to 99.5th percentile);
* homopolymer runs: the longest run found in at least 3 coding sequences (a single odd gene does
  not set the limit);
* repeats: not measured here. The direct-repeat limit comes from the minimal length for efficient
  RecBC-dependent recombination in E. coli, 23 to 27 bp (Shen & Huang 1986, Genetics 112:441), so
  repeats of 23 bp or more are avoided. No threshold is cited for inverted repeats, so that check is
  left off by default (200 is the largest value the check accepts).

    python tools/build_limits.py ecoli_bl21 genome.gb            # print
    python tools/build_limits.py ecoli_bl21 genome.gb --write    # update the host profile

The genome is a RefSeq GenBank record, e.g. NC_012971.2 (E. coli BL21(DE3)) or NC_000964.3
(B. subtilis 168).  Only coding sequences of 150 bp or more that start with ATG, GTG or TTG and
contain no ambiguous base are used.
"""

from __future__ import annotations

import argparse
import re
from pathlib import Path

from Bio import SeqIO

HOSTS_DIR = Path(__file__).resolve().parents[1] / "src" / "p2d" / "hosts"
MARKER = "# ---- limits measured on the host's genome (tools/build_limits.py): do not edit ----"
GENOMES = {"ecoli_bl21": "NC_012971.2", "bsubtilis": "NC_000964.3"}
LOW, HIGH = 0.5, 99.5
WINDOW = 50
STEP = 10
MIN_GENES_FOR_A_RUN = 3
DIRECT_REPEAT_MAX = 22  # repeats of 23 bp or more are avoided (Shen & Huang 1986: 23 to 27 bp)
INVERTED_OFF = 200


def percentile(values: list[float], p: float) -> float:
    xs = sorted(values)
    k = (len(xs) - 1) * p / 100
    lo = int(k)
    hi = min(lo + 1, len(xs) - 1)
    return xs[lo] + (xs[hi] - xs[lo]) * (k - lo)


def coding_sequences(record, min_len: int = 150) -> list[str]:
    out = []
    for f in record.features:
        if f.type != "CDS" or "pseudo" in f.qualifiers or "pseudogene" in f.qualifiers:
            continue
        s = str(f.extract(record.seq)).upper()
        clean = len(s) >= min_len and len(s) % 3 == 0 and set(s) <= set("ACGT")
        if clean and s[:3] in ("ATG", "GTG", "TTG"):
            out.append(s)
    return out


def _gc(s: str) -> float:
    return (s.count("G") + s.count("C")) / len(s)


def derive(cds: list[str]) -> dict:
    if len(cds) < 500:
        raise ValueError(f"only {len(cds)} coding sequences; this is not a genome record")
    gene = [_gc(s) for s in cds]
    window = [_gc(s[i : i + WINDOW]) for s in cds for i in range(0, len(s) - WINDOW + 1, STEP)]
    runs: dict[str, list[int]] = {b: [] for b in "ACGT"}
    for s in cds:
        for m in re.finditer(r"A+|C+|G+|T+", s):
            runs[m.group()[0]].append(len(m.group()))
    longest = {}
    for base, lengths in runs.items():
        longest[base] = max(
            n for n in set(lengths) if sum(x >= n for x in lengths) >= MIN_GENES_FOR_A_RUN
        )
    return {
        "n_cds": len(cds),
        "gc": [round(percentile(gene, LOW), 2), round(percentile(gene, HIGH), 2)],
        "window": [round(percentile(window, LOW), 2), round(percentile(window, HIGH), 2)],
        "homopolymer": longest,
    }


def render(slug: str, accession: str, description: str, d: dict) -> str:
    h = d["homopolymer"]
    return (
        f"{MARKER}\n"
        "limits_source:\n"
        f'  genome: "{accession}, {description}; {d["n_cds"]} coding sequences of 150 bp or more"\n'
        f'  gc: "central 99% of gene GC and of GC in {WINDOW} bp windows '
        f'({LOW}th to {HIGH}th percentile)"\n'
        f'  homopolymer: "longest run found in at least {MIN_GENES_FOR_A_RUN} coding sequences"\n'
        '  repeats: "direct: Shen & Huang 1986, Genetics 112:441 (minimal length for efficient '
        'RecBC-dependent recombination, 23 to 27 bp); inverted: no cited threshold, left off"\n'
        f"homopolymer_max: {{A: {h['A']}, T: {h['T']}, G: {h['G']}, C: {h['C']}}}\n\n"
        "gc:\n"
        f"  global: [{d['gc'][0]}, {d['gc'][1]}]\n"
        f"  window: {{size: {WINDOW}, range: [{d['window'][0]}, {d['window'][1]}]}}\n\n"
        "repeats:\n"
        f"  direct_max: {DIRECT_REPEAT_MAX}\n"
        f"  inverted_max: {INVERTED_OFF}\n"
    )


def splice(text: str, block: str) -> str:
    if MARKER not in text:
        raise SystemExit("the profile has no limits marker; see CONTRIBUTING.md")
    return text[: text.index(MARKER)] + block


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("host", choices=sorted(GENOMES))
    ap.add_argument("genome", type=Path, help="RefSeq GenBank record of the host's genome")
    ap.add_argument("--write", action="store_true")
    args = ap.parse_args(argv)

    record = SeqIO.read(str(args.genome), "genbank")
    wanted = GENOMES[args.host]
    if record.id != wanted:
        raise SystemExit(f"expected {wanted} for {args.host}, got {record.id}")
    block = render(
        args.host, record.id, record.description.rstrip("."), derive(coding_sequences(record))
    )
    path = HOSTS_DIR / f"{args.host}.yaml"
    if args.write:
        path.write_text(splice(path.read_text(), block))
        print(f"updated {path}")
    else:
        print(block, end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
