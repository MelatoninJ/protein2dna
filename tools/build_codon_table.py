"""Generate the ``codon_usage`` block of a host profile from ``python_codon_tables``.

Frequencies come from ``python_codon_tables`` (Edinburgh Genome Foundry, derived from Kazusa
CUTG) and are never typed by hand, so the provenance string always names the exact package
version that produced the numbers.

    python tools/build_codon_table.py                    # print the E. coli block
    python tools/build_codon_table.py --write            # splice it into ecoli_bl21.yaml
    python tools/build_codon_table.py --host hsapiens    # print a whole profile
    python tools/build_codon_table.py --host hsapiens --write   # create or refresh the file
    python tools/build_codon_table.py --all --write      # every host below

A profile that already exists keeps everything outside its ``codon_usage`` block, so edits to
its motifs and limits survive a refresh.  A new profile starts from the template at the bottom.
"""

from __future__ import annotations

import argparse
import re
from importlib.metadata import version
from pathlib import Path

import python_codon_tables as pct

HOSTS_DIR = Path(__file__).resolve().parents[1] / "src" / "p2d" / "hosts"

ECOLI_NOTE = (
    "K-12 table used for a BL21 (B strain) profile (standard practice); "
    "revisit if B-strain-specific data becomes available"
)
GENERIC_NOTE = "Whole-genome usage for the species, not for a cell line, strain or expression level"

# slug -> (display name, python_codon_tables table, taxid, organism text for the source line)
HOSTS = {
    "ecoli_bl21": ("Escherichia coli BL21(DE3)", "e_coli_316407", 469008, "E. coli K-12"),
    "hsapiens": ("Homo sapiens (human)", "h_sapiens_9606", 9606, "Homo sapiens"),
    "mmusculus": ("Mus musculus (mouse)", "m_musculus_10090", 10090, "Mus musculus"),
    "scerevisiae": (
        "Saccharomyces cerevisiae (baker's yeast)",
        "s_cerevisiae_4932",
        4932,
        "Saccharomyces cerevisiae",
    ),
    "bsubtilis": ("Bacillus subtilis", "b_subtilis_1423", 1423, "Bacillus subtilis"),
    "dmelanogaster": (
        "Drosophila melanogaster (fruit fly)",
        "d_melanogaster_7227",
        7227,
        "Drosophila melanogaster",
    ),
    "celegans": ("Caenorhabditis elegans", "c_elegans_6239", 6239, "Caenorhabditis elegans"),
}

# Textbook sequences of mammalian mRNAs, listed as soft preferences.  They are NOT checked against
# a reference here (see the comment written into each profile).
MAMMALIAN_MOTIFS = [
    ("AATAAA", "polyadenylation signal (premature cleavage inside the coding region)"),
    ("ATTAAA", "common variant of the polyadenylation signal"),
    ("ATTTA", "AU-rich element pentamer (mRNA instability)"),
]
MOTIFS = {"hsapiens": MAMMALIAN_MOTIFS, "mmusculus": MAMMALIAN_MOTIFS}

# Bacteria are released; every eukaryotic profile is still in development, because nothing about
# the vector-side workflow for them has been validated beyond one mammalian vector
STATUS = {
    "ecoli_bl21": "released",
    "bsubtilis": "released",
    "hsapiens": "development",
    "mmusculus": "development",
    "scerevisiae": "development",
    "dmelanogaster": "development",
    "celegans": "development",
}

RARE_BELOW = 0.10  # a codon used for less than this share of its amino acid's codons


def _flow(codons: dict[str, float]) -> str:
    ranked = sorted(codons.items(), key=lambda kv: (-kv[1], kv[0]))
    return "{" + ", ".join(f"{c}: {f:g}" for c, f in ranked) + "}"


def build_block(slug: str) -> str:
    _, table_name, _, organism = HOSTS[slug]
    table = pct.get_codons_table(table_name)
    source = (
        f"python_codon_tables {version('python_codon_tables')}, table '{table_name}' "
        f"({organism}, Kazusa CUTG)"
    )
    note = ECOLI_NOTE if slug == "ecoli_bl21" else GENERIC_NOTE
    lines = [
        "codon_usage:",
        f'  source: "{source}"',
        f'  note: "{note}"',
        "  verified: true",
        "  # fraction of each amino acid's codons, per codon",
        "  table:",
    ]
    lines += [f"    {aa}: {_flow(table[aa])}" for aa in sorted(k for k in table if k != "*")]
    lines.append(f"  stop: {_flow(table['*'])}")
    return "\n".join(lines) + "\n"


def rare_codons(slug: str) -> list[str]:
    """Codons used for under 10% of their amino acid's codons (only where there is a choice)."""
    table = pct.get_codons_table(HOSTS[slug][1])
    out = [
        c
        for aa, codons in table.items()
        if aa != "*" and len(codons) > 1
        for c, f in codons.items()
        if f < RARE_BELOW
    ]
    return sorted(out)


def render_new(slug: str) -> str:
    name, _, taxid, _ = HOSTS[slug]
    motifs = MOTIFS.get(slug, [])
    motif_lines = (
        "\n".join(f"  - {{pattern: {p}, reason: {why}}}" for p, why in motifs) if motifs else None
    )
    return (
        f"name: {name}\n"
        f"slug: {slug}\n"
        f"taxid: {taxid}\n"
        f"status: {STATUS[slug]}  # released, or development (shown as such and warned about)\n\n"
        "# The codon_usage block is generated; do not edit by hand.\n"
        f"# Regenerate with: python tools/build_codon_table.py --host {slug} --write\n"
        f"{build_block(slug)}\n"
        f"# Derived from the table above: codons used for under {RARE_BELOW:.0%} of their amino "
        "acid's codons.\n"
        "# This is a statement about the table, not a claim about tRNA supply.\n"
        f"rare_codons: [{', '.join(rare_codons(slug))}]\n\n"
        "# Only the codon table above is taken from a cited source. Everything below is a project\n"
        "# default that has not been checked against a reference. The motifs and homopolymer\n"
        "# limits are soft preferences: the optimiser avoids them when it can and names any that\n"
        "# stay. The GC and repeat limits steer a local search and can be changed per design.\n"
        + (
            f"forbidden_motifs:\n{motif_lines}\n"
            if motif_lines
            else "forbidden_motifs: []  # none listed for this host\n"
        )
        + "\nhomopolymer_max: {A: 8, T: 8, G: 5, C: 5}\n\n"
        "gc:\n"
        "  global: [0.35, 0.65]\n"
        "  window: {size: 50, range: [0.25, 0.75]}\n\n"
        "repeats:\n"
        "  direct_max: 20\n"
        "  inverted_max: 20\n\n"
        "required_elements: []\n"
    )


def splice(text: str, block: str) -> str:
    """Replace the existing top-level ``codon_usage`` mapping, leaving the rest alone."""
    pattern = re.compile(r"^codon_usage:\n(?:[ \t]+.*\n|\n)*", re.MULTILINE)
    if not pattern.search(text):
        raise SystemExit("no codon_usage block found in host YAML")
    return pattern.sub(lambda _: block + "\n", text, count=1)


def write_host(slug: str) -> Path:
    path = HOSTS_DIR / f"{slug}.yaml"
    if path.exists():
        path.write_text(splice(path.read_text(), build_block(slug)))
    else:
        path.write_text(render_new(slug))
    return path


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--host", choices=sorted(HOSTS), default="ecoli_bl21")
    parser.add_argument("--all", action="store_true", help="every host in this file")
    parser.add_argument("--write", action="store_true", help="write the host YAML in place")
    args = parser.parse_args()

    slugs = sorted(HOSTS) if args.all else [args.host]
    for slug in slugs:
        if args.write:
            print(f"updated {write_host(slug)}")
        elif (HOSTS_DIR / f"{slug}.yaml").exists() or slug == "ecoli_bl21":
            print(build_block(slug), end="")
        else:
            print(render_new(slug), end="")


if __name__ == "__main__":
    main()
