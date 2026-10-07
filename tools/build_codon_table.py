"""Regenerate the ``codon_usage`` block of ``hosts/ecoli_bl21.yaml``.

Frequencies come from ``python_codon_tables`` (Edinburgh Genome Foundry, derived
from Kazusa CUTG) and are never typed by hand, so the provenance string always
names the exact package version that produced the numbers.

    python tools/build_codon_table.py            # print the block
    python tools/build_codon_table.py --write    # splice it into the host YAML
"""

from __future__ import annotations

import argparse
import re
from importlib.metadata import version
from pathlib import Path

import python_codon_tables as pct

TABLE_NAME = "e_coli_316407"
HOST_YAML = Path(__file__).resolve().parents[1] / "src" / "p2d" / "hosts" / "ecoli_bl21.yaml"

NOTE = (
    "K-12 table used for a BL21 (B strain) profile -- standard practice; "
    "revisit if B-strain-specific data becomes available"
)


def _flow(codons: dict[str, float]) -> str:
    ranked = sorted(codons.items(), key=lambda kv: (-kv[1], kv[0]))
    return "{" + ", ".join(f"{c}: {f:g}" for c, f in ranked) + "}"


def build_block() -> str:
    table = pct.get_codons_table(TABLE_NAME)
    source = (
        f"python_codon_tables {version('python_codon_tables')}, table '{TABLE_NAME}' "
        "(E. coli K-12, Kazusa CUTG)"
    )
    lines = [
        "codon_usage:",
        f'  source: "{source}"',
        f'  note: "{NOTE}"',
        "  verified: true",
        "  # fraction of each amino acid's codons, per codon",
        "  table:",
    ]
    lines += [f"    {aa}: {_flow(table[aa])}" for aa in sorted(k for k in table if k != "*")]
    lines.append(f"  stop: {_flow(table['*'])}")
    return "\n".join(lines) + "\n"


def splice(text: str, block: str) -> str:
    """Replace the existing top-level ``codon_usage`` mapping, leaving the rest alone."""
    pattern = re.compile(r"^codon_usage:\n(?:[ \t]+.*\n|\n)*", re.MULTILINE)
    if not pattern.search(text):
        raise SystemExit("no codon_usage block found in host YAML")
    return pattern.sub(lambda _: block + "\n", text, count=1)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--write", action="store_true", help="update the host YAML in place")
    args = parser.parse_args()

    block = build_block()
    if args.write:
        HOST_YAML.write_text(splice(HOST_YAML.read_text(), block))
        print(f"updated {HOST_YAML}")
    else:
        print(block, end="")


if __name__ == "__main__":
    main()
