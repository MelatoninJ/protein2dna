"""Command line interface.

Deliberately thin: every command is a few lines that call the library.  If
logic starts accumulating here, it belongs in the library instead.
"""

from __future__ import annotations

from pathlib import Path

import typer
from Bio import SeqIO

from . import __version__
from .assemble import assemble
from .frame import Severity, design_insert
from .host import available_hosts, load_host
from .model import ConstructPlan
from .vector import Vector
from .vectors import available_bundled, load_bundled

app = typer.Typer(add_completion=False, help="protein-to-DNA reverse translation, vector-aware")

COLOUR = {
    Severity.ERROR: typer.colors.RED,
    Severity.WARNING: typer.colors.YELLOW,
    Severity.INFO: typer.colors.BLUE,
}


@app.command()
def version() -> None:
    typer.echo(f"p2d {__version__}")


@app.command()
def hosts() -> None:
    """List available host profiles."""
    for slug in available_hosts():
        h = load_host(slug)
        mark = "" if h.codon_table_verified else "  (unverified placeholder table)"
        typer.echo(f"{slug:<16} {h.name}{mark}")


@app.command()
def vectors() -> None:
    """List bundled vectors and their insertion sites."""
    for name in available_bundled():
        v = load_bundled(name)
        typer.echo(f"{name}  ({len(v)} bp)")
        for key, site in v.sites.items():
            typer.echo(f"    {key:<14} {site.note}")


@app.command()
def design(
    protein: str = typer.Argument(..., help="payload amino-acid sequence, or a FASTA path"),
    site: str = typer.Option("NdeI-XhoI", "--site", "-s"),
    vector_path: Path | None = typer.Option(None, "--vector", "-v", help="GenBank file"),
    host: str = typer.Option("ecoli_bl21", "--host"),
    gb: Path | None = typer.Option(None, "--gb", help="write the assembled plasmid here"),
    fasta: Path | None = typer.Option(None, "--fasta", help="write the insert here"),
) -> None:
    """Design an insert and validate it against the vector."""
    payload = _read_payload(protein)

    if vector_path:
        v = Vector.from_genbank(vector_path)
        if site not in v.sites:
            typer.secho(
                f"vector {v.name} has no registered insertion site {site!r}. "
                "v0.1 needs the site defined in code (see p2d/vectors.py); "
                "automatic enzyme-pair search lands in v0.4.",
                fg=typer.colors.RED,
            )
            raise typer.Exit(2)
    else:
        v = load_bundled()

    plan = ConstructPlan.single_chain(payload, host=host)
    a = assemble(design_insert(plan, v, v.sites[site], load_host(host)))

    typer.secho(f"\n{v.name} / {site}\n", bold=True)
    typer.echo(a.fusion_report())
    typer.echo("")
    for issue in a.issues:
        typer.secho(f"  {issue}", fg=COLOUR[issue.severity])

    typer.echo(f"\n  insert: {len(a.insert_dna)} bp\n  plasmid: {len(a.seq)} bp")
    typer.echo(f"\n{a.insert_dna}\n")

    if gb:
        SeqIO.write(a.to_record(), str(gb), "genbank")
        typer.echo(f"wrote {gb}")
    if fasta:
        fasta.write_text(f">{v.name}_{plan.chains[0].name}_insert\n{a.insert_dna}\n")
        typer.echo(f"wrote {fasta}")

    if not a.ok:
        raise typer.Exit(1)


def _read_payload(arg: str) -> str:
    p = Path(arg)
    if p.exists():
        return str(next(SeqIO.parse(str(p), "fasta")).seq)
    return arg


if __name__ == "__main__":
    app()
