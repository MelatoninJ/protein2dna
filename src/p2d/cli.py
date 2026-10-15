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
from .frame import CTerm, Severity, design_insert
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
        if h.status == "development":
            mark += "  (in development)"
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
                "The command line needs the site defined in code (see p2d/vectors.py); "
                "the website (`p2d ui`) lets you pick any two sites of your own plasmid.",
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


@app.command()
def batch(
    fasta: Path = typer.Argument(..., help="FASTA file with one protein per record"),
    out: Path = typer.Option(..., "--out", "-o", help="folder to write the results into"),
    vector_path: Path | None = typer.Option(
        None, "--vector", "-v", help="plasmid file (GenBank, SnapGene .dna or FASTA)"
    ),
    start: int | None = typer.Option(None, "--start", help="1-based position of the start codon"),
    upstream: str = typer.Option("NdeI", "--upstream", "-u", help="5' site, ex. NcoI or XhoI:1367"),
    downstream: str = typer.Option("XhoI", "--downstream", "-d", help="3' site"),
    host: str = typer.Option("ecoli_bl21", "--host"),
    cterm: str = typer.Option("stop", "--cterm", help="stop, tag or vector_frame"),
    strategy: str = typer.Option(
        "weighted_sample", "--strategy", help="weighted_sample or max_cai"
    ),
    seed: int | None = typer.Option(None, "--seed"),
    clamp_5: int = typer.Option(6, "--clamp-5", help="extra bases outside the 5' site"),
    clamp_3: int = typer.Option(6, "--clamp-3", help="extra bases outside the 3' site"),
    n_part: list[str] = typer.Option([], "--n-part", help="tag or linker before the protein"),
    c_part: list[str] = typer.Option([], "--c-part", help="tag or linker after the protein"),
) -> None:
    """Design every protein in a FASTA file into one vector and pair of sites."""
    from .analysis import analyze_vector, confirm
    from .batch import design_batch, read_proteins, write_batch
    from .enzymes import resolve_site
    from .model import resolve_parts
    from .optimize import CodonStrategy
    from .vectorio import read_vector_file

    def fail(message: str) -> typer.Exit:
        typer.secho(message, fg=typer.colors.RED, err=True)
        return typer.Exit(2)

    try:
        proteins = read_proteins(fasta.read_text())
    except OSError as exc:
        raise fail(f"cannot read {fasta}: {exc}") from exc
    if not proteins:
        raise fail(f"no protein records found in {fasta}")
    try:
        cterm_value = CTerm(cterm)
        strategy_value = CodonStrategy(strategy)
        parts = (tuple(resolve_parts(n_part)), tuple(resolve_parts(c_part)))
        host_profile = load_host(host)
        if vector_path:
            analysis = analyze_vector(read_vector_file(vector_path))
            if start is None:
                best = analysis.best_start
                if best is None or best.confidence != "high":
                    options = (
                        ", ".join(str(s.position + 1) for s in analysis.starts) or "none found"
                    )
                    raise ValueError(f"choose the start codon with --start (candidates: {options})")
                begin = best.position
            else:
                begin = start - 1
            v = confirm(analysis, start=begin)
        else:
            v = load_bundled()
        up, down = resolve_site(v, upstream), resolve_site(v, downstream)
        result = design_batch(
            proteins,
            v,
            up,
            down,
            host_profile,
            cterm=cterm_value,
            strategy=strategy_value,
            seed=seed,
            clamp_5=clamp_5,
            clamp_3=clamp_3,
            n_parts=parts[0],
            c_parts=parts[1],
        )
    except ValueError as exc:
        raise fail(str(exc)) from exc

    write_batch(result, out)
    typer.secho(f"{v.name}: {up.label} to {down.label}, host {host_profile.slug}", bold=True)
    for e in result.entries:
        if e.ok:
            typer.echo(f"  ok      {e.name}  ({len(e.protein)} aa)")
        else:
            typer.secho(f"  failed  {e.name}  {e.message}", fg=typer.colors.RED)
    typer.echo(f"\n{result.n_ok} of {len(result.entries)} designed; results are in {out}")
    if not result.ok:
        raise typer.Exit(1)


@app.command()
def ui(
    port: int = typer.Option(8765, "--port", "-p", help="local port to serve on"),
    no_browser: bool = typer.Option(False, "--no-browser", help="do not open a browser tab"),
) -> None:
    """Open the local web UI (a shell over this library; runs only on this computer)."""
    from p2d_ui import serve  # imported here so the core CLI never depends on the UI package

    serve(port=port, open_browser=not no_browser)


def _read_payload(arg: str) -> str:
    p = Path(arg)
    if p.exists():
        return str(next(SeqIO.parse(str(p), "fasta")).seq)
    return arg


if __name__ == "__main__":
    app()
