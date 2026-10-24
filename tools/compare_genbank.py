"""Rebuild a published construct from its empty vector with p2d and compare the protein it makes.

A GenBank record of a finished expression plasmid is the ground truth: somebody cloned a gene into a
known vector and the record holds the result.  Given the empty vector, this script

1. lines the construct up with the vector and finds which piece of the vector was replaced,
2. names the two cloning sites that bound it (they stay in the construct, so they can be read off),
3. reads the fusion protein the construct really makes, from the vector's own start codon,
4. gives p2d only the vector, the two sites and the inserted protein, and
5. checks that p2d's plasmid translates to exactly the same fusion protein.

Codons differ (p2d recodes the gene); the protein, the reading frame and the junction residues
must not.
The empty vector and the construct are read from files; nothing is written.
"""

from __future__ import annotations

from dataclasses import dataclass, field

from Bio import SeqIO
from Bio.Seq import Seq

import p2d
from p2d.enzymes import list_cut_sites, plan_cut

ANCHOR = 40
MIN_COVERAGE = (
    0.8  # share of the insert the real protein must cover for the comparison to mean anything
)
JUNCTION = 8  # bases a site may lie from the point where the construct leaves the vector


def rc(s: str) -> str:
    return str(Seq(s).reverse_complement())


@dataclass
class Comparison:
    name: str
    # "match": same protein; "mismatch": p2d made a different protein; "declined": p2d refuses the
    # pair the real construct used; "not comparable": the record cannot test p2d (reason given)
    status: str
    reason: str = ""
    up: str = ""
    down: str = ""
    cterm: str = ""
    real_protein: str = ""
    p2d_protein: str = ""
    payload: str = ""
    real_bp: int = 0
    p2d_bp: int = 0
    coverage: float = 0.0
    notes: list[str] = field(default_factory=list)

    @property
    def ok(self) -> bool:
        return self.status == "match"


def _align(vector_seq: str, construct: str) -> tuple[str, int, int] | None:
    """The construct in the vector's coordinates, plus the shared prefix and suffix lengths."""
    n = len(vector_seq)
    for cand in (construct, rc(construct)):
        for offset in range(0, n - ANCHOR, max(1, (n - ANCHOR) // 7)):
            anchor = vector_seq[offset : offset + ANCHOR]
            doubled = cand + cand
            k = doubled.find(anchor)
            if k < 0:
                continue
            c = doubled[k : k + len(cand)]  # begins where the vector's own base ``offset`` is
            v = vector_seq[offset:] + vector_seq[:offset]
            pre = 0
            while pre < min(len(c), len(v)) and c[pre] == v[pre]:
                pre += 1
            suf = 0
            while suf < min(len(c), len(v)) - pre and c[-1 - suf] == v[-1 - suf]:
                suf += 1
            back = c[-offset:] + c[:-offset] if offset else c  # to the vector's own origin
            return back, (pre + offset) % n, suf
    return None


def _translate_to_stop(dna: str) -> str:
    return str(Seq(dna[: len(dna) - len(dna) % 3]).translate(to_stop=True))


def compare(name: str, vector: p2d.Vector, construct_seq: str, host: p2d.HostProfile) -> Comparison:
    vs = vector.seq.upper()
    out = Comparison(name, "not comparable")
    aligned = _align(vs, construct_seq.upper())
    if aligned is None:
        out.reason = "the construct does not contain this vector's backbone"
        return out
    c, pre_end, suf = aligned
    n = len(vs)
    if suf >= n or pre_end + suf > n + 0:
        out.reason = "the construct and the vector are the same sequence"
        return out
    start = vector.expression_start
    diverge_end = n - suf  # first vector base that is no longer shared on the 3' side
    # a restriction-cloned insert leaves its two sites right at the junctions; a site far from the
    # point where the construct leaves the vector was not used (seamless cloning leaves none)
    ups = [s for s in list_cut_sites(vector) if s.end <= pre_end and pre_end - s.end <= JUNCTION]
    downs = [
        s
        for s in list_cut_sites(vector)
        if s.start >= diverge_end and s.start - diverge_end <= JUNCTION
    ]
    if not ups or not downs:
        out.reason = "no cloning site at the junction (cloned without restriction sites)"
        return out
    up = max(ups, key=lambda s: s.end)
    down = min(downs, key=lambda s: s.start)
    out.up, out.down = up.id, down.id
    shift = len(c) - n
    real = _translate_to_stop(c[start:])
    out.real_protein, out.real_bp = real, len(c)

    pad_up = (3 - (up.end - start) % 3) % 3
    first = (up.end + pad_up - start) // 3
    down_in_c = down.start + shift
    if 3 * len(real) + start <= down_in_c:  # the real protein stops inside the insert
        modes = [p2d.CTerm.STOP]
        payload = real[first:]
    else:  # it runs on into the vector after the 3' site
        modes = [p2d.CTerm.TAG, p2d.CTerm.VECTOR_FRAME]
        payload = real[first : (down_in_c - start) // 3]
    out.payload = payload
    if not payload:
        out.reason = "nothing between the two sites is translated"
        return out
    insert_bp = (
        len(c)
        - n
        + (down.end - up.start)
        - len(vs[up.start : up.end])
        - len(vs[down.start : down.end])
    )
    out.coverage = 3 * len(payload) / max(insert_bp, 1)
    if out.coverage < MIN_COVERAGE:
        out.reason = (
            f"only {len(payload)} residues are made in frame from the vector's start codon, "
            f"out of an insert of about {insert_bp // 3} codons: the gene is not read through"
        )
        return out
    try:
        site = plan_cut(vector, up, down)
    except ValueError as exc:
        out.reason = f"p2d refuses this pair: {exc}"
        out.status = "declined"
        return out

    best = None
    for mode in modes:
        try:
            a = p2d.assemble(
                p2d.design_insert(
                    p2d.ConstructPlan.single_chain(payload),
                    vector,
                    site,
                    host,
                    strategy=p2d.CodonStrategy.MAX_CAI,
                    cterm=mode,
                )
            )
        except ValueError as exc:
            out.notes.append(f"{mode.value}: {exc}")
            continue
        got = a.protein
        if got == real and (best is None or best[0] != real) or best is None:
            best = (got, mode, len(a.seq), a)
    if best is None:
        out.status, out.reason = "declined", "p2d could not design it: " + "; ".join(out.notes)
        return out
    out.p2d_protein, out.cterm, out.p2d_bp = best[0], best[1].value, best[2]
    out.status = "match" if best[0] == real else "mismatch"
    return out


def read_construct(path: str) -> tuple[str, str]:
    rec = SeqIO.read(path, "genbank")
    return rec.id, str(rec.seq)


def first_difference(a: str, b: str) -> int | None:
    for i, (x, y) in enumerate(zip(a, b, strict=False)):
        if x != y:
            return i
    return None if len(a) == len(b) else min(len(a), len(b))


def main(argv: list[str] | None = None) -> int:
    import argparse

    from p2d.analysis import analyze_vector, confirm
    from p2d.vectorio import read_vector_file

    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("vector", help="the empty vector (GenBank, SnapGene .dna or FASTA)")
    ap.add_argument("constructs", nargs="+", help="GenBank records of finished constructs")
    ap.add_argument("--start", type=int, help="1-based start codon, if the vector does not show it")
    ap.add_argument("--host", default="ecoli_bl21")
    args = ap.parse_args(argv)

    analysis = analyze_vector(read_vector_file(args.vector))
    vector = confirm(analysis, start=None if args.start is None else args.start - 1)
    host = p2d.load_host(args.host)
    counts: dict[str, int] = {}
    print(f"{'record':12} {'result':15} {'sites':22} {'cterm':6} {'aa':>9}  note")
    for path in args.constructs:
        name, seq = read_construct(path)
        r = compare(name, vector, seq, host)
        counts[r.status] = counts.get(r.status, 0) + 1
        sites = f"{r.up} to {r.down}" if r.up else ""
        aa = f"{len(r.real_protein)}" if r.status in ("match", "mismatch") else ""
        note = r.reason
        if r.status == "mismatch":
            i = first_difference(r.real_protein, r.p2d_protein)
            note = f"first difference at residue {i + 1 if i is not None else '?'}"
        print(f"{r.name:12} {r.status:15} {sites:22} {r.cterm:6} {aa:>9}  {note}")
    print("\n" + ", ".join(f"{n} {k}" for k, n in sorted(counts.items())))
    return 1 if counts.get("mismatch") else 0


if __name__ == "__main__":
    raise SystemExit(main())
