"""Synthetic plasmids for tests: a hand-built SnapGene ``.dna`` writer and pTEST1 annotations."""

from __future__ import annotations

import struct

from p2d.vector import Vector


def snapgene_bytes(seq: str, features: list[dict], circular: bool = True) -> bytes:
    """Minimal ``.dna`` writer: cookie, DNA and feature packets (Biopython's format notes)."""

    def packet(kind: int, data: bytes) -> bytes:
        return struct.pack(">BI", kind, len(data)) + data

    xml = [f'<?xml version="1.0"?><Features nextValidID="{len(features)}">']
    for f in features:
        name, kind = f["name"], f["type"]
        attrs = f'name="{name}" type="{kind}"'
        if f.get("reverse"):
            attrs += ' directionality="2"'
        notes = "".join(
            f'<Q name="{k}"><V text="&lt;html&gt;&lt;body&gt;{v}&lt;/body&gt;&lt;/html&gt;"/></Q>'
            for k, v in f.get("notes", {}).items()
        )
        xml.append(f'<Feature {attrs}><Segment range="{f["start"]}-{f["end"]}"/>{notes}</Feature>')
    xml.append("</Features>")
    return (
        packet(0x09, struct.pack(">8sHHH", b"SnapGene", 1, 15, 15))
        + packet(0x00, bytes([1 if circular else 0]) + seq.encode())
        + packet(0x0A, "".join(xml).encode())
    )


def ptest1_features(vector: Vector) -> list[dict]:
    """The annotations a SnapGene map of pTEST1 would plausibly carry (1-based ranges)."""
    return [
        {
            "name": "T7 promoter",
            "type": "promoter",
            "start": 1,
            "end": 19,
            "notes": {"note": "promoter for T7 RNA polymerase"},
        },
        {"name": "RBS", "type": "RBS", "start": 73, "end": 78},
        {"name": "ATG", "type": "CDS", "start": 88, "end": 90, "notes": {"product": "start codon"}},
        {
            "name": "MCS",
            "type": "misc_feature",
            "start": 145,
            "end": 180,
            "notes": {"note": "multiple cloning site"},
        },
        {
            "name": "KanR",
            "type": "CDS",
            "start": 670,
            "end": 1149,
            "notes": {"note": "confers resistance to kanamycin"},
        },
    ]


def mirror(features: list[dict], length: int) -> list[dict]:
    """The same features as stored in a file that holds the opposite strand."""
    return [
        {
            **f,
            "start": length - f["end"] + 1,
            "end": length - f["start"] + 1,
            "reverse": not f.get("reverse"),
        }
        for f in features
    ]


def promoter_only_cassette(seed: int = 3) -> str:
    """A plasmid whose only sign of an expression start is a lac operator, an RBS and an ATG.

    No NcoI or NdeI site holds the ATG, so the start is a medium-confidence guess.
    """
    import random

    rng = random.Random(seed)
    codons = ["GCT", "AAA", "GAA", "CTG", "GGC", "GTT", "ACC", "CAG"]  # no stop among them

    def flank(n):
        return "".join(rng.choice("ACGT") for _ in range(n)).replace("ATG", "ACG")

    orf = "ATG" + "".join(rng.choice(codons) for _ in range(60)) + "TAA"
    return flank(300) + "TTGTGAGCGGATAACAA" + flank(30) + "AGGAGG" + "TATATAC" + orf + flank(300)
