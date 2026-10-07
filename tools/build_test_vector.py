"""Generate pTEST1, a compact synthetic vector that mimics pET-28a topology.

Deliberately synthetic: it is small enough to read in a test failure message,
and it makes no claim to be a real plasmid.  Real vectors come from user
GenBank upload or the Addgene API.
"""
import random

from Bio import SeqIO
from Bio.Restriction import RestrictionBatch
from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord

ENZYMES = ["NdeI", "BamHI", "EcoRI", "XhoI", "NcoI", "HindIII", "SalI", "NotI"]

parts = []


def add(label, seq, ftype="misc_feature", note=None):
    start = sum(len(s) for _, s, _, _ in parts)
    parts.append((label, seq.upper(), ftype, note))
    return start


utr = "TAATACGACTCACTATAGGGGAATTGTGAGCGGATAACAATTCCCCTCTAGAAATAATTTTGTTTAACTTTAAGAAGGAGATATACC"
add("T7 promoter / lac operator / RBS", utr)
atg_pos = add("start codon", "ATG")
add("His6 tag", "GGCAGCAGCCATCATCATCATCATCAC")
add("thrombin site", "AGCAGCGGCCTGGTGCCGCGCGGCAGC")
nde_pos = add("NdeI", "CATATG")
bam_pos = add("BamHI", "GGATCC")
add("stuffer", "GGTACCGGCAGC")
eco_pos = add("EcoRI", "GAATTC")
xho_pos = add("XhoI", "CTCGAG")
add("His6 tag (C-term)", "CACCACCACCACCACCAC")
add("stop", "TGA")
add("T7 terminator", "CTAGCATAACCCCTTGGGGCCTCTAAACGGGTCTTGAGGGGTTTTTTG")

# ---- filler backbone, screened against every enzyme we care about ----------
rng = random.Random(20261007)
sites = [RestrictionBatch([e]) for e in ENZYMES]
site_seqs = [list(b)[0].site for b in sites]


def clean(n):
    while True:
        s = "".join(rng.choice("ACGT") for _ in range(n))
        if any(x in s for x in site_seqs):
            continue
        if any(b * 6 in s for b in "ACGT"):
            continue
        gc = (s.count("G") + s.count("C")) / len(s)
        if not 0.42 <= gc <= 0.58:
            continue
        return s


add("ori (placeholder)", clean(420))
add("KanR (placeholder)", clean(480))

seq = "".join(s for _, s, _, _ in parts)
rec = SeqRecord(Seq(seq), id="pTEST1", name="pTEST1",
                description="synthetic pET-like test vector for p2d (not a real plasmid)")
rec.annotations["molecule_type"] = "ds-DNA"
rec.annotations["topology"] = "circular"
rec.annotations["organism"] = "synthetic construct"

pos = 0
for label, s, ftype, note in parts:
    q = {"label": [label]}
    if note:
        q["note"] = [note]
    rec.features.append(SeqFeature(FeatureLocation(pos, pos + len(s)), type=ftype, qualifiers=q))
    pos += len(s)

# ---- verification ----------------------------------------------------------
print(f"length: {len(seq)} bp")
for e in ENZYMES:
    n = len(list(RestrictionBatch([e]))[0].search(Seq(seq), linear=False))
    print(f"  {e:<8} sites: {n}")

print(f"\natg_pos={atg_pos} nde={nde_pos} bam={bam_pos} eco={eco_pos} xho={xho_pos}")
for name, p in [("NdeI", nde_pos), ("BamHI", bam_pos), ("EcoRI", eco_pos), ("XhoI", xho_pos)]:
    print(f"  {name}: (start - ATG) % 3 = {(p - atg_pos) % 3}")

orf = seq[atg_pos:]
prot = Seq(orf).translate(to_stop=True)
print(f"\nempty-vector ORF: {prot}  ({len(prot)} aa)")

SeqIO.write(rec, "src/p2d/data/pTEST1.gb", "genbank")
print("\nwritten")
