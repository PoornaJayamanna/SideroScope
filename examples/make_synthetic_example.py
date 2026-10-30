"""Generate the synthetic demo files in this folder.

The sequence is RANDOMLY GENERATED for demonstration only. It is not a real organism
and the "genes" are not real genes. Replace it with a record from NCBI (for example a
bacterial contig or a gene cluster you are interested in) to explore real data.

Run from the repository root:  python examples/make_synthetic_example.py
"""

import random
from pathlib import Path

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord

random.seed(42)

HERE = Path(__file__).parent
STOPS = {"TAA", "TAG", "TGA"}
SENSE_CODONS = [a + b + c for a in "ACGT" for b in "ACGT" for c in "ACGT" if a + b + c not in STOPS]


def background(n: int, gc: float = 0.55) -> str:
    """Random bases with the given GC fraction."""
    return "".join(
        random.choice("GC") if random.random() < gc else random.choice("AT") for _ in range(n)
    )


def make_orf(n_codons: int) -> str:
    """ATG + (n_codons - 1) random sense codons + a stop codon."""
    body = "".join(random.choice(SENSE_CODONS) for _ in range(n_codons - 1))
    return "ATG" + body + random.choice(sorted(STOPS))


# (strand, number of codons including the start codon)
GENES = [("+", 220), ("+", 310), ("-", 180), ("+", 150), ("-", 260)]

sequence = background(200)
features = []
for k, (strand, n_codons) in enumerate(GENES, start=1):
    orf = make_orf(n_codons)
    segment = orf if strand == "+" else str(Seq(orf).reverse_complement())
    start0 = len(sequence)
    sequence += segment
    features.append(
        SeqFeature(
            FeatureLocation(start0, start0 + len(segment), strand=1 if strand == "+" else -1),
            type="CDS",
            qualifiers={
                "gene": [f"demo{k}"],
                "locus_tag": [f"SYNTH_{k:04d}"],
                "product": ["synthetic hypothetical protein"],
            },
        )
    )
    sequence += background(random.randint(150, 300))

record = SeqRecord(
    Seq(sequence),
    id="SYNTH_DEMO.1",
    name="SYNTH_DEMO",
    description="Synthetic demo sequence (randomly generated, not a real organism)",
    features=features,
    annotations={"molecule_type": "DNA"},
)

SeqIO.write(record, HERE / "example.fasta", "fasta")
SeqIO.write(record, HERE / "example.gb", "genbank")
print(f"Wrote {len(sequence):,} bp with {len(features)} annotated CDS to {HERE}")
