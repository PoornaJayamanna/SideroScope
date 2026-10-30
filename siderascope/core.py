"""Core sequence analysis for SideroScope.

Plain Python / NumPy / pandas / Biopython logic with no plotting or UI code, so it
can be tested on its own and reused elsewhere (for example in a notebook).
"""

from __future__ import annotations

import io
import re
import textwrap
from dataclasses import dataclass, field

import numpy as np
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq

IUPAC_DNA = frozenset("ACGTNRYSWKMBDHV")
STOP_CODONS = frozenset({"TAA", "TAG", "TGA"})
ATG_ONLY = frozenset({"ATG"})
BACTERIAL_STARTS = frozenset({"ATG", "GTG", "TTG"})
SUPPORTED_TABLES = (1, 11)  # both use the same three stop codons

ORF_COLUMNS = [
    "orf_id",
    "strand",
    "frame",
    "start",
    "end",
    "length_nt",
    "length_aa",
    "start_codon",
    "stop_codon",
    "protein",
]
CDS_COLUMNS = ["start", "end", "strand", "gene", "locus_tag", "product"]

_NON_SEQUENCE = re.compile(r"[\s\d\-\*]+")


@dataclass(eq=False)
class SequenceRecord:
    """One sequence plus, for GenBank input, its annotated CDS features."""

    id: str
    description: str
    sequence: str
    cds: pd.DataFrame = field(default_factory=lambda: pd.DataFrame(columns=CDS_COLUMNS))

    @property
    def length(self) -> int:
        return len(self.sequence)


# --------------------------------------------------------------------------- #
# Reading sequences
# --------------------------------------------------------------------------- #
def clean_sequence(text: str) -> str:
    """Upper-case a sequence, drop whitespace/digits/gaps, and convert U to T.

    Raises ValueError if anything that is not an IUPAC DNA code remains, which
    usually means the input is a protein sequence.
    """
    seq = _NON_SEQUENCE.sub("", text).upper().replace("U", "T")
    bad = sorted(set(seq) - IUPAC_DNA)
    if bad:
        raise ValueError(
            "The sequence contains characters that are not valid DNA codes "
            f"({', '.join(bad[:8])}). Is this a protein sequence?"
        )
    return seq


def _first(qualifiers: dict, key: str) -> str:
    values = qualifiers.get(key, [""])
    return values[0] if values else ""


def _cds_table(record) -> pd.DataFrame:
    """Annotated CDS features of a Biopython GenBank record as a DataFrame.

    For features made of several parts (joins), the whole span is reported.
    """
    rows = []
    for feature in record.features:
        if feature.type != "CDS":
            continue
        q = feature.qualifiers
        rows.append(
            {
                "start": int(feature.location.start) + 1,
                "end": int(feature.location.end),
                "strand": "-" if feature.location.strand == -1 else "+",
                "gene": _first(q, "gene"),
                "locus_tag": _first(q, "locus_tag"),
                "product": _first(q, "product"),
            }
        )
    return pd.DataFrame(rows, columns=CDS_COLUMNS)


def parse_sequences(text: str) -> list[SequenceRecord]:
    """Read FASTA, GenBank, or a raw pasted sequence into SequenceRecord objects."""
    stripped = text.lstrip("﻿").strip()
    if not stripped:
        raise ValueError("No sequence was provided.")

    if stripped.startswith(">"):
        fmt = "fasta"
    elif stripped.upper().startswith("LOCUS"):
        fmt = "genbank"
    else:
        return [SequenceRecord("pasted_sequence", "Pasted sequence", clean_sequence(stripped))]

    try:
        parsed = list(SeqIO.parse(io.StringIO(stripped), fmt))
    except Exception as exc:  # Biopython raises several different exception types
        raise ValueError(f"Could not read the {fmt.upper()} data: {exc}") from exc

    records = []
    for rec in parsed:
        try:
            raw = str(rec.seq)
        except Exception as exc:  # e.g. a GenBank record with no sequence
            raise ValueError(f"Record '{rec.id}' has no readable sequence.") from exc
        seq = clean_sequence(raw)
        if not seq:
            continue
        cds = _cds_table(rec) if fmt == "genbank" else pd.DataFrame(columns=CDS_COLUMNS)
        records.append(SequenceRecord(rec.id, rec.description, seq, cds))

    if not records:
        raise ValueError("No sequences were found in the input.")
    return records


# --------------------------------------------------------------------------- #
# Composition, GC content and GC skew
# --------------------------------------------------------------------------- #
def base_composition(seq: str) -> dict[str, int]:
    """Counts of A, C, G, T and 'other' (N and other ambiguity codes)."""
    counts = {base: seq.count(base) for base in "ACGT"}
    counts["other"] = len(seq) - sum(counts.values())
    return counts


def gc_content(seq: str) -> float:
    """GC percentage over unambiguous bases (A, C, G, T); NaN if there are none."""
    counts = base_composition(seq)
    unambiguous = sum(counts[b] for b in "ACGT")
    if unambiguous == 0:
        return float("nan")
    return 100.0 * (counts["G"] + counts["C"]) / unambiguous


def sequence_summary(seq: str) -> dict:
    """Headline numbers for the overview tab."""
    counts = base_composition(seq)
    return {"length": len(seq), "gc": gc_content(seq), "counts": counts, "ambiguous": counts["other"]}


def windowed_stats(seq: str, window: int = 1000, step: int = 250) -> pd.DataFrame:
    """Sliding-window GC content, GC skew and cumulative GC skew.

    Columns: start, end (1-based, inclusive), mid, gc (%), skew, cum_skew.

    * GC skew per window is (G - C) / (G + C); NaN when a window has no G or C.
    * cum_skew is the running sum of the per-window skew, scaled by step / window so
      the curve has the same shape whatever step is used. With step == window it is
      the textbook cumulative GC skew over non-overlapping windows.
    * If the sequence is shorter than the window, one window covers everything.
    """
    n = len(seq)
    if n == 0:
        raise ValueError("The sequence is empty.")
    if window < 1 or step < 1:
        raise ValueError("window and step must be at least 1.")
    window = min(window, n)

    arr = np.frombuffer(seq.encode("ascii"), dtype=np.uint8)

    def cumulative(base: str) -> np.ndarray:
        out = np.zeros(n + 1, dtype=np.int64)
        np.cumsum(arr == ord(base), dtype=np.int64, out=out[1:])
        return out

    cum_a, cum_c, cum_g, cum_t = (cumulative(b) for b in "ACGT")

    starts = np.arange(0, n - window + 1, step)
    ends = starts + window
    a = cum_a[ends] - cum_a[starts]
    c = cum_c[ends] - cum_c[starts]
    g = cum_g[ends] - cum_g[starts]
    t = cum_t[ends] - cum_t[starts]

    unambiguous = a + c + g + t
    gc_sum = g + c
    gc = np.divide(100.0 * gc_sum, unambiguous, out=np.full(len(starts), np.nan), where=unambiguous > 0)
    skew = np.divide(
        (g - c).astype(float), gc_sum, out=np.full(len(starts), np.nan), where=gc_sum > 0
    )
    cum_skew = np.nancumsum(skew) * step / window

    return pd.DataFrame(
        {
            "start": starts + 1,
            "end": ends,
            "mid": (starts + 1 + ends) / 2,
            "gc": gc,
            "skew": skew,
            "cum_skew": cum_skew,
        }
    )


# --------------------------------------------------------------------------- #
# Open reading frames
# --------------------------------------------------------------------------- #
def find_orfs(
    seq: str,
    min_aa: int = 100,
    table: int = 11,
    alt_starts: bool = False,
) -> pd.DataFrame:
    """Find open reading frames on all six frames.

    An ORF runs from a start codon (ATG; also GTG and TTG if alt_starts) to the
    first in-frame stop codon (TAA, TAG, TGA). When several start codons share one
    stop, the outermost (longest) ORF is reported. ORFs that run off the end of the
    sequence without a stop codon are not reported.

    Coordinates are 1-based, inclusive, on the forward strand, and include the stop
    codon. Frames +1..+3 start at forward-strand bases 1..3; frames -1..-3 start at
    the last three bases of the forward strand (the 3' end), reading the reverse
    complement.

    These are candidate ORFs, not gene predictions.
    """
    if table not in SUPPORTED_TABLES:
        raise ValueError(f"Unsupported genetic code table {table}; use one of {SUPPORTED_TABLES}.")
    if min_aa < 1:
        raise ValueError("min_aa must be at least 1.")

    seq = seq.upper()
    n = len(seq)
    start_codons = BACTERIAL_STARTS if alt_starts else ATG_ONLY
    rows = []

    for strand, s in (("+", seq), ("-", str(Seq(seq).reverse_complement()))):
        for offset in range(3):
            start_idx = None
            for i in range(offset, len(s) - 2, 3):
                codon = s[i : i + 3]
                if start_idx is None:
                    if codon in start_codons:
                        start_idx = i
                elif codon in STOP_CODONS:
                    length_aa = (i - start_idx) // 3
                    if length_aa >= min_aa:
                        orf_seq = s[start_idx : i + 3]
                        protein = str(Seq(orf_seq).translate(table=table, to_stop=True))
                        protein = "M" + protein[1:]  # alternative starts are read as Met
                        if strand == "+":
                            start, end = start_idx + 1, i + 3
                        else:
                            start, end = n - (i + 3) + 1, n - start_idx
                        rows.append(
                            {
                                "strand": strand,
                                "frame": (offset + 1) if strand == "+" else -(offset + 1),
                                "start": start,
                                "end": end,
                                "length_nt": i + 3 - start_idx,
                                "length_aa": length_aa,
                                "start_codon": s[start_idx : start_idx + 3],
                                "stop_codon": codon,
                                "protein": protein,
                            }
                        )
                    start_idx = None

    if not rows:
        return pd.DataFrame(columns=ORF_COLUMNS)

    orfs = pd.DataFrame(rows).sort_values(["start", "end"], kind="stable").reset_index(drop=True)
    orfs.insert(0, "orf_id", [f"ORF_{k:04d}" for k in range(1, len(orfs) + 1)])
    return orfs[ORF_COLUMNS]


def orfs_to_fasta(orfs: pd.DataFrame, record_id: str = "sequence") -> str:
    """Protein translations of the ORFs as FASTA text (60 characters per line)."""
    lines = []
    for r in orfs.itertuples(index=False):
        lines.append(f">{record_id}|{r.orf_id}|strand={r.strand}|{r.start}-{r.end}|{r.length_aa}aa")
        lines.extend(textwrap.wrap(r.protein, 60))
    return "\n".join(lines) + ("\n" if lines else "")


# --------------------------------------------------------------------------- #
# Helpers used by the interface
# --------------------------------------------------------------------------- #
def fmt_bp(n: float) -> str:
    """Human-readable length: 4809 -> '4.8 kb', 6827737 -> '6.83 Mb'."""
    n = float(n)
    if n >= 1_000_000:
        return f"{n / 1e6:.2f} Mb"
    if n >= 1_000:
        return f"{n / 1e3:.1f} kb"
    return f"{int(n)} bp"


def recommended_settings(length: int) -> dict:
    """Starting settings that suit the size of the sequence.

    Long sequences get larger GC windows (so plots stay readable) and a higher minimum
    ORF length (so chance ORFs do not swamp the real genes).
    """
    if length >= 1_000_000:
        window = max(1000, int(round(length / 500, -3)))
        return {"scale": "genome", "label": "genome-scale sequence", "window": window, "min_aa": 300}
    if length >= 100_000:
        return {"scale": "contig", "label": "large contig", "window": 1000, "min_aa": 150}
    window = int(min(1000, max(20, length // 20)))
    return {"scale": "short", "label": "short sequence", "window": window, "min_aa": 100}


def window_step(length: int, window: int) -> int:
    """Step between GC windows: a quarter window, but never more than ~4000 windows."""
    return max(1, window // 4, length // 4000)


def annotation_agreement(cds: pd.DataFrame, orfs: pd.DataFrame) -> tuple[int, int]:
    """How many annotated CDS are matched by a predicted ORF.

    A CDS counts as matched when an ORF on the same strand ends at the same stop codon
    and covers the whole CDS (the ORF may start earlier, see find_orfs). Returns
    (matched, total).
    """
    total = len(cds)
    if total == 0 or orfs.empty:
        return 0, total
    # earliest ORF start per stop position (+ strand), latest ORF end per stop position (- strand)
    plus = orfs[orfs["strand"] == "+"].groupby("end")["start"].min().to_dict()
    minus = orfs[orfs["strand"] == "-"].groupby("start")["end"].max().to_dict()
    matched = 0
    for row in cds.itertuples():
        if row.strand == "+":
            first = plus.get(row.end)
            matched += first is not None and first <= row.start
        else:
            last = minus.get(row.start)
            matched += last is not None and last >= row.end
    return int(matched), total
