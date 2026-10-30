"""Plain-language reading of the numbers shown in the app.

Everything here is a rule of thumb, written to help a student read the plots, not to
replace annotation. Each note says what was measured and how sure it can be.

NOTE: this module was rebuilt after the previous version was overwritten by mistake.
It keeps the interface the app relies on: describe(...) -> list[Insight].
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np
import pandas as pd

from .core import fmt_bp

GENOME_SCALE_BP = 100_000   # below this, genome-level comments make no sense
BP_PER_GENE = 1_000         # bacteria have roughly one gene per kb
ISLAND_DROP = 8.0           # GC points below the average that mark an AT-rich stretch
ISLAND_MIN_SPAN = 5_000     # shorter dips are too noisy to call
SKEW_MIN_STRENGTH = 0.006   # below this the cumulative skew looks like noise
LEVELS = ("info", "good", "warning")


@dataclass(frozen=True)
class Insight:
    level: str  # one of LEVELS
    title: str
    text: str

    def __post_init__(self):
        if self.level not in LEVELS:
            raise ValueError(f"level must be one of {LEVELS}, not {self.level!r}")


def gc_class(gc: float) -> str:
    """'high', 'moderate', 'low' (or 'unknown' for NaN)."""
    if gc is None or math.isnan(gc):
        return "unknown"
    if gc >= 60:
        return "high"
    if gc >= 40:
        return "moderate"
    return "low"


def _window_and_step(stats: pd.DataFrame) -> tuple[int, int]:
    window = int(stats["end"].iloc[0] - stats["start"].iloc[0] + 1)
    step = int(stats["start"].iloc[1] - stats["start"].iloc[0]) if len(stats) > 1 else window
    return window, max(step, 1)


def skew_origin_terminus(stats: pd.DataFrame, length: int) -> dict | None:
    """Estimate replication origin and terminus from the cumulative GC skew.

    The origin is where the cumulative skew is lowest and the terminus where it is
    highest. `strength` is the swing of the curve per base covered (about half the
    average skew), so a flat, noisy curve scores near zero. Returns None when there
    are too few windows to say anything.
    """
    if len(stats) < 20 or length <= 0:
        return None
    window, _ = _window_and_step(stats)
    cum = stats["cum_skew"].to_numpy()
    lo, hi = int(np.argmin(cum)), int(np.argmax(cum))
    strength = float((cum[hi] - cum[lo]) / (length / window))
    return {
        "origin": float(stats["mid"].iloc[lo]),
        "terminus": float(stats["mid"].iloc[hi]),
        "strength": strength,
    }


def low_gc_stretches(stats: pd.DataFrame, mean_gc: float, min_span: int = ISLAND_MIN_SPAN) -> list[tuple[int, int]]:
    """Runs of consecutive windows at least ISLAND_DROP points below the average GC.

    Returns (start, end) pairs in bp, longest first, for runs at least min_span long.
    """
    low = (stats["gc"] < mean_gc - ISLAND_DROP).to_numpy()
    starts = stats["start"].to_numpy()
    ends = stats["end"].to_numpy()
    runs, begin = [], None
    for i, flag in enumerate(low):
        if flag and begin is None:
            begin = i
        if begin is not None and (not flag or i == len(low) - 1):
            last = i if flag else i - 1
            span = (int(starts[begin]), int(ends[last]))
            if span[1] - span[0] + 1 >= min_span:
                runs.append(span)
            begin = None
    return sorted(runs, key=lambda r: r[1] - r[0], reverse=True)


def describe(
    length: int,
    gc: float,
    ambiguous: int,
    n_orfs: int,
    min_aa: int,
    stats: pd.DataFrame,
    agreement: tuple[int, int] | None = None,
) -> list[Insight]:
    """The short list of notes shown next to the genome map.

    length      sequence length in bp
    gc          whole-sequence GC percentage (NaN if there are no A/C/G/T bases)
    ambiguous   number of ambiguous bases (N, R, Y ...)
    n_orfs      ORFs found at the current settings
    min_aa      the minimum ORF length used, in amino acids
    stats       sliding-window table from core.windowed_stats
    agreement   (matched, total) annotated CDS recovered as ORFs, for GenBank input
    """
    notes: list[Insight] = []

    # --- GC content -------------------------------------------------------
    cls = gc_class(gc)
    if cls != "unknown":
        examples = {
            "high": "Genomes such as Pseudomonas, Streptomyces and Burkholderia sit in this range.",
            "moderate": "This is the middle range, where many soil and gut bacteria sit.",
            "low": "Genomes such as Staphylococcus, Clostridium and many Bacillus sit in this range.",
        }[cls]
        notes.append(Insight("info", f"{gc:.1f}% GC is {cls}", examples))

    # --- ambiguous bases ---------------------------------------------------
    if ambiguous > 0 and length > 0:
        share = 100 * ambiguous / length
        notes.append(
            Insight(
                "warning" if share >= 1 else "info",
                f"{ambiguous:,} ambiguous bases ({share:.2f}%)",
                "Codons that contain an ambiguous base cannot be read, so ORFs can break or "
                "show an X in the protein. GC content ignores these bases.",
            )
        )

    # --- ORF count ----------------------------------------------------------
    if length >= GENOME_SCALE_BP:
        expected = length / BP_PER_GENE
        if n_orfs > 1.3 * expected:
            notes.append(
                Insight(
                    "warning",
                    f"{n_orfs:,} ORFs is more than this sequence has genes",
                    f"A bacterial genome of {fmt_bp(length)} holds about {expected:,.0f} genes. "
                    "ORFs are only stretches without a stop codon, and GC-rich DNA has few stop "
                    "codons, so chance ORFs appear on the wrong strand and frame. "
                    f"Raise the minimum length above {min_aa} aa to thin them out, and use a "
                    "gene finder such as Prodigal for real gene calls.",
                )
            )
        else:
            notes.append(
                Insight(
                    "good",
                    f"{n_orfs:,} ORFs of {min_aa}+ aa",
                    f"That is in the range expected for {fmt_bp(length)} (about {expected:,.0f} "
                    "genes). These are still candidates, not confirmed genes.",
                )
            )
    elif n_orfs == 0:
        notes.append(
            Insight(
                "info",
                f"No ORFs of {min_aa}+ aa",
                "Try a lower minimum length. Short pieces of DNA often hold only part of a gene, "
                "and ORFs without a stop codon are not reported.",
            )
        )

    # --- annotated genes recovered ---------------------------------------------
    if agreement is not None and agreement[1] > 0:
        matched, total = agreement
        if matched == total:
            notes.append(
                Insight("good", f"All {total} annotated genes are matched by an ORF",
                        "Each annotated CDS ends at the same stop codon as a predicted ORF.")
            )
        else:
            notes.append(
                Insight(
                    "info",
                    f"{matched} of {total} annotated genes are matched by an ORF",
                    "Misses usually come from genes shorter than the minimum length, genes made "
                    "of several parts, or ORFs that start earlier than the annotated start.",
                )
            )

    # --- replication origin / terminus ---------------------------------------------
    if length >= GENOME_SCALE_BP:
        est = skew_origin_terminus(stats, length)
        if est is not None:
            if est["strength"] >= SKEW_MIN_STRENGTH:
                notes.append(
                    Insight(
                        "info",
                        f"Possible replication origin near {fmt_bp(est['origin'])}",
                        f"The cumulative GC skew bottoms out there and peaks near "
                        f"{fmt_bp(est['terminus'])}, where replication would end. This works for a "
                        "complete circular chromosome. Check it against the dnaA gene.",
                    )
                )
            else:
                notes.append(
                    Insight(
                        "info",
                        "GC skew gives no clear origin",
                        "The cumulative skew is nearly flat. That is normal for plasmids and "
                        "fragmented assemblies, so no origin is suggested.",
                    )
                )

    # --- AT-rich stretches (possible genomic islands) ---------------------------------
    if length >= 50_000 and cls != "unknown":
        runs = low_gc_stretches(stats, gc)
        if runs:
            start, end = runs[0]
            notes.append(
                Insight(
                    "info",
                    f"{len(runs)} AT-rich stretch{'es' if len(runs) > 1 else ''} found",
                    f"The longest runs from {fmt_bp(start)} to {fmt_bp(end)}, with GC at least "
                    f"{ISLAND_DROP:.0f} points below the average. Regions like this are often DNA "
                    "acquired from another organism, and can carry biosynthesis gene clusters. "
                    "They are worth a closer look.",
                )
            )

    return notes
