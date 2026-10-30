"""Plotly figures for SideroScope.

Base colours are from the Okabe-Ito palette, which stays distinguishable for people
with common forms of colour-vision deficiency.
"""

from __future__ import annotations

import math

import numpy as np
import pandas as pd
import plotly.graph_objects as go
from plotly.subplots import make_subplots

BASE_ORDER = ["A", "C", "G", "T", "other"]
# Light tints behind the letters in the sequence view (dark text stays readable).
BASE_FILL = {"A": "#8FD3B6", "C": "#8DBDE6", "G": "#F5CD7A", "T": "#F0A27C", "other": "#D4D4D4"}
# Stronger colours for bars and lines.
BASE_STRONG = {"A": "#009E73", "C": "#0072B2", "G": "#E69F00", "T": "#D55E00", "other": "#8C8C8C"}
STRAND_COLORS = {"+": "#0072B2", "-": "#CC79A7"}
BRAND_TEAL = "#0F9D8F"
NEUTRAL = "#9AA5B1"
FRAME_LABELS = ["+1", "+2", "+3", "-1", "-2", "-3"]
ANNOTATION_LABEL = "Annotated CDS"
MAX_ORF_BARS = 6000  # drawing more than this makes the browser sluggish


def _fmt_tick(pos: float, seq_len: int) -> str:
    if pos == 0:
        return "0"
    if seq_len >= 1_000_000:
        return f"{pos / 1e6:g} Mb"
    if seq_len >= 1_000:
        return f"{pos / 1e3:g} kb"
    return f"{int(pos)} bp"


def _nice_step(seq_len: int, target: int = 10) -> int:
    """A round tick spacing (1, 2, 5 x 10^k) that gives about `target` ticks."""
    raw = seq_len / target
    magnitude = 10 ** math.floor(math.log10(max(raw, 1)))
    for m in (1, 2, 5, 10):
        if raw <= m * magnitude:
            return max(1, int(m * magnitude))
    return max(1, int(10 * magnitude))


def base_view_figure(seq: str, start: int = 1, length: int = 600, line_width: int = 60) -> go.Figure:
    """Coloured, wrapped view of one stretch of the sequence (start is 1-based)."""
    n = len(seq)
    start = max(1, min(start, n))
    end = min(n, start + length - 1)
    window = seq[start - 1 : end]
    rows = -(-len(window) // line_width)  # ceiling division

    code = {"A": 0, "C": 1, "G": 2, "T": 3}
    z, text, positions = [], [], []
    for r in range(rows):
        chunk = window[r * line_width : (r + 1) * line_width]
        pad = line_width - len(chunk)
        z.append([code.get(ch, 4) for ch in chunk] + [None] * pad)
        text.append(list(chunk) + [""] * pad)
        first = start + r * line_width
        positions.append([first + c for c in range(line_width)])

    colors = [BASE_FILL[b] for b in BASE_ORDER]
    colorscale = []
    for i, color in enumerate(colors):
        colorscale += [[i / len(colors), color], [(i + 1) / len(colors), color]]

    fig = go.Figure(
        go.Heatmap(
            z=z,
            text=text,
            customdata=positions,
            x=list(range(line_width)),
            y=list(range(rows)),
            zmin=-0.5,
            zmax=len(colors) - 0.5,
            colorscale=colorscale,
            showscale=False,
            xgap=2,
            ygap=2,
            texttemplate="%{text}",
            textfont=dict(family="JetBrainsMono, monospace", size=13, color="#1a1a1a"),
            hovertemplate="Position %{customdata:,}<br>Base %{text}<extra></extra>",
        )
    )
    fig.update_layout(
        height=min(900, 90 + 30 * rows),
        margin=dict(l=70, r=10, t=10, b=10),
        xaxis=dict(visible=False),
        yaxis=dict(
            autorange="reversed",
            tickmode="array",
            tickvals=list(range(rows)),
            ticktext=[f"{start + r * line_width:,}" for r in range(rows)],
            showgrid=False,
            zeroline=False,
        ),
    )
    return fig


def gc_figure(stats: pd.DataFrame, mean_gc: float | None = None) -> go.Figure:
    """GC content, GC skew, and cumulative GC skew along the sequence."""
    fig = make_subplots(rows=3, cols=1, shared_xaxes=True, vertical_spacing=0.07)
    x = stats["mid"]

    fig.add_trace(
        go.Scatter(
            x=x,
            y=stats["gc"],
            mode="lines",
            line=dict(color=BRAND_TEAL, width=2),
            name="GC (%)",
            hovertemplate="Position %{x:,.0f}<br>GC %{y:.1f}%<extra></extra>",
        ),
        row=1,
        col=1,
    )
    if mean_gc is not None and not np.isnan(mean_gc):
        fig.add_hline(y=mean_gc, line_dash="dash", line_color="#888888", row=1, col=1)

    skew = stats["skew"]
    fig.add_trace(
        go.Bar(
            x=x,
            y=skew,
            marker_color=np.where(skew >= 0, BASE_STRONG["G"], BASE_STRONG["C"]),
            name="GC skew",
            hovertemplate="Position %{x:,.0f}<br>Skew %{y:.3f}<extra></extra>",
        ),
        row=2,
        col=1,
    )
    fig.add_trace(
        go.Scatter(
            x=x,
            y=stats["cum_skew"],
            mode="lines",
            line=dict(color=BASE_STRONG["T"], width=2),
            name="Cumulative skew",
            hovertemplate="Position %{x:,.0f}<br>Cumulative skew %{y:.2f}<extra></extra>",
        ),
        row=3,
        col=1,
    )

    fig.update_yaxes(title_text="GC (%)", row=1, col=1)
    fig.update_yaxes(title_text="GC skew", row=2, col=1)
    fig.update_yaxes(title_text="Cumulative skew", row=3, col=1)
    fig.update_xaxes(title_text="Position (bp)", row=3, col=1)
    fig.update_layout(height=620, showlegend=False, margin=dict(l=70, r=20, t=20, b=50), bargap=0)
    return fig


def orf_map_figure(
    orfs: pd.DataFrame,
    seq_len: int,
    cds: pd.DataFrame | None = None,
    min_aa: int | None = None,
    max_bars: int = MAX_ORF_BARS,
) -> go.Figure:
    """Six-frame ORF map; annotated CDS (GenBank input) are drawn as outlines underneath.

    If there are more than max_bars ORFs, only the longest ones are drawn.
    """
    has_cds = cds is not None and len(cds) > 0
    categories = FRAME_LABELS + ([ANNOTATION_LABEL] if has_cds else [])

    total = len(orfs)
    capped = total > max_bars
    if capped:
        orfs = orfs.nlargest(max_bars, "length_aa")

    fig = go.Figure()
    # Invisible trace so every frame row is drawn even when it has no ORFs.
    fig.add_trace(
        go.Scatter(
            x=[0] * len(categories),
            y=categories,
            mode="markers",
            marker=dict(opacity=0),
            hoverinfo="skip",
            showlegend=False,
        )
    )

    for strand in ("+", "-"):
        sub = orfs[orfs["strand"] == strand]
        if sub.empty:
            continue
        fig.add_trace(
            go.Bar(
                orientation="h",
                y=[f"{int(f):+d}" for f in sub["frame"]],
                x=sub["length_nt"],
                base=sub["start"] - 1,
                marker_color=STRAND_COLORS[strand],
                name=f"ORFs, {strand} strand",
                customdata=np.column_stack([sub["orf_id"], sub["start"], sub["end"], sub["length_aa"]]),
                hovertemplate="%{customdata[0]}<br>%{customdata[1]}-%{customdata[2]} bp<br>%{customdata[3]} aa<extra></extra>",
            )
        )

    if has_cds:
        for strand in ("+", "-"):
            sub = cds[cds["strand"] == strand]
            if sub.empty:
                continue
            label = sub["gene"].where(sub["gene"] != "", sub["locus_tag"])
            fig.add_trace(
                go.Bar(
                    orientation="h",
                    y=[ANNOTATION_LABEL] * len(sub),
                    x=sub["end"] - sub["start"] + 1,
                    base=sub["start"] - 1,
                    marker=dict(color="rgba(0,0,0,0)", line=dict(color=STRAND_COLORS[strand], width=2)),
                    name=f"Annotated CDS, {strand} strand",
                    customdata=np.column_stack([label, sub["product"], sub["start"], sub["end"]]),
                    hovertemplate="%{customdata[0]}<br>%{customdata[1]}<br>%{customdata[2]}-%{customdata[3]} bp<extra></extra>",
                )
            )

    if total == 0:
        note = "No ORFs found" + (f" at {min_aa} aa or longer" if min_aa else "")
        fig.add_annotation(text=note, x=0.5, y=0.5, xref="paper", yref="paper", showarrow=False)
    elif capped:
        fig.add_annotation(
            text=f"Drawing the longest {max_bars:,} of {total:,} ORFs. Raise the minimum length to see all of them.",
            x=0, y=-0.2, xref="paper", yref="paper", showarrow=False, xanchor="left", font=dict(size=11),
        )

    fig.update_layout(
        barmode="overlay",
        bargap=0.35,
        height=130 + 38 * len(categories),
        margin=dict(l=100, r=20, t=40, b=60),
        xaxis=dict(title="Position (bp)", range=[0, seq_len]),
        # type="category" is essential: labels such as "+1" and "-2" would otherwise be
        # auto-detected as numbers and the axis would show -2, 0, 2, 4.
        yaxis=dict(type="category", categoryorder="array", categoryarray=categories, autorange="reversed"),
        legend=dict(orientation="h", y=1.15, x=0),
        hovermode="closest",
    )
    return fig


def genome_ring_figure(
    stats: pd.DataFrame,
    orfs: pd.DataFrame,
    seq_len: int,
    gc: float,
    record_id: str = "",
    max_orfs: int = MAX_ORF_BARS,
) -> go.Figure:
    """Circular genome map. From the outside in: ORFs on the + strand, ORFs on the -
    strand, GC content, and GC skew (orange where G exceeds C, blue where C exceeds G).

    The sequence is drawn clockwise from the top. If there are more than max_orfs ORFs,
    only the longest ones are drawn.
    """
    deg = 360.0 / seq_len  # degrees per base
    fig = go.Figure()

    # --- ORF rings (outer two)
    shown = orfs.nlargest(max_orfs, "length_aa") if len(orfs) > max_orfs else orfs
    for strand, base, label in (("+", 8.7, "ORFs, + strand"), ("-", 7.7, "ORFs, - strand")):
        sub = shown[shown["strand"] == strand]
        if sub.empty:
            continue
        length_nt = sub["length_nt"].to_numpy(dtype=float)
        fig.add_trace(
            go.Barpolar(
                r=np.full(len(sub), 0.9),
                base=base,
                theta=((sub["start"].to_numpy(dtype=float) - 1) + length_nt / 2) * deg,
                width=np.maximum(length_nt * deg, 0.12),  # keep tiny ORFs visible
                marker=dict(color=STRAND_COLORS[strand], line=dict(width=0)),
                name=label,
                customdata=np.column_stack([sub["orf_id"], sub["start"], sub["end"], sub["length_aa"]]),
                hovertemplate="%{customdata[0]}<br>%{customdata[1]}-%{customdata[2]} bp<br>%{customdata[3]} aa<extra></extra>",
            )
        )

    # --- window rings (GC content and skew)
    if len(stats) > 1:
        step = float(stats["start"].iloc[1] - stats["start"].iloc[0])
    else:
        step = float(stats["end"].iloc[0] - stats["start"].iloc[0] + 1)
    bar_width = step * deg

    gc_windows = stats.dropna(subset=["gc"])
    if len(gc_windows) > 0:
        lo, hi = float(gc_windows["gc"].min()), float(gc_windows["gc"].max())
        span = max(hi - lo, 1e-9)
        for above, color, label in (
            (True, BRAND_TEAL, "GC above average"),
            (False, NEUTRAL, "GC below average"),
        ):
            sub = gc_windows[(gc_windows["gc"] >= gc) == above] if not np.isnan(gc) else gc_windows.iloc[0:0]
            if sub.empty:
                continue
            fig.add_trace(
                go.Barpolar(
                    r=0.15 + 1.55 * (sub["gc"].to_numpy(dtype=float) - lo) / span,
                    base=5.6,
                    theta=(sub["mid"].to_numpy(dtype=float) - 0.5) * deg,
                    width=bar_width,
                    marker=dict(color=color, line=dict(width=0)),
                    name=label,
                    customdata=np.column_stack([sub["mid"], sub["gc"]]),
                    hovertemplate="Position %{customdata[0]:,.0f}<br>GC %{customdata[1]:.1f}%<extra></extra>",
                )
            )

    skew_windows = stats.dropna(subset=["skew"])
    if len(skew_windows) > 0:
        amplitude = float(skew_windows["skew"].abs().max()) or 1.0
        for positive, color, label in (
            (True, BASE_STRONG["G"], "GC skew, G > C"),
            (False, BASE_STRONG["C"], "GC skew, C > G"),
        ):
            sub = skew_windows[(skew_windows["skew"] >= 0) == positive]
            if sub.empty:
                continue
            fig.add_trace(
                go.Barpolar(
                    r=sub["skew"].to_numpy(dtype=float) / amplitude * 0.85,
                    base=3.95,
                    theta=(sub["mid"].to_numpy(dtype=float) - 0.5) * deg,
                    width=bar_width,
                    marker=dict(color=color, line=dict(width=0)),
                    name=label,
                    customdata=np.column_stack([sub["mid"], sub["skew"]]),
                    hovertemplate="Position %{customdata[0]:,.0f}<br>Skew %{customdata[1]:.3f}<extra></extra>",
                )
            )
        fig.add_trace(
            go.Scatterpolar(
                theta=np.linspace(0, 360, 181),
                r=np.full(181, 3.95),
                mode="lines",
                line=dict(color="rgba(128,128,128,0.55)", width=1, dash="dot"),
                hoverinfo="skip",
                showlegend=False,
            )
        )

    # --- axis ticks in bp / kb / Mb
    step_bp = _nice_step(seq_len)
    ticks = list(range(0, seq_len, step_bp))
    center = f"<b>{record_id}</b><br>{seq_len:,} bp" + ("" if np.isnan(gc) else f"<br>GC {gc:.1f}%")

    fig.update_layout(
        polar=dict(
            hole=0.28,
            bgcolor="rgba(0,0,0,0)",
            radialaxis=dict(range=[0, 10], visible=False),
            angularaxis=dict(
                rotation=90,
                direction="clockwise",
                tickmode="array",
                tickvals=[t * deg for t in ticks],
                ticktext=[_fmt_tick(t, seq_len) for t in ticks],
                tickfont=dict(size=11),
                gridcolor="rgba(128,128,128,0.18)",
                linecolor="rgba(128,128,128,0.35)",
            ),
        ),
        legend=dict(orientation="h", y=-0.06, x=0.5, xanchor="center", font=dict(size=11)),
        height=640,
        margin=dict(l=40, r=40, t=30, b=70),
        annotations=[
            dict(text=center, x=0.5, y=0.5, xref="paper", yref="paper", showarrow=False, font=dict(size=14))
        ],
    )
    return fig


def orf_length_histogram(orfs: pd.DataFrame) -> go.Figure:
    """Distribution of ORF lengths in amino acids."""
    fig = go.Figure(
        go.Histogram(
            x=orfs["length_aa"],
            nbinsx=40,
            marker=dict(color=BRAND_TEAL, line=dict(width=0)),
            hovertemplate="%{x} aa: %{y:,} ORFs<extra></extra>",
        )
    )
    fig.update_layout(
        height=300,
        margin=dict(l=60, r=10, t=10, b=50),
        xaxis_title="ORF length (aa)",
        yaxis_title="Number of ORFs",
        bargap=0.05,
    )
    return fig


def composition_figure(counts: dict[str, int], gc: float | None = None) -> go.Figure:
    """Donut chart of A, C, G, T (and other, if present) with the GC content in the middle."""
    keys = [k for k in BASE_ORDER if k in counts and (counts[k] > 0 or k != "other")]
    fig = go.Figure(
        go.Pie(
            labels=keys,
            values=[counts[k] for k in keys],
            hole=0.62,
            sort=False,
            direction="clockwise",
            marker=dict(colors=[BASE_STRONG[k] for k in keys], line=dict(width=0)),
            textinfo="label+percent",
            hovertemplate="%{label}: %{value:,} bp (%{percent})<extra></extra>",
            showlegend=False,
        )
    )
    annotations = []
    if gc is not None and not np.isnan(gc):
        annotations.append(
            dict(text=f"GC<br><b>{gc:.1f}%</b>", x=0.5, y=0.5, showarrow=False, font=dict(size=16))
        )
    fig.update_layout(height=300, margin=dict(l=10, r=10, t=10, b=10), annotations=annotations)
    return fig
