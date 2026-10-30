"""SideroScope: a Streamlit app for exploring DNA sequences.

Run with:  streamlit run app.py
"""

from __future__ import annotations

import math
import textwrap
from pathlib import Path

import streamlit as st

from siderascope import __version__, plots, ui
from siderascope.core import (
    annotation_agreement,
    find_orfs,
    fmt_bp,
    orfs_to_fasta,
    parse_sequences,
    recommended_settings,
    sequence_summary,
    window_step,
    windowed_stats,
)
from siderascope.insights import describe

ROOT = Path(__file__).parent
EXAMPLES = {
    "Synthetic demo (FASTA)": ROOT / "examples" / "example.fasta",
    "Synthetic demo with annotations (GenBank)": ROOT / "examples" / "example.gb",
}
UPLOAD_TYPES = ["fasta", "fa", "fna", "fas", "gb", "gbk", "genbank", "txt"]
LARGE_SEQUENCE = 2_000_000

st.set_page_config(page_title="SideroScope", page_icon="🧬", layout="wide")
ui.inject_css()

ss = st.session_state
ss.setdefault("mode", None)


def set_mode(mode):
    ss["mode"] = mode


@st.cache_data(show_spinner=False)
def cached_parse(text: str):
    return parse_sequences(text)


@st.cache_data(show_spinner=False)
def cached_windows(seq: str, window: int, step: int):
    return windowed_stats(seq, window, step)


@st.cache_data(show_spinner="Scanning the six reading frames...")
def cached_orfs(seq: str, min_aa: int, table: int, alt_starts: bool):
    return find_orfs(seq, min_aa=min_aa, table=table, alt_starts=alt_starts)


# --------------------------------------------------------------------------- #
# Landing page
# --------------------------------------------------------------------------- #
if ss["mode"] is None:
    ui.hero()
    cols = st.columns(3, gap="medium")
    with cols[0].container(border=True):
        ui.action_card("spark", "Start with the example", "A small synthetic sequence with annotated genes. See every feature in seconds.", ui.TEAL)
        st.button("Open the example", key="go_example", type="primary", on_click=set_mode, args=("example",), width="stretch")
    with cols[1].container(border=True):
        ui.action_card("upload", "Upload your sequence", "A FASTA or GenBank file, from a single gene up to a whole bacterial genome.", ui.BLUE)
        st.button("Choose a file", key="go_upload", on_click=set_mode, args=("upload",), width="stretch")
    with cols[2].container(border=True):
        ui.action_card("paste", "Paste a sequence", "Copy bases straight from a paper, a database or your notes.", ui.PURPLE)
        st.button("Paste text", key="go_paste", on_click=set_mode, args=("paste",), width="stretch")
    ui.feature_tiles()
    ui.footer(__version__)
    st.stop()

mode = ss["mode"]

# --------------------------------------------------------------------------- #
# Input
# --------------------------------------------------------------------------- #
with st.sidebar:
    ui.sidebar_brand(__version__)
    st.button("Back to start", key="back", icon=":material/arrow_back:", on_click=set_mode, args=(None,), width="stretch")

text = None
if mode == "example":
    with st.sidebar:
        choice = st.selectbox("Example sequence", list(EXAMPLES), key="example_choice")
    text = EXAMPLES[choice].read_text()
elif mode == "upload":
    loaded = ss.get("uploader") is not None
    with st.expander("Input file", expanded=not loaded, icon=":material/upload_file:"):
        uploaded = st.file_uploader("Drop a FASTA or GenBank file here", type=UPLOAD_TYPES, key="uploader")
    if uploaded is not None:
        text = uploaded.getvalue().decode("utf-8", errors="replace")
else:
    has_text = bool(ss.get("pasted", "").strip())
    with st.expander("Pasted sequence", expanded=not has_text, icon=":material/content_paste:"):
        text = st.text_area(
            "Paste FASTA or a raw DNA sequence",
            key="pasted",
            height=180,
            placeholder=">my_sequence\nATGGCT...",
        )

if not text or not text.strip():
    ui.waiting_card(mode)
    st.stop()

try:
    records = cached_parse(text)
except ValueError as err:
    st.error(str(err), icon=":material/error:")
    st.stop()

# --------------------------------------------------------------------------- #
# Sidebar: record and settings
# --------------------------------------------------------------------------- #
with st.sidebar:
    if len(records) > 1:
        labels = [f"{r.id} ({r.length:,} bp)" for r in records]
        record = records[st.selectbox("Record", range(len(records)), format_func=lambda i: labels[i])]
    else:
        record = records[0]

    seq = record.sequence
    n = len(seq)
    suggested = recommended_settings(n)
    tag = f"{record.id}_{n}"  # new defaults whenever the record changes

    st.markdown("#### Settings")
    st.caption(f"Starting values are suggested for a {suggested['label']} ({fmt_bp(n)}).")
    max_window = max(n, 10)
    window = st.number_input(
        "GC window (bp)",
        min_value=10,
        max_value=max_window,
        value=min(suggested["window"], max_window),
        step=10 ** max(1, len(str(suggested["window"])) - 2),
        key=f"window_{tag}",
        help="Width of the sliding window for GC content and skew. Larger windows give smoother curves.",
    )
    min_aa = st.number_input(
        "Minimum ORF length (aa)",
        min_value=10,
        max_value=2000,
        value=suggested["min_aa"],
        step=10,
        key=f"min_aa_{tag}",
        help="Shorter ORFs are ignored. Higher values remove more chance ORFs.",
    )
    table = st.selectbox(
        "Genetic code",
        [11, 1],
        format_func=lambda t: {11: "11 - Bacterial / archaeal", 1: "1 - Standard"}[t],
        help="Table 11 is the bacterial, archaeal and plastid code. Both tables use the same stop codons.",
    )
    alt_starts = st.checkbox(
        "Also allow GTG and TTG start codons",
        value=False,
        help="Many bacterial genes start with GTG or TTG, but allowing them finds many more chance ORFs.",
    )

if n > LARGE_SEQUENCE:
    st.toast(f"{fmt_bp(n)} is a large sequence, so the ORF scan may take a few seconds.", icon=":material/hourglass_top:")

# --------------------------------------------------------------------------- #
# Analysis
# --------------------------------------------------------------------------- #
summary = sequence_summary(seq)
step = window_step(n, int(window))
stats = cached_windows(seq, int(window), step)
orfs = cached_orfs(seq, int(min_aa), int(table), bool(alt_starts))
agreement = annotation_agreement(record.cds, orfs) if len(record.cds) else None
insights = describe(
    length=n,
    gc=summary["gc"],
    ambiguous=summary["ambiguous"],
    n_orfs=len(orfs),
    min_aa=int(min_aa),
    stats=stats,
    agreement=agreement,
)
gc_text = "n/a" if math.isnan(summary["gc"]) else f"{summary['gc']:.1f}%"

chips = [fmt_bp(n), f"GC {gc_text}"]
if len(records) > 1:
    chips.append(f"{len(records)} records")
if len(record.cds):
    chips.append(f"{len(record.cds):,} annotated CDS")
description = record.description
if description.startswith(record.id):  # FASTA headers repeat the ID before the description
    description = description[len(record.id) :].strip(" |:-")
ui.banner(record.id, description, chips)

tab_overview, tab_sequence, tab_gc, tab_orfs, tab_about = st.tabs(
    ["Overview", "Sequence", "GC content & skew", "ORFs", "About"]
)

# --------------------------------------------------------------------------- #
# Overview
# --------------------------------------------------------------------------- #
with tab_overview:
    per_kb = len(orfs) / (n / 1000)
    cards = [
        ("Length", fmt_bp(n), f"{n:,} bp", ui.TEAL),
        ("GC content", gc_text, "of unambiguous bases", ui.ORANGE),
        (f"Candidate ORFs, {int(min_aa)}+ aa", f"{len(orfs):,}", f"{per_kb:.2f} per kb, not a gene count", ui.BLUE),
    ]
    if len(record.cds):
        cards.append(("Annotated CDS", f"{len(record.cds):,}", "from the GenBank file", ui.PURPLE))
    else:
        cards.append(("Ambiguous bases", f"{summary['ambiguous']:,}", "N, R, Y and similar codes", ui.PURPLE))
    ui.stat_cards(cards)

    left, right = st.columns([3, 2], gap="large")
    with left:
        ui.section("Genome map", "From the outside in: ORFs on each strand, GC content, GC skew")
        st.plotly_chart(plots.genome_ring_figure(stats, orfs, n, summary["gc"], record.id))
        st.caption(
            "GC skew is scaled to its largest value. Tick marks show position along the sequence, "
            "clockwise from the top."
        )
    with right:
        ui.section("What this means", "A plain-language reading of your numbers")
        ui.insight_cards(insights)
        ui.section("Base composition")
        st.plotly_chart(plots.composition_figure(summary["counts"], summary["gc"]))

# --------------------------------------------------------------------------- #
# Sequence view
# --------------------------------------------------------------------------- #
with tab_sequence:
    ui.section("Sequence view", "Each base is a coloured tile. Hover for its exact position.")
    st.markdown(
        ":green[**A**] &nbsp; :blue[**C**] &nbsp; :orange[**G**] &nbsp; :red[**T**] &nbsp; "
        "grey = ambiguous base (N, R, Y ...)"
    )

    orf_starts = dict(zip(orfs["orf_id"], orfs["start"]))
    orf_choices = ["-"] + list(orf_starts)[:5000]

    def jump_to_orf():
        pick = ss.get("jump_orf")
        if pick in orf_starts:
            ss["view_start"] = max(1, int(orf_starts[pick]) - 30)

    ss["view_start"] = min(max(1, int(ss.get("view_start", 1))), n)
    c_a, c_b, c_c, c_d = st.columns(4)
    c_a.number_input("Start position (bp)", min_value=1, max_value=n, step=60, key="view_start")
    c_b.selectbox("Jump to an ORF", orf_choices, key="jump_orf", on_change=jump_to_orf)
    view_len = c_c.select_slider("Bases to show", options=[120, 300, 600, 1200, 2400], value=600)
    line_width = c_d.select_slider("Bases per line", options=[30, 45, 60, 90, 120], value=60)
    st.plotly_chart(plots.base_view_figure(seq, int(ss["view_start"]), int(view_len), int(line_width)))

# --------------------------------------------------------------------------- #
# GC content and skew
# --------------------------------------------------------------------------- #
with tab_gc:
    ui.section("GC content and skew", f"Window {min(int(window), n):,} bp, step {step:,} bp")
    st.caption(
        "The dashed line is the whole-sequence GC content. Dips can mark regions acquired from other "
        "organisms. In many bacterial chromosomes the cumulative skew changes direction near the "
        "origin and the terminus of replication."
    )
    st.plotly_chart(plots.gc_figure(stats, summary["gc"]))
    with st.expander("Window data"):
        st.dataframe(stats, hide_index=True)
        st.download_button(
            "Download window data (CSV)",
            stats.to_csv(index=False),
            file_name=f"{record.id}_gc_windows.csv",
            mime="text/csv",
        )

# --------------------------------------------------------------------------- #
# ORFs
# --------------------------------------------------------------------------- #
with tab_orfs:
    ui.section("Open reading frames", "Candidates on all six frames. These are not gene predictions.")

    if orfs.empty:
        st.info(f"No ORFs of {int(min_aa)} aa or longer were found. Try lowering the minimum length.")
        st.plotly_chart(plots.orf_map_figure(orfs, n, record.cds, int(min_aa)))
    else:
        f_strand, f_len = st.columns([1, 2], gap="large")
        strand_pick = f_strand.radio("Strand", ["Both", "+", "-"], horizontal=True)
        lo, hi = int(orfs["length_aa"].min()), int(orfs["length_aa"].max())
        len_range = f_len.slider("ORF length (aa)", lo, hi, (lo, hi)) if hi > lo else (lo, hi)

        mask = orfs["length_aa"].between(len_range[0], len_range[1])
        if strand_pick != "Both":
            mask &= orfs["strand"] == strand_pick
        view = orfs[mask].reset_index(drop=True)
        st.caption(f"Showing {len(view):,} of {len(orfs):,} ORFs")

        st.plotly_chart(plots.orf_map_figure(view, n, record.cds, int(min_aa)))

        with st.expander("Length distribution"):
            st.plotly_chart(plots.orf_length_histogram(view if len(view) else orfs))

        st.markdown("**Click a row to see its protein sequence.**")
        event = st.dataframe(
            view.drop(columns="protein"),
            hide_index=True,
            on_select="rerun",
            selection_mode="single-row",
            key="orf_table",
            column_config={
                "length_aa": st.column_config.ProgressColumn(
                    "length_aa", min_value=0, max_value=int(hi), format="%d"
                )
            },
        )
        picked = event.selection.rows if event is not None else []
        if picked and picked[0] < len(view):
            row = view.iloc[picked[0]]
            ui.section(
                f"{row['orf_id']}: {int(row['length_aa'])} aa",
                f"strand {row['strand']}, {int(row['start']):,}-{int(row['end']):,} bp, "
                f"starts with {row['start_codon']}, ends with {row['stop_codon']}",
            )
            st.code(textwrap.fill(row["protein"], 60), language=None)

        d1, d2 = st.columns(2)
        d1.download_button(
            "Download ORF table (CSV)",
            view.to_csv(index=False),
            file_name=f"{record.id}_orfs.csv",
            mime="text/csv",
            icon=":material/download:",
            width="stretch",
        )
        d2.download_button(
            "Download proteins (FASTA)",
            orfs_to_fasta(view, record.id),
            file_name=f"{record.id}_orf_proteins.fasta",
            mime="text/plain",
            icon=":material/download:",
            width="stretch",
        )

    if len(record.cds):
        with st.expander(f"Annotated CDS from the GenBank file ({len(record.cds):,})"):
            st.dataframe(record.cds, hide_index=True)

# --------------------------------------------------------------------------- #
# About
# --------------------------------------------------------------------------- #
with tab_about:
    ui.section("About SideroScope", f"Version {__version__}")
    st.markdown(
        """
**What the numbers mean**

- **GC content** is (G + C) / (A + C + G + T) x 100. Ambiguous bases (N, R, Y ...) are ignored.
- **GC skew** is (G - C) / (G + C) in each window. **Cumulative skew** is its running sum; in many
  bacterial chromosomes the slope changes near the origin and terminus of replication.
- **ORFs** run from a start codon (ATG, optionally GTG and TTG) to the first in-frame stop codon
  (TAA, TAG, TGA) on any of the six frames. When several starts share a stop, the longest ORF is kept.
  ORFs that run off the end of the sequence without a stop are not reported.
- Frames **+1 to +3** start at bases 1 to 3 of the forward strand. Frames **-1 to -3** read the
  reverse complement, starting from the last bases of the forward strand.

**Limits**

- ORFs are *candidates*. Random sequence also contains ORFs, especially short ones, and G/C-rich genomes
  have many long ones because stop codons are AT-rich. For real gene calling use a dedicated tool such as
  Prodigal, Glimmer or Bakta.
- The "What this means" notes are rules of thumb for typical bacterial sequences, not diagnoses.
- Outlined boxes on the ORF map are CDS annotated in a GenBank file; filled boxes are predicted here.
- Sequences with ambiguity codes are supported, but ambiguous bases break codons that contain them.

SideroScope is open source (MIT). Ideas and bug reports are welcome as GitHub issues.
Fonts (Inter, Space Grotesk, JetBrains Mono) are bundled under the SIL Open Font License.
"""
    )

ui.footer(__version__)
