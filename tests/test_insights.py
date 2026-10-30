"""Tests for siderascope.insights (rebuilt) and the size-aware helpers in core."""

import math
import random

import pandas as pd
import pytest

from siderascope.core import (
    annotation_agreement,
    fmt_bp,
    recommended_settings,
    window_step,
    windowed_stats,
)
from siderascope.insights import (
    Insight,
    describe,
    gc_class,
    low_gc_stretches,
    skew_origin_terminus,
)


# ------------------------------------------------------------- core helpers
def test_fmt_bp():
    assert fmt_bp(850) == "850 bp"
    assert fmt_bp(4_809) == "4.8 kb"
    assert fmt_bp(6_827_737) == "6.83 Mb"


def test_recommended_settings_scale_with_size():
    short, contig, genome = (recommended_settings(n) for n in (4_809, 250_000, 6_800_000))
    assert (short["scale"], contig["scale"], genome["scale"]) == ("short", "contig", "genome")
    assert short["min_aa"] == 100 and contig["min_aa"] == 150 and genome["min_aa"] == 300
    assert genome["window"] >= 1000
    assert 20 <= short["window"] <= 1000


def test_window_step_keeps_window_count_bounded():
    assert window_step(10_000, 1000) == 250
    assert window_step(6_800_000, 1000) >= 6_800_000 // 4000


def _orf_table(rows):
    return pd.DataFrame(rows, columns=["strand", "start", "end"])


def _cds_table(rows):
    return pd.DataFrame(rows, columns=["strand", "start", "end"])


def test_annotation_agreement_matches_on_stop_codon_and_coverage():
    orfs = _orf_table([("+", 90, 500), ("-", 1000, 1600)])
    cds = _cds_table([
        ("+", 100, 500),    # same stop, ORF starts earlier: matched
        ("-", 1000, 1500),  # same start coordinate, ORF extends further: matched
        ("+", 700, 900),    # nothing there: missed
    ])
    assert annotation_agreement(cds, orfs) == (2, 3)


def test_annotation_agreement_with_no_orfs_or_no_cds():
    empty = _orf_table([])
    assert annotation_agreement(_cds_table([("+", 1, 30)]), empty) == (0, 1)
    assert annotation_agreement(_cds_table([]), _orf_table([("+", 1, 30)])) == (0, 0)


# ---------------------------------------------------------------- insights
def _skewed_genome(n_bins=100, bin_size=1000):
    """Skew is negative in the first half of the genome and positive in the second."""
    half = n_bins // 2
    first = ("C" * 600 + "G" * 400) * (bin_size // 1000)
    second = ("G" * 600 + "C" * 400) * (bin_size // 1000)
    return first * half + second * (n_bins - half)


def _stats(seq, window=1000, step=250):
    return windowed_stats(seq, window, step)


def test_insight_rejects_unknown_levels():
    with pytest.raises(ValueError):
        Insight("danger", "x", "y")


def test_gc_class_boundaries():
    assert gc_class(66.1) == "high"
    assert gc_class(60.0) == "high"
    assert gc_class(50.0) == "moderate"
    assert gc_class(39.9) == "low"
    assert gc_class(float("nan")) == "unknown"


def test_skew_origin_and_terminus_sit_at_the_turning_points():
    seq = _skewed_genome()
    est = skew_origin_terminus(_stats(seq), len(seq))
    assert est["strength"] > 0.05
    assert est["origin"] == pytest.approx(50_000, abs=1_500)  # curve falls, then rises
    assert est["terminus"] > 95_000


def test_random_sequence_has_a_weak_skew_signal():
    random.seed(0)
    seq = "".join(random.choice("ACGT") for _ in range(300_000))
    assert skew_origin_terminus(_stats(seq), len(seq))["strength"] < 0.006


def test_too_few_windows_gives_no_origin_estimate():
    assert skew_origin_terminus(_stats("ACGT" * 50, 20, 20), 200) is None


def test_low_gc_stretches_are_longest_first_and_respect_min_span():
    gc_rich, at_rich = "GC" * 500, "AT" * 500
    seq = gc_rich * 5 + at_rich * 3 + gc_rich * 2 + at_rich * 1 + gc_rich * 4
    stats = _stats(seq, 1000, 1000)
    mean = 100 * sum(ch in "GC" for ch in seq) / len(seq)
    runs = low_gc_stretches(stats, mean, min_span=1000)
    assert runs == [(5001, 8000), (10001, 11000)]
    assert low_gc_stretches(stats, mean, min_span=2000) == [(5001, 8000)]


def test_a_low_gc_stretch_at_the_very_end_is_kept():
    seq = ("GC" * 500) * 5 + ("AT" * 500) * 2
    stats = _stats(seq, 1000, 1000)
    mean = 100 * sum(ch in "GC" for ch in seq) / len(seq)
    assert low_gc_stretches(stats, mean, min_span=1000) == [(5001, 7000)]


def _describe(length, gc, n_orfs, seq=None, **kw):
    seq = seq or _skewed_genome()
    return describe(length=length, gc=gc, ambiguous=kw.pop("ambiguous", 0), n_orfs=n_orfs,
                    min_aa=kw.pop("min_aa", 300), stats=_stats(seq), **kw)


def test_inflated_orf_count_is_a_warning_that_names_a_next_step():
    notes = _describe(6_800_000, 66.1, 22_000)
    warn = [n for n in notes if n.level == "warning"]
    assert len(warn) == 1
    assert "22,000 ORFs" in warn[0].title
    assert "Prodigal" in warn[0].text


def test_reasonable_orf_count_is_reported_as_good():
    notes = _describe(5_000_000, 55.0, 4_800)
    assert not [n for n in notes if n.level == "warning"]
    assert any(n.level == "good" and "4,800 ORFs" in n.title for n in notes)


def test_small_sequence_gets_no_genome_level_comments():
    seq = "ACGT" * 1000
    notes = describe(length=4_000, gc=50.0, ambiguous=0, n_orfs=3, min_aa=100, stats=_stats(seq, 100, 100))
    assert [n.title for n in notes] == ["50.0% GC is moderate"]


def test_no_orfs_suggests_a_lower_minimum():
    seq = "ACGT" * 1000
    notes = describe(length=4_000, gc=50.0, ambiguous=0, n_orfs=0, min_aa=100, stats=_stats(seq, 100, 100))
    assert any("No ORFs" in n.title and "lower" in n.text for n in notes)


def test_all_n_sequence_has_no_gc_note():
    seq = "N" * 1000
    notes = describe(length=1_000, gc=math.nan, ambiguous=1_000, n_orfs=0, min_aa=100, stats=_stats(seq, 100, 100))
    assert not any("GC is" in n.title for n in notes)
    assert any("ambiguous" in n.title for n in notes)


def test_ambiguous_bases_are_a_warning_only_when_common():
    seq = "ACGT" * 1000
    few = describe(length=4_000, gc=50.0, ambiguous=3, n_orfs=1, min_aa=100, stats=_stats(seq, 100, 100))
    many = describe(length=4_000, gc=50.0, ambiguous=400, n_orfs=1, min_aa=100, stats=_stats(seq, 100, 100))
    level = lambda notes: next(n.level for n in notes if "ambiguous" in n.title)
    assert (level(few), level(many)) == ("info", "warning")


def test_annotation_agreement_notes():
    seq = "ACGT" * 1000
    kw = dict(length=4_000, gc=50.0, ambiguous=0, n_orfs=3, min_aa=100, stats=_stats(seq, 100, 100))
    full = describe(agreement=(5, 5), **kw)
    part = describe(agreement=(3, 5), **kw)
    assert any(n.level == "good" and "All 5" in n.title for n in full)
    assert any("3 of 5" in n.title for n in part)


def test_origin_note_appears_for_a_skewed_genome_and_not_for_a_flat_one():
    skewed = _describe(100_000, 50.0, 90)
    assert any("Possible replication origin" in n.title for n in skewed)
    flat_seq = ("ACGT" * 250) * 100  # perfectly balanced, no skew
    flat = _describe(100_000, 50.0, 90, seq=flat_seq)
    assert any("no clear origin" in n.title for n in flat)


def test_at_rich_stretch_note():
    seq = ("GC" * 500) * 60 + ("AT" * 500) * 10 + ("GC" * 500) * 60
    stats = _stats(seq, 1000, 1000)
    gc = 100 * sum(ch in "GC" for ch in seq) / len(seq)
    notes = describe(length=len(seq), gc=gc, ambiguous=0, n_orfs=100, min_aa=300, stats=stats)
    assert any("AT-rich stretch" in n.title for n in notes)
