"""End-to-end checks on the bundled synthetic example files."""

from pathlib import Path

import pytest

from siderascope import plots
from siderascope.core import find_orfs, parse_sequences, windowed_stats

EXAMPLES = Path(__file__).resolve().parent.parent / "examples"


@pytest.fixture(scope="module")
def genbank_record():
    return parse_sequences((EXAMPLES / "example.gb").read_text())[0]


def test_fasta_and_genbank_examples_hold_the_same_sequence(genbank_record):
    fasta = parse_sequences((EXAMPLES / "example.fasta").read_text())[0]
    assert fasta.sequence == genbank_record.sequence
    assert len(fasta.cds) == 0
    assert len(genbank_record.cds) == 5


def test_every_annotated_cds_is_recovered_as_an_orf(genbank_record):
    """The stop codon must match exactly and the ORF must cover the annotated gene.

    The start can differ: SideroScope reports the outermost in-frame start, which may
    sit upstream of the annotated one (a known limitation of simple ORF finding).
    """
    orfs = find_orfs(genbank_record.sequence, min_aa=100)
    for cds in genbank_record.cds.itertuples():
        same_strand = orfs[orfs["strand"] == cds.strand]
        if cds.strand == "+":
            hits = same_strand[(same_strand["end"] == cds.end) & (same_strand["start"] <= cds.start)]
        else:
            hits = same_strand[(same_strand["start"] == cds.start) & (same_strand["end"] >= cds.end)]
        assert len(hits) == 1, f"annotated CDS {cds.locus_tag} was not recovered"


def test_figures_build_for_example(genbank_record):
    seq = genbank_record.sequence
    orfs = find_orfs(seq, min_aa=100)
    stats = windowed_stats(seq, 500, 125)
    assert len(plots.base_view_figure(seq, 1, 600, 60).data) == 1
    assert len(plots.gc_figure(stats, 55.0).data) == 3
    assert len(plots.orf_map_figure(orfs, len(seq), genbank_record.cds, 100).data) >= 3
    assert len(plots.composition_figure({"A": 1, "C": 2, "G": 3, "T": 4, "other": 0}).data) == 1


def test_orf_map_handles_no_orfs_and_no_annotations():
    empty = find_orfs("ACGT" * 10, min_aa=100)
    fig = plots.orf_map_figure(empty, 40, None, 100)
    assert any("No ORFs" in a.text for a in fig.layout.annotations)


def test_base_view_handles_short_sequences_and_ambiguity():
    fig = plots.base_view_figure("ACGTNNRY", start=1, length=600, line_width=60)
    assert len(fig.data) == 1
