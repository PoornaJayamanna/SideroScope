import io
import math

import numpy as np
import pytest
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord

from siderascope.core import (
    base_composition,
    clean_sequence,
    find_orfs,
    gc_content,
    orfs_to_fasta,
    parse_sequences,
    sequence_summary,
    windowed_stats,
)


def revcomp(s: str) -> str:
    return str(Seq(s).reverse_complement())


# ---------------------------------------------------------------- reading input
def test_clean_sequence_strips_whitespace_numbers_and_uppercases():
    assert clean_sequence("acg t\n  1 gcA") == "ACGTGCA"


def test_clean_sequence_converts_rna_to_dna():
    assert clean_sequence("ACGU") == "ACGT"


def test_clean_sequence_rejects_protein():
    with pytest.raises(ValueError, match="protein"):
        clean_sequence("MKVLEEPQ")


def test_parse_fasta_with_two_records():
    recs = parse_sequences(">a first\nACGT\nACGT\n>b second\nTTTT\n")
    assert [r.id for r in recs] == ["a", "b"]
    assert recs[0].sequence == "ACGTACGT"
    assert recs[1].description == "b second"


def test_parse_raw_sequence():
    recs = parse_sequences("  acgt acgt \n")
    assert len(recs) == 1
    assert recs[0].sequence == "ACGTACGT"
    assert recs[0].id == "pasted_sequence"


def test_parse_empty_raises():
    with pytest.raises(ValueError):
        parse_sequences("   \n ")


def test_parse_genbank_reads_cds_features():
    rec = SeqRecord(
        Seq("ATG" + "GCT" * 10 + "TAA" + "CCCC"),
        id="TEST.1",
        name="TEST",
        description="test record",
        annotations={"molecule_type": "DNA"},
        features=[
            SeqFeature(
                FeatureLocation(0, 36, strand=-1),
                type="CDS",
                qualifiers={"gene": ["abc"], "locus_tag": ["T_0001"], "product": ["a protein"]},
            )
        ],
    )
    handle = io.StringIO()
    SeqIO.write(rec, handle, "genbank")

    parsed = parse_sequences(handle.getvalue())
    assert len(parsed) == 1
    cds = parsed[0].cds
    assert len(cds) == 1
    row = cds.iloc[0]
    assert (row["start"], row["end"], row["strand"]) == (1, 36, "-")
    assert (row["gene"], row["locus_tag"], row["product"]) == ("abc", "T_0001", "a protein")


# ------------------------------------------------------- composition, GC, skew
def test_composition_and_gc_ignore_ambiguous_bases():
    seq = "AACCGGTTNN"
    assert base_composition(seq) == {"A": 2, "C": 2, "G": 2, "T": 2, "other": 2}
    assert gc_content(seq) == 50.0
    assert sequence_summary(seq)["ambiguous"] == 2


def test_gc_of_all_n_is_nan():
    assert math.isnan(gc_content("NNNN"))


def test_windowed_stats_non_overlapping():
    stats = windowed_stats("GGGGCCCC", window=4, step=4)
    assert list(stats["start"]) == [1, 5]
    assert list(stats["end"]) == [4, 8]
    assert list(stats["gc"]) == [100.0, 100.0]
    assert list(stats["skew"]) == [1.0, -1.0]
    assert list(stats["cum_skew"]) == [1.0, 0.0]


def test_windowed_stats_overlapping_scales_cumulative_skew():
    stats = windowed_stats("GGGGCCCC", window=4, step=2)
    assert list(stats["skew"]) == [1.0, 0.0, -1.0]
    assert np.allclose(stats["cum_skew"], [0.5, 0.5, 0.0])


def test_windowed_stats_window_larger_than_sequence_uses_whole_sequence():
    stats = windowed_stats("GGCCAT", window=1000, step=250)
    assert len(stats) == 1
    assert (stats.loc[0, "start"], stats.loc[0, "end"]) == (1, 6)
    assert stats.loc[0, "gc"] == pytest.approx(100 * 4 / 6)


def test_windowed_stats_no_gc_gives_nan_skew():
    stats = windowed_stats("AAAATTTT", window=4, step=4)
    assert stats["skew"].isna().all()
    assert list(stats["gc"]) == [0.0, 0.0]


# ------------------------------------------------------------------------ ORFs
ORF = "ATG" + "GCT" * 10 + "TAA"  # 11 codons + stop, 36 nt


def test_forward_orf_coordinates_frame_and_protein():
    seq = "CC" + ORF + "GG"
    orfs = find_orfs(seq, min_aa=5)
    assert len(orfs) == 1
    row = orfs.iloc[0]
    assert row["strand"] == "+"
    assert (row["start"], row["end"]) == (3, 38)
    assert row["frame"] == 3
    assert row["length_nt"] == 36
    assert row["length_aa"] == 11
    assert row["protein"] == "M" + "A" * 10
    assert (row["start_codon"], row["stop_codon"]) == ("ATG", "TAA")


def test_reverse_strand_orf_maps_back_to_forward_coordinates():
    forward = "CC" + ORF + "GG"
    seq = revcomp(forward)
    n = len(seq)
    row = find_orfs(seq, min_aa=5).iloc[0]
    assert row["strand"] == "-"
    assert (row["start"], row["end"]) == (n - 38 + 1, n - 3 + 1)
    assert row["protein"] == "M" + "A" * 10
    assert row["frame"] < 0


def test_min_aa_filters_and_empty_result_keeps_columns():
    orfs = find_orfs("CC" + ORF + "GG", min_aa=50)
    assert orfs.empty
    assert "protein" in orfs.columns


def test_orf_without_stop_codon_is_not_reported():
    assert find_orfs("ATG" + "GCT" * 10, min_aa=5).empty


def test_outermost_start_is_used_when_starts_are_nested():
    seq = "ATG" + "ATG" + "GCT" * 10 + "TAA"
    orfs = find_orfs(seq, min_aa=5)
    assert len(orfs) == 1
    assert orfs.iloc[0]["start"] == 1
    assert orfs.iloc[0]["length_aa"] == 12


def test_alternative_start_codons_are_opt_in():
    seq = "GTG" + "GCT" * 10 + "TAA"
    assert find_orfs(seq, min_aa=5).empty
    row = find_orfs(seq, min_aa=5, alt_starts=True).iloc[0]
    assert row["start_codon"] == "GTG"
    assert row["protein"].startswith("M")


def test_unsupported_genetic_code_raises():
    with pytest.raises(ValueError):
        find_orfs("ATGTAA", table=4)


def test_orfs_to_fasta_format():
    orfs = find_orfs("CC" + ORF + "GG", min_aa=5)
    text = orfs_to_fasta(orfs, "rec1")
    lines = text.strip().splitlines()
    assert lines[0] == ">rec1|ORF_0001|strand=+|3-38|11aa"
    assert lines[1] == "M" + "A" * 10
    assert orfs_to_fasta(find_orfs("ACGT", min_aa=5)) == ""
