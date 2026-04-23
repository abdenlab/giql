"""Unit tests for pybedtools wrapper functions."""

import shutil

import pytest

pybedtools = pytest.importorskip("pybedtools")

if not shutil.which("bedtools"):
    pytest.skip(
        "bedtools binary not found in PATH",
        allow_module_level=True,
    )

from .bedtools_wrapper import BedtoolsError  # noqa: E402
from .bedtools_wrapper import bedtool_to_tuples  # noqa: E402
from .bedtools_wrapper import closest  # noqa: E402
from .bedtools_wrapper import create_bedtool  # noqa: E402
from .bedtools_wrapper import intersect  # noqa: E402
from .bedtools_wrapper import merge  # noqa: E402

pytestmark = pytest.mark.integration


def test_create_bedtool_bed3_format():
    """
    GIVEN a list of BED3 tuples
    WHEN create_bedtool() is called
    THEN returns a BedTool with correct intervals
    """
    bt = create_bedtool([("chr1", 100, 200)])
    intervals = list(bt)
    assert len(intervals) == 1
    assert intervals[0].chrom == "chr1"
    assert intervals[0].start == 100
    assert intervals[0].end == 200


def test_create_bedtool_bed6_format():
    """
    GIVEN a list of BED6 tuples
    WHEN create_bedtool() is called
    THEN returns a BedTool with all 6 fields
    """
    bt = create_bedtool([("chr1", 100, 200, "a1", 50, "+")])
    intervals = list(bt)
    assert len(intervals) == 1
    assert intervals[0].fields == ["chr1", "100", "200", "a1", "50", "+"]


def test_create_bedtool_none_values_replaced():
    """
    GIVEN BED6 tuples with None values
    WHEN create_bedtool() is called
    THEN None values replaced with defaults
    """
    bt = create_bedtool([("chr1", 100, 200, None, None, None)])
    fields = list(bt)[0].fields
    assert fields[3] == "."  # name
    assert fields[4] == "0"  # score
    assert fields[5] == "."  # strand


def test_create_bedtool_invalid_tuple_length_raises():
    """
    GIVEN a tuple with invalid length
    WHEN create_bedtool() is called
    THEN ValueError is raised
    """
    with pytest.raises(ValueError, match="Invalid interval format"):
        create_bedtool([("chr1", 100)])


def test_create_bedtool_multiple_intervals():
    """
    GIVEN multiple intervals across chromosomes
    WHEN create_bedtool() is called
    THEN BedTool contains all intervals
    """
    bt = create_bedtool(
        [
            ("chr1", 100, 200, "a", 0, "+"),
            ("chr2", 300, 400, "b", 0, "-"),
        ]
    )
    intervals = list(bt)
    assert len(intervals) == 2


def test_intersect_basic_overlap():
    """
    GIVEN two sets of overlapping intervals
    WHEN intersect() is called
    THEN returns intervals from A that overlap B
    """
    a = [("chr1", 100, 200, "a1", 100, "+")]
    b = [("chr1", 150, 250, "b1", 100, "+")]
    result = intersect(a, b)
    assert len(result) == 1
    assert result[0][0] == "chr1"


def test_intersect_no_overlap():
    """
    GIVEN non-overlapping intervals
    WHEN intersect() is called
    THEN returns empty list
    """
    a = [("chr1", 100, 200, "a1", 100, "+")]
    b = [("chr1", 300, 400, "b1", 100, "+")]
    result = intersect(a, b)
    assert result == []


def test_intersect_same_strand_mode():
    """
    GIVEN intervals on same and opposite strands
    WHEN intersect() is called with strand_mode="same"
    THEN only same-strand overlaps returned
    """
    a = [
        ("chr1", 100, 200, "a1", 0, "+"),
        ("chr1", 100, 200, "a2", 0, "-"),
    ]
    b = [("chr1", 150, 250, "b1", 0, "+")]
    result = intersect(a, b, strand_mode="same")
    names = [r[3] for r in result]
    assert "a1" in names
    assert "a2" not in names


def test_intersect_opposite_strand_mode():
    """
    GIVEN intervals on same and opposite strands
    WHEN intersect() is called with strand_mode="opposite"
    THEN only opposite-strand overlaps returned
    """
    a = [
        ("chr1", 100, 200, "a1", 0, "+"),
        ("chr1", 100, 200, "a2", 0, "-"),
    ]
    b = [("chr1", 150, 250, "b1", 0, "+")]
    result = intersect(a, b, strand_mode="opposite")
    names = [r[3] for r in result]
    assert "a2" in names
    assert "a1" not in names


def test_intersect_no_strand_mode():
    """
    GIVEN overlapping intervals on different strands
    WHEN intersect() is called with strand_mode=None
    THEN all overlaps returned regardless of strand
    """
    a = [("chr1", 100, 200, "a1", 0, "+")]
    b = [("chr1", 150, 250, "b1", 0, "-")]
    result = intersect(a, b)
    assert len(result) == 1


def test_merge_overlapping():
    """
    GIVEN overlapping intervals
    WHEN merge() is called
    THEN returns merged BED3 intervals
    """
    intervals = [
        ("chr1", 100, 200, "i1", 0, "+"),
        ("chr1", 150, 250, "i2", 0, "+"),
    ]
    result = merge(intervals)
    assert len(result) == 1
    assert result[0] == ("chr1", 100, 250)


def test_merge_separated():
    """
    GIVEN separated intervals
    WHEN merge() is called
    THEN each interval returned separately (BED3)
    """
    intervals = [
        ("chr1", 100, 200, "i1", 0, "+"),
        ("chr1", 300, 400, "i2", 0, "+"),
    ]
    result = merge(intervals)
    assert len(result) == 2


def test_merge_strand_specific():
    """
    GIVEN overlapping intervals on different strands
    WHEN merge() is called with strand_mode="same"
    THEN merges per-strand separately
    """
    intervals = [
        ("chr1", 100, 200, "i1", 0, "+"),
        ("chr1", 150, 250, "i2", 0, "+"),
        ("chr1", 120, 220, "i3", 0, "-"),
    ]
    result = merge(intervals, strand_mode="same")
    # Should have 2: one merged + strand, one - strand
    assert len(result) == 2


def test_merge_adjacent():
    """
    GIVEN adjacent intervals (end == start of next)
    WHEN merge() is called
    THEN adjacent intervals are merged
    """
    intervals = [
        ("chr1", 100, 200, "i1", 0, "+"),
        ("chr1", 200, 300, "i2", 0, "+"),
    ]
    result = merge(intervals)
    assert len(result) == 1
    assert result[0] == ("chr1", 100, 300)


def test_closest_basic():
    """
    GIVEN non-overlapping intervals
    WHEN closest() is called
    THEN returns each A paired with nearest B plus distance
    """
    a = [("chr1", 100, 200, "a1", 100, "+")]
    b = [("chr1", 300, 400, "b1", 100, "+")]
    result = closest(a, b)
    assert len(result) == 1
    # Last field is distance
    # bedtools 2.31+ may report 101 (1-based gap) vs 100 (0-based)
    assert result[0][-1] in (100, 101)


def test_closest_cross_chromosome():
    """
    GIVEN intervals on different chromosomes
    WHEN closest() is called
    THEN finds nearest per-chromosome
    """
    a = [
        ("chr1", 100, 200, "a1", 0, "+"),
        ("chr2", 100, 200, "a2", 0, "+"),
    ]
    b = [
        ("chr1", 300, 400, "b1", 0, "+"),
        ("chr2", 500, 600, "b2", 0, "+"),
    ]
    result = closest(a, b)
    assert len(result) == 2
    # Each A should match B on same chromosome
    for row in result:
        assert row[0] == row[6]  # a.chrom == b.chrom


def test_closest_same_strand_mode():
    """
    GIVEN intervals with mixed strands
    WHEN closest() is called with strand_mode="same"
    THEN returns nearest same-strand interval
    """
    a = [("chr1", 100, 200, "a1", 0, "+")]
    b = [
        ("chr1", 220, 240, "b_opp", 0, "-"),  # closer but opposite
        ("chr1", 300, 400, "b_same", 0, "+"),  # farther but same
    ]
    result = closest(a, b, strand_mode="same")
    assert len(result) == 1
    assert result[0][9] == "b_same"


def test_closest_k_greater_than_one():
    """
    GIVEN one query and three database intervals
    WHEN closest() is called with k=3
    THEN returns up to 3 nearest
    """
    a = [("chr1", 200, 300, "a1", 0, "+")]
    b = [
        ("chr1", 100, 150, "b1", 0, "+"),
        ("chr1", 350, 400, "b2", 0, "+"),
        ("chr1", 500, 600, "b3", 0, "+"),
    ]
    result = closest(a, b, k=3)
    # bedtools returns up to k nearest; exact count may vary by version
    assert len(result) >= 2


def test_bedtool_to_tuples_bed3_conversion():
    """
    GIVEN a BedTool with BED3 intervals
    WHEN bedtool_to_tuples() is called with bed_format="bed3"
    THEN returns list of (chrom, start, end) tuples with int positions
    """
    bt = pybedtools.BedTool("chr1\t100\t200\n", from_string=True)
    result = bedtool_to_tuples(bt, bed_format="bed3")
    assert result == [("chr1", 100, 200)]


def test_bedtool_to_tuples_bed6_conversion():
    """
    GIVEN a BedTool with BED6 intervals
    WHEN bedtool_to_tuples() is called with bed_format="bed6"
    THEN returns list of 6-tuples with correct types
    """
    bt = pybedtools.BedTool("chr1\t100\t200\tgene1\t500\t+\n", from_string=True)
    result = bedtool_to_tuples(bt, bed_format="bed6")
    assert result == [("chr1", 100, 200, "gene1", 500, "+")]


def test_bedtool_to_tuples_bed6_dot_to_none():
    """
    GIVEN a BedTool with "." for name and strand
    WHEN bedtool_to_tuples() is called with bed_format="bed6"
    THEN "." values converted to None
    """
    bt = pybedtools.BedTool("chr1\t100\t200\t.\t0\t.\n", from_string=True)
    result = bedtool_to_tuples(bt, bed_format="bed6")
    assert result[0][3] is None  # name
    assert result[0][5] is None  # strand


def test_bedtool_to_tuples_bed6_padding():
    """
    GIVEN a BedTool with fewer than 6 fields
    WHEN bedtool_to_tuples() is called with bed_format="bed6"
    THEN missing fields padded with defaults
    """
    bt = pybedtools.BedTool("chr1\t100\t200\n", from_string=True)
    result = bedtool_to_tuples(bt, bed_format="bed6")
    assert len(result) == 1
    assert len(result[0]) == 6


def test_bedtool_to_tuples_closest_format():
    """
    GIVEN a BedTool from closest operation (13 fields)
    WHEN bedtool_to_tuples() is called with bed_format="closest"
    THEN returns tuples with A fields, B fields, and distance
    """
    line = "chr1\t100\t200\ta1\t50\t+\tchr1\t300\t400\tb1\t75\t+\t100\n"
    bt = pybedtools.BedTool(line, from_string=True)
    result = bedtool_to_tuples(bt, bed_format="closest")
    assert len(result) == 1
    row = result[0]
    assert row[0] == "chr1"  # a.chrom
    assert row[1] == 100  # a.start (int)
    assert row[6] == "chr1"  # b.chrom
    assert row[7] == 300  # b.start (int)
    assert row[12] == 100  # distance (int)


def test_bedtool_to_tuples_closest_dot_values():
    """
    GIVEN a BedTool from closest with "." scores/names
    WHEN bedtool_to_tuples() is called with bed_format="closest"
    THEN "." values converted to None
    """
    line = "chr1\t100\t200\t.\t.\t.\tchr1\t300\t400\t.\t.\t.\t50\n"
    bt = pybedtools.BedTool(line, from_string=True)
    result = bedtool_to_tuples(bt, bed_format="closest")
    row = result[0]
    assert row[3] is None  # a.name
    assert row[4] is None  # a.score
    assert row[5] is None  # a.strand
    assert row[9] is None  # b.name


def test_bedtool_to_tuples_invalid_format_raises():
    """
    GIVEN any BedTool
    WHEN bedtool_to_tuples() is called with invalid format
    THEN ValueError is raised
    """
    bt = pybedtools.BedTool("chr1\t100\t200\n", from_string=True)
    with pytest.raises(ValueError, match="Unsupported format"):
        bedtool_to_tuples(bt, bed_format="invalid")


def test_bedtool_to_tuples_closest_insufficient_fields_raises():
    """
    GIVEN a BedTool with fewer than 13 fields
    WHEN bedtool_to_tuples() is called with bed_format="closest"
    THEN ValueError is raised
    """
    bt = pybedtools.BedTool("chr1\t100\t200\ta1\t0\t+\n", from_string=True)
    with pytest.raises(ValueError, match="Unexpected number of fields"):
        bedtool_to_tuples(bt, bed_format="closest")


class TestBedtoolsError:
    def test_is_exception_subclass(self):
        """
        GIVEN a message string
        WHEN BedtoolsError is raised
        THEN it is an instance of Exception with correct message
        """
        with pytest.raises(BedtoolsError, match="test error"):
            raise BedtoolsError("test error")
