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


def test_create_bedtool_should_parse_bed3():
    """Test that create_bedtool constructs a BedTool from BED3 tuples.

    Given:
        A list of BED3 tuples
    When:
        create_bedtool() is called
    Then:
        It should return a BedTool with correct intervals
    """
    # Arrange / Act
    bt = create_bedtool([("chr1", 100, 200)])
    intervals = list(bt)

    # Assert
    assert len(intervals) == 1
    assert intervals[0].chrom == "chr1"
    assert intervals[0].start == 100
    assert intervals[0].end == 200


def test_create_bedtool_should_parse_bed6():
    """Test that create_bedtool constructs a BedTool from BED6 tuples.

    Given:
        A list of BED6 tuples
    When:
        create_bedtool() is called
    Then:
        It should return a BedTool with all 6 fields
    """
    # Arrange / Act
    bt = create_bedtool([("chr1", 100, 200, "a1", 50, "+")])
    intervals = list(bt)

    # Assert
    assert len(intervals) == 1
    assert intervals[0].fields == ["chr1", "100", "200", "a1", "50", "+"]


def test_create_bedtool_should_replace_none_with_defaults():
    """Test that create_bedtool substitutes defaults for None values.

    Given:
        BED6 tuples with None values
    When:
        create_bedtool() is called
    Then:
        It should replace None values with defaults
    """
    # Arrange / Act
    bt = create_bedtool([("chr1", 100, 200, None, None, None)])
    fields = list(bt)[0].fields

    # Assert
    assert fields[3] == "."  # name
    assert fields[4] == "0"  # score
    assert fields[5] == "."  # strand


def test_create_bedtool_should_raise_when_tuple_length_invalid():
    """Test that create_bedtool rejects tuples with wrong arity.

    Given:
        A tuple with invalid length
    When:
        create_bedtool() is called
    Then:
        It should raise ValueError
    """
    # Arrange / Act / Assert
    with pytest.raises(ValueError, match="Invalid interval format"):
        create_bedtool([("chr1", 100)])


def test_create_bedtool_should_accept_multiple_intervals():
    """Test that create_bedtool handles multiple intervals across chromosomes.

    Given:
        Multiple intervals across chromosomes
    When:
        create_bedtool() is called
    Then:
        It should return a BedTool containing all intervals
    """
    # Arrange / Act
    bt = create_bedtool(
        [
            ("chr1", 100, 200, "a", 0, "+"),
            ("chr2", 300, 400, "b", 0, "-"),
        ]
    )
    intervals = list(bt)

    # Assert
    assert len(intervals) == 2


def test_intersect_should_return_overlapping_intervals():
    """Test that intersect returns A intervals overlapping B.

    Given:
        Two sets of overlapping intervals
    When:
        intersect() is called
    Then:
        It should return intervals from A that overlap B
    """
    # Arrange
    a = [("chr1", 100, 200, "a1", 100, "+")]
    b = [("chr1", 150, 250, "b1", 100, "+")]

    # Act
    result = intersect(a, b)

    # Assert
    assert len(result) == 1
    assert result[0][0] == "chr1"


def test_intersect_should_return_empty_when_no_overlap():
    """Test that intersect returns no rows when intervals disjoint.

    Given:
        Non-overlapping intervals
    When:
        intersect() is called
    Then:
        It should return an empty list
    """
    # Arrange
    a = [("chr1", 100, 200, "a1", 100, "+")]
    b = [("chr1", 300, 400, "b1", 100, "+")]

    # Act
    result = intersect(a, b)

    # Assert
    assert result == []


def test_intersect_should_filter_same_strand_only_when_strand_mode_same():
    """Test that intersect in same-strand mode keeps only same-strand hits.

    Given:
        Intervals on same and opposite strands
    When:
        intersect() is called with strand_mode="same"
    Then:
        It should return only same-strand overlaps
    """
    # Arrange
    a = [
        ("chr1", 100, 200, "a1", 0, "+"),
        ("chr1", 100, 200, "a2", 0, "-"),
    ]
    b = [("chr1", 150, 250, "b1", 0, "+")]

    # Act
    result = intersect(a, b, strand_mode="same")
    names = [r[3] for r in result]

    # Assert
    assert "a1" in names
    assert "a2" not in names


def test_intersect_should_filter_opposite_strand_only_when_strand_mode_opposite():
    """Test that intersect in opposite-strand mode keeps only opposite-strand hits.

    Given:
        Intervals on same and opposite strands
    When:
        intersect() is called with strand_mode="opposite"
    Then:
        It should return only opposite-strand overlaps
    """
    # Arrange
    a = [
        ("chr1", 100, 200, "a1", 0, "+"),
        ("chr1", 100, 200, "a2", 0, "-"),
    ]
    b = [("chr1", 150, 250, "b1", 0, "+")]

    # Act
    result = intersect(a, b, strand_mode="opposite")
    names = [r[3] for r in result]

    # Assert
    assert "a2" in names
    assert "a1" not in names


def test_intersect_should_ignore_strand_when_strand_mode_none():
    """Test that intersect ignores strand when strand_mode is None.

    Given:
        Overlapping intervals on different strands
    When:
        intersect() is called with strand_mode=None
    Then:
        It should return all overlaps regardless of strand
    """
    # Arrange
    a = [("chr1", 100, 200, "a1", 0, "+")]
    b = [("chr1", 150, 250, "b1", 0, "-")]

    # Act
    result = intersect(a, b)

    # Assert
    assert len(result) == 1


def test_merge_should_combine_overlapping_intervals():
    """Test that merge collapses overlapping intervals into one.

    Given:
        Overlapping intervals
    When:
        merge() is called
    Then:
        It should return merged BED3 intervals
    """
    # Arrange
    intervals = [
        ("chr1", 100, 200, "i1", 0, "+"),
        ("chr1", 150, 250, "i2", 0, "+"),
    ]

    # Act
    result = merge(intervals)

    # Assert
    assert len(result) == 1
    assert result[0] == ("chr1", 100, 250)


def test_merge_should_preserve_separated_intervals():
    """Test that merge keeps non-overlapping intervals separate.

    Given:
        Separated intervals
    When:
        merge() is called
    Then:
        It should return each interval separately as BED3
    """
    # Arrange
    intervals = [
        ("chr1", 100, 200, "i1", 0, "+"),
        ("chr1", 300, 400, "i2", 0, "+"),
    ]

    # Act
    result = merge(intervals)

    # Assert
    assert len(result) == 2


def test_merge_should_merge_per_strand_when_strand_mode_same():
    """Test that merge segregates intervals by strand in same-strand mode.

    Given:
        Overlapping intervals on different strands
    When:
        merge() is called with strand_mode="same"
    Then:
        It should merge per-strand separately
    """
    # Arrange
    intervals = [
        ("chr1", 100, 200, "i1", 0, "+"),
        ("chr1", 150, 250, "i2", 0, "+"),
        ("chr1", 120, 220, "i3", 0, "-"),
    ]

    # Act
    result = merge(intervals, strand_mode="same")

    # Assert
    # Should have 2: one merged + strand, one - strand
    assert len(result) == 2


def test_merge_should_combine_adjacent_intervals():
    """Test that merge joins intervals where one ends at the next's start.

    Given:
        Adjacent intervals (end == start of next)
    When:
        merge() is called
    Then:
        It should merge adjacent intervals
    """
    # Arrange
    intervals = [
        ("chr1", 100, 200, "i1", 0, "+"),
        ("chr1", 200, 300, "i2", 0, "+"),
    ]

    # Act
    result = merge(intervals)

    # Assert
    assert len(result) == 1
    assert result[0] == ("chr1", 100, 300)


def test_closest_should_pair_a_with_nearest_b_and_distance():
    """Test that closest pairs each A interval with the nearest B and a distance.

    Given:
        Non-overlapping intervals
    When:
        closest() is called
    Then:
        It should return each A paired with nearest B plus distance
    """
    # Arrange
    a = [("chr1", 100, 200, "a1", 100, "+")]
    b = [("chr1", 300, 400, "b1", 100, "+")]

    # Act
    result = closest(a, b)

    # Assert
    assert len(result) == 1
    # bedtools >= 2.31 reports N+1 for an N-base half-open gap between
    # intervals (here 300 - 200 = 100, so expected distance is 101).
    # The project pins bedtools >= 2.31.0 via pixi.
    assert result[0][-1] == 101


def test_closest_should_match_per_chromosome():
    """Test that closest restricts neighbor search to the same chromosome.

    Given:
        Intervals on different chromosomes
    When:
        closest() is called
    Then:
        It should find the nearest per-chromosome
    """
    # Arrange
    a = [
        ("chr1", 100, 200, "a1", 0, "+"),
        ("chr2", 100, 200, "a2", 0, "+"),
    ]
    b = [
        ("chr1", 300, 400, "b1", 0, "+"),
        ("chr2", 500, 600, "b2", 0, "+"),
    ]

    # Act
    result = closest(a, b)

    # Assert
    assert len(result) == 2
    # Each A should match B on same chromosome
    for row in result:
        assert row[0] == row[6]  # a.chrom == b.chrom


def test_closest_should_return_nearest_same_strand_when_strand_mode_same():
    """Test that closest in same-strand mode picks the nearest same-strand B.

    Given:
        Intervals with mixed strands
    When:
        closest() is called with strand_mode="same"
    Then:
        It should return the nearest same-strand interval
    """
    # Arrange
    a = [("chr1", 100, 200, "a1", 0, "+")]
    b = [
        ("chr1", 220, 240, "b_opp", 0, "-"),  # closer but opposite
        ("chr1", 300, 400, "b_same", 0, "+"),  # farther but same
    ]

    # Act
    result = closest(a, b, strand_mode="same")

    # Assert
    assert len(result) == 1
    assert result[0][9] == "b_same"


def test_closest_should_return_k_neighbors():
    """Test that closest returns up to k nearest neighbors when k > 1.

    Given:
        One query and three database intervals
    When:
        closest() is called with k=3
    Then:
        It should return up to 3 nearest
    """
    # Arrange
    a = [("chr1", 200, 300, "a1", 0, "+")]
    b = [
        ("chr1", 100, 150, "b1", 0, "+"),
        ("chr1", 350, 400, "b2", 0, "+"),
        ("chr1", 500, 600, "b3", 0, "+"),
    ]

    # Act
    result = closest(a, b, k=3)

    # Assert
    # bedtools 2.31 with -t first collapses tied-distance candidates
    # (b1 and b2 are both distance 51 from a1), so k=3 returns 2 rows
    # for this specific fixture rather than 3.
    assert len(result) == 2


def test_bedtool_to_tuples_should_parse_bed3():
    """Test that bedtool_to_tuples converts BED3 intervals to 3-tuples.

    Given:
        A BedTool with BED3 intervals
    When:
        bedtool_to_tuples() is called with bed_format="bed3"
    Then:
        It should return a list of (chrom, start, end) tuples with int positions
    """
    # Arrange
    bt = pybedtools.BedTool("chr1\t100\t200\n", from_string=True)

    # Act
    result = bedtool_to_tuples(bt, bed_format="bed3")

    # Assert
    assert result == [("chr1", 100, 200)]


def test_bedtool_to_tuples_should_parse_bed6():
    """Test that bedtool_to_tuples converts BED6 intervals to 6-tuples.

    Given:
        A BedTool with BED6 intervals
    When:
        bedtool_to_tuples() is called with bed_format="bed6"
    Then:
        It should return a list of 6-tuples with correct types
    """
    # Arrange
    bt = pybedtools.BedTool("chr1\t100\t200\tgene1\t500\t+\n", from_string=True)

    # Act
    result = bedtool_to_tuples(bt, bed_format="bed6")

    # Assert
    assert result == [("chr1", 100, 200, "gene1", 500, "+")]


def test_bedtool_to_tuples_should_convert_dot_to_none_for_bed6():
    """Test that bedtool_to_tuples maps "." placeholders to None in BED6.

    Given:
        A BedTool with "." for name and strand
    When:
        bedtool_to_tuples() is called with bed_format="bed6"
    Then:
        It should convert "." values to None
    """
    # Arrange
    bt = pybedtools.BedTool("chr1\t100\t200\t.\t0\t.\n", from_string=True)

    # Act
    result = bedtool_to_tuples(bt, bed_format="bed6")

    # Assert
    assert result[0][3] is None  # name
    assert result[0][5] is None  # strand


def test_bedtool_to_tuples_should_pad_missing_bed6_fields():
    """Test that bedtool_to_tuples pads short rows to 6 fields.

    Given:
        A BedTool with fewer than 6 fields
    When:
        bedtool_to_tuples() is called with bed_format="bed6"
    Then:
        It should pad missing fields with defaults
    """
    # Arrange
    bt = pybedtools.BedTool("chr1\t100\t200\n", from_string=True)

    # Act
    result = bedtool_to_tuples(bt, bed_format="bed6")

    # Assert
    assert len(result) == 1
    assert len(result[0]) == 6


def test_bedtool_to_tuples_should_parse_closest_format():
    """Test that bedtool_to_tuples parses the 13-field closest format.

    Given:
        A BedTool from closest operation (13 fields)
    When:
        bedtool_to_tuples() is called with bed_format="closest"
    Then:
        It should return tuples with A fields, B fields, and distance
    """
    # Arrange
    line = "chr1\t100\t200\ta1\t50\t+\tchr1\t300\t400\tb1\t75\t+\t100\n"
    bt = pybedtools.BedTool(line, from_string=True)

    # Act
    result = bedtool_to_tuples(bt, bed_format="closest")

    # Assert
    assert len(result) == 1
    row = result[0]
    assert row[0] == "chr1"  # a.chrom
    assert row[1] == 100  # a.start (int)
    assert row[6] == "chr1"  # b.chrom
    assert row[7] == 300  # b.start (int)
    assert row[12] == 100  # distance (int)


def test_bedtool_to_tuples_should_convert_dot_to_none_for_closest():
    """Test that bedtool_to_tuples maps "." placeholders to None in closest rows.

    Given:
        A BedTool from closest with "." scores/names
    When:
        bedtool_to_tuples() is called with bed_format="closest"
    Then:
        It should convert "." values to None
    """
    # Arrange
    line = "chr1\t100\t200\t.\t.\t.\tchr1\t300\t400\t.\t.\t.\t50\n"
    bt = pybedtools.BedTool(line, from_string=True)

    # Act
    result = bedtool_to_tuples(bt, bed_format="closest")

    # Assert
    row = result[0]
    assert row[3] is None  # a.name
    assert row[4] is None  # a.score
    assert row[5] is None  # a.strand
    assert row[9] is None  # b.name


def test_bedtool_to_tuples_should_raise_when_format_invalid():
    """Test that bedtool_to_tuples rejects unknown bed_format values.

    Given:
        Any BedTool
    When:
        bedtool_to_tuples() is called with invalid format
    Then:
        It should raise ValueError
    """
    # Arrange
    bt = pybedtools.BedTool("chr1\t100\t200\n", from_string=True)

    # Act / Assert
    with pytest.raises(ValueError, match="Unsupported format"):
        bedtool_to_tuples(bt, bed_format="invalid")


def test_bedtool_to_tuples_should_raise_when_closest_fields_insufficient():
    """Test that bedtool_to_tuples rejects closest rows with too few fields.

    Given:
        A BedTool with fewer than 13 fields
    When:
        bedtool_to_tuples() is called with bed_format="closest"
    Then:
        It should raise ValueError
    """
    # Arrange
    bt = pybedtools.BedTool("chr1\t100\t200\ta1\t0\t+\n", from_string=True)

    # Act / Assert
    with pytest.raises(ValueError, match="Unexpected number of fields"):
        bedtool_to_tuples(bt, bed_format="closest")


class TestBedtoolsError:
    def test___init___should_create_exception_with_message(self):
        """Test that BedtoolsError behaves as an Exception carrying its message.

        Given:
            A message string
        When:
            BedtoolsError is raised
        Then:
            It should be an instance of Exception with the correct message
        """
        # Arrange / Act / Assert
        with pytest.raises(BedtoolsError, match="test error"):
            raise BedtoolsError("test error")
