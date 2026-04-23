"""Unit tests for bedtools integration test data models."""

import pytest
from hypothesis import given
from hypothesis import settings
from hypothesis import strategies as st

from .data_models import ComparisonResult
from .data_models import GenomicInterval

pytestmark = pytest.mark.integration


class TestGenomicInterval:
    def test___init___should_succeed_when_minimal_args_supplied(self):
        """Test that minimal instantiation populates required fields and defaults.

        Given:
            Valid chrom, start, end values
        When:
            GenomicInterval is instantiated
        Then:
            It should create an object with correct attributes and None defaults
        """
        # Arrange / Act
        gi = GenomicInterval("chr1", 100, 200)

        # Assert
        assert gi.chrom == "chr1"
        assert gi.start == 100
        assert gi.end == 200
        assert gi.name is None
        assert gi.score is None
        assert gi.strand is None

    def test___init___should_populate_optional_fields_when_supplied(self):
        """Test that all fields are set when provided to the constructor.

        Given:
            All fields provided
        When:
            GenomicInterval is instantiated
        Then:
            It should set all attributes correctly
        """
        # Arrange / Act
        gi = GenomicInterval("chrX", 500, 1000, "gene1", 800, "+")

        # Assert
        assert gi.chrom == "chrX"
        assert gi.start == 500
        assert gi.end == 1000
        assert gi.name == "gene1"
        assert gi.score == 800
        assert gi.strand == "+"

    def test___post_init___should_raise_when_start_equals_end(self):
        """Test that a zero-length interval is rejected.

        Given:
            start equals end
        When:
            GenomicInterval is instantiated
        Then:
            It should raise ValueError
        """
        # Arrange / Act / Assert
        with pytest.raises(ValueError, match="start .* >= end"):
            GenomicInterval("chr1", 200, 200)

    def test___post_init___should_raise_when_start_greater_than_end(self):
        """Test that an inverted interval is rejected.

        Given:
            start > end
        When:
            GenomicInterval is instantiated
        Then:
            It should raise ValueError
        """
        # Arrange / Act / Assert
        with pytest.raises(ValueError, match="start .* >= end"):
            GenomicInterval("chr1", 300, 200)

    def test___post_init___should_raise_when_start_is_negative(self):
        """Test that a negative start coordinate is rejected.

        Given:
            start < 0
        When:
            GenomicInterval is instantiated
        Then:
            It should raise ValueError
        """
        # Arrange / Act / Assert
        with pytest.raises(ValueError, match="start .* < 0"):
            GenomicInterval("chr1", -1, 200)

    def test___post_init___should_raise_when_strand_is_invalid(self):
        """Test that an invalid strand value is rejected.

        Given:
            An invalid strand value
        When:
            GenomicInterval is instantiated
        Then:
            It should raise ValueError
        """
        # Arrange / Act / Assert
        with pytest.raises(ValueError, match="Invalid strand"):
            GenomicInterval("chr1", 100, 200, strand="X")

    def test___post_init___should_raise_when_score_below_range(self):
        """Test that a score below the BED range is rejected.

        Given:
            score < 0
        When:
            GenomicInterval is instantiated
        Then:
            It should raise ValueError
        """
        # Arrange / Act / Assert
        with pytest.raises(ValueError, match="Invalid score"):
            GenomicInterval("chr1", 100, 200, score=-1)

    def test___post_init___should_raise_when_score_above_range(self):
        """Test that a score above the BED range is rejected.

        Given:
            score > 1000
        When:
            GenomicInterval is instantiated
        Then:
            It should raise ValueError
        """
        # Arrange / Act / Assert
        with pytest.raises(ValueError, match="Invalid score"):
            GenomicInterval("chr1", 100, 200, score=1001)

    @pytest.mark.parametrize("strand", ["+", "-", "."])
    def test___post_init___should_accept_when_strand_is_valid(self, strand):
        """Test that each allowed strand value is accepted.

        Given:
            A valid strand value
        When:
            GenomicInterval is instantiated
        Then:
            It should create the object successfully
        """
        # Arrange / Act
        gi = GenomicInterval("chr1", 100, 200, strand=strand)

        # Assert
        assert gi.strand == strand

    def test___post_init___should_accept_when_score_is_zero(self):
        """Test that the lower boundary score is accepted.

        Given:
            score = 0
        When:
            GenomicInterval is instantiated
        Then:
            It should create the object successfully
        """
        # Arrange / Act
        gi = GenomicInterval("chr1", 100, 200, score=0)

        # Assert
        assert gi.score == 0

    def test___post_init___should_accept_when_score_is_thousand(self):
        """Test that the upper boundary score is accepted.

        Given:
            score = 1000
        When:
            GenomicInterval is instantiated
        Then:
            It should create the object successfully
        """
        # Arrange / Act
        gi = GenomicInterval("chr1", 100, 200, score=1000)

        # Assert
        assert gi.score == 1000

    def test_to_tuple_should_return_all_fields_when_fully_populated(self):
        """Test that to_tuple returns every field in order.

        Given:
            A GenomicInterval with all fields
        When:
            to_tuple() is called
        Then:
            It should return a 6-element tuple with all field values
        """
        # Arrange
        gi = GenomicInterval("chr1", 100, 200, "a1", 500, "+")

        # Act
        result = gi.to_tuple()

        # Assert
        assert result == ("chr1", 100, 200, "a1", 500, "+")

    def test_to_tuple_should_include_none_when_optional_fields_missing(self):
        """Test that to_tuple preserves None for unset optional fields.

        Given:
            A GenomicInterval with optional fields as None
        When:
            to_tuple() is called
        Then:
            It should return a tuple containing None for optional fields
        """
        # Arrange
        gi = GenomicInterval("chr1", 100, 200)

        # Act
        result = gi.to_tuple()

        # Assert
        assert result == ("chr1", 100, 200, None, None, None)

    @settings(max_examples=50)
    @given(
        chrom=st.sampled_from(["chr1", "chr2", "chrX", "chrM"]),
        start=st.integers(min_value=0, max_value=999_999),
        size=st.integers(min_value=1, max_value=10_000),
        strand=st.sampled_from(["+", "-", "."]),
        score=st.integers(min_value=0, max_value=1000),
    )
    def test_to_tuple_should_roundtrip_when_any_valid_interval(
        self, chrom, start, size, strand, score
    ):
        """Test that to_tuple reflects the exact constructor inputs.

        Given:
            Any valid GenomicInterval
        When:
            to_tuple() is called
        Then:
            It should return a tuple that matches the interval's key fields
        """
        # Arrange
        end = start + size
        gi = GenomicInterval(chrom, start, end, "name", score, strand)

        # Act
        t = gi.to_tuple()

        # Assert
        assert t == (chrom, start, end, "name", score, strand)


class TestComparisonResult:
    def test___init___should_populate_attributes_when_match_is_true(self):
        """Test that a matching result stores its fields correctly.

        Given:
            match=True with equal row counts
        When:
            ComparisonResult is instantiated
        Then:
            It should set attributes correctly with an empty differences list
        """
        # Arrange / Act
        cr = ComparisonResult(match=True, giql_row_count=5, bedtools_row_count=5)

        # Assert
        assert cr.match is True
        assert cr.giql_row_count == 5
        assert cr.bedtools_row_count == 5
        assert cr.differences == []

    def test___init___should_populate_attributes_when_match_is_false(self):
        """Test that a mismatching result stores its differences.

        Given:
            match=False with differences
        When:
            ComparisonResult is instantiated
        Then:
            It should set attributes correctly including the differences list
        """
        # Arrange
        diffs = ["Row 0: mismatch"]

        # Act
        cr = ComparisonResult(
            match=False,
            giql_row_count=3,
            bedtools_row_count=4,
            differences=diffs,
        )

        # Assert
        assert cr.match is False
        assert cr.differences == diffs

    def test___bool___should_return_true_when_match_is_true(self):
        """Test truthiness of a matching result.

        Given:
            A matching ComparisonResult
        When:
            Used in a boolean context
        Then:
            It should evaluate to True
        """
        # Arrange
        cr = ComparisonResult(match=True, giql_row_count=1, bedtools_row_count=1)

        # Act / Assert
        assert cr

    def test___bool___should_return_false_when_match_is_false(self):
        """Test falsiness of a non-matching result.

        Given:
            A non-matching ComparisonResult
        When:
            Used in a boolean context
        Then:
            It should evaluate to False
        """
        # Arrange
        cr = ComparisonResult(match=False, giql_row_count=1, bedtools_row_count=2)

        # Act / Assert
        assert not cr

    def test_failure_message_should_return_success_when_match_is_true(self):
        """Test the message for a matching result.

        Given:
            A matching ComparisonResult
        When:
            failure_message() is called
        Then:
            It should return a success message
        """
        # Arrange
        cr = ComparisonResult(match=True, giql_row_count=1, bedtools_row_count=1)

        # Act
        msg = cr.failure_message()

        # Assert
        assert "match" in msg.lower()

    def test_failure_message_should_include_counts_and_diffs_when_mismatch(self):
        """Test the message formatting for a mismatching result.

        Given:
            A non-matching ComparisonResult with differences
        When:
            failure_message() is called
        Then:
            It should return a formatted message with row counts and differences
        """
        # Arrange
        cr = ComparisonResult(
            match=False,
            giql_row_count=3,
            bedtools_row_count=5,
            differences=["Row 0: val mismatch", "Row 1: missing"],
        )

        # Act
        msg = cr.failure_message()

        # Assert
        assert "3" in msg
        assert "5" in msg
        assert "Row 0: val mismatch" in msg
        assert "Row 1: missing" in msg

    def test_failure_message_should_truncate_when_over_ten_differences(self):
        """Test that the message truncates the differences list at ten.

        Given:
            A ComparisonResult with more than 10 differences
        When:
            failure_message() is called
        Then:
            It should show only the first 10 with a count of the remainder
        """
        # Arrange
        diffs = [f"diff_{i}" for i in range(15)]
        cr = ComparisonResult(
            match=False,
            giql_row_count=0,
            bedtools_row_count=15,
            differences=diffs,
        )

        # Act
        msg = cr.failure_message()

        # Assert
        assert "diff_9" in msg
        assert "diff_10" not in msg
        assert "5 more" in msg

    def test___init___should_default_metadata_when_not_supplied(self):
        """Test that comparison_metadata defaults to an empty dict.

        Given:
            No comparison_metadata provided
        When:
            ComparisonResult is instantiated
        Then:
            It should default metadata to an empty dict
        """
        # Arrange / Act
        cr = ComparisonResult(match=True, giql_row_count=0, bedtools_row_count=0)

        # Assert
        assert cr.comparison_metadata == {}
