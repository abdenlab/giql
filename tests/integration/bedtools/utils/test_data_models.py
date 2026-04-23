"""Unit tests for bedtools integration test data models."""

import pytest
from hypothesis import given
from hypothesis import strategies as st

from .data_models import ComparisonResult
from .data_models import GenomicInterval

pytestmark = pytest.mark.integration


class TestGenomicInterval:
    def test_basic_instantiation(self):
        """
        GIVEN valid chrom, start, end values
        WHEN GenomicInterval is instantiated
        THEN object is created with correct attributes
        """
        gi = GenomicInterval("chr1", 100, 200)
        assert gi.chrom == "chr1"
        assert gi.start == 100
        assert gi.end == 200
        assert gi.name is None
        assert gi.score is None
        assert gi.strand is None

    def test_full_instantiation(self):
        """
        GIVEN all fields provided
        WHEN GenomicInterval is instantiated
        THEN all attributes are set correctly
        """
        gi = GenomicInterval("chrX", 500, 1000, "gene1", 800, "+")
        assert gi.chrom == "chrX"
        assert gi.start == 500
        assert gi.end == 1000
        assert gi.name == "gene1"
        assert gi.score == 800
        assert gi.strand == "+"

    def test_start_equals_end_raises(self):
        """
        GIVEN start equals end
        WHEN GenomicInterval is instantiated
        THEN ValueError is raised
        """
        with pytest.raises(ValueError, match="start .* >= end"):
            GenomicInterval("chr1", 200, 200)

    def test_start_greater_than_end_raises(self):
        """
        GIVEN start > end
        WHEN GenomicInterval is instantiated
        THEN ValueError is raised
        """
        with pytest.raises(ValueError, match="start .* >= end"):
            GenomicInterval("chr1", 300, 200)

    def test_negative_start_raises(self):
        """
        GIVEN start < 0
        WHEN GenomicInterval is instantiated
        THEN ValueError is raised
        """
        with pytest.raises(ValueError, match="start .* < 0"):
            GenomicInterval("chr1", -1, 200)

    def test_invalid_strand_raises(self):
        """
        GIVEN an invalid strand value
        WHEN GenomicInterval is instantiated
        THEN ValueError is raised
        """
        with pytest.raises(ValueError, match="Invalid strand"):
            GenomicInterval("chr1", 100, 200, strand="X")

    def test_score_below_range_raises(self):
        """
        GIVEN score < 0
        WHEN GenomicInterval is instantiated
        THEN ValueError is raised
        """
        with pytest.raises(ValueError, match="Invalid score"):
            GenomicInterval("chr1", 100, 200, score=-1)

    def test_score_above_range_raises(self):
        """
        GIVEN score > 1000
        WHEN GenomicInterval is instantiated
        THEN ValueError is raised
        """
        with pytest.raises(ValueError, match="Invalid score"):
            GenomicInterval("chr1", 100, 200, score=1001)

    @pytest.mark.parametrize("strand", ["+", "-", "."])
    def test_valid_strand_values(self, strand):
        """
        GIVEN a valid strand value
        WHEN GenomicInterval is instantiated
        THEN object is created successfully
        """
        gi = GenomicInterval("chr1", 100, 200, strand=strand)
        assert gi.strand == strand

    def test_score_boundary_zero(self):
        """
        GIVEN score = 0
        WHEN GenomicInterval is instantiated
        THEN object is created successfully
        """
        gi = GenomicInterval("chr1", 100, 200, score=0)
        assert gi.score == 0

    def test_score_boundary_thousand(self):
        """
        GIVEN score = 1000
        WHEN GenomicInterval is instantiated
        THEN object is created successfully
        """
        gi = GenomicInterval("chr1", 100, 200, score=1000)
        assert gi.score == 1000

    def test_to_tuple(self):
        """
        GIVEN a GenomicInterval with all fields
        WHEN to_tuple() is called
        THEN returns 6-element tuple with all field values
        """
        gi = GenomicInterval("chr1", 100, 200, "a1", 500, "+")
        assert gi.to_tuple() == ("chr1", 100, 200, "a1", 500, "+")

    def test_to_tuple_with_nones(self):
        """
        GIVEN a GenomicInterval with optional fields as None
        WHEN to_tuple() is called
        THEN tuple contains None for optional fields
        """
        gi = GenomicInterval("chr1", 100, 200)
        assert gi.to_tuple() == ("chr1", 100, 200, None, None, None)

    @given(
        chrom=st.sampled_from(["chr1", "chr2", "chrX", "chrM"]),
        start=st.integers(min_value=0, max_value=999_999),
        size=st.integers(min_value=1, max_value=10_000),
        strand=st.sampled_from(["+", "-", "."]),
        score=st.integers(min_value=0, max_value=1000),
    )
    def test_to_tuple_roundtrip(self, chrom, start, size, strand, score):
        """
        GIVEN any valid GenomicInterval
        WHEN to_tuple() is called
        THEN the tuple can be used to reconstruct the interval's key fields
        """
        end = start + size
        gi = GenomicInterval(chrom, start, end, "name", score, strand)
        t = gi.to_tuple()
        assert t == (chrom, start, end, "name", score, strand)


class TestComparisonResult:
    def test_matching_result(self):
        """
        GIVEN match=True with equal row counts
        WHEN ComparisonResult is instantiated
        THEN attributes are set correctly
        """
        cr = ComparisonResult(match=True, giql_row_count=5, bedtools_row_count=5)
        assert cr.match is True
        assert cr.giql_row_count == 5
        assert cr.bedtools_row_count == 5
        assert cr.differences == []

    def test_mismatching_result(self):
        """
        GIVEN match=False with differences
        WHEN ComparisonResult is instantiated
        THEN attributes are set correctly
        """
        diffs = ["Row 0: mismatch"]
        cr = ComparisonResult(
            match=False,
            giql_row_count=3,
            bedtools_row_count=4,
            differences=diffs,
        )
        assert cr.match is False
        assert cr.differences == diffs

    def test_bool_true(self):
        """
        GIVEN a matching ComparisonResult
        WHEN used in boolean context
        THEN evaluates to True
        """
        cr = ComparisonResult(match=True, giql_row_count=1, bedtools_row_count=1)
        assert cr

    def test_bool_false(self):
        """
        GIVEN a non-matching ComparisonResult
        WHEN used in boolean context
        THEN evaluates to False
        """
        cr = ComparisonResult(match=False, giql_row_count=1, bedtools_row_count=2)
        assert not cr

    def test_failure_message_match(self):
        """
        GIVEN a matching ComparisonResult
        WHEN failure_message() is called
        THEN returns success message
        """
        cr = ComparisonResult(match=True, giql_row_count=1, bedtools_row_count=1)
        assert "match" in cr.failure_message().lower()

    def test_failure_message_mismatch(self):
        """
        GIVEN a non-matching ComparisonResult with differences
        WHEN failure_message() is called
        THEN returns formatted message with row counts and differences
        """
        cr = ComparisonResult(
            match=False,
            giql_row_count=3,
            bedtools_row_count=5,
            differences=["Row 0: val mismatch", "Row 1: missing"],
        )
        msg = cr.failure_message()
        assert "3" in msg
        assert "5" in msg
        assert "Row 0: val mismatch" in msg
        assert "Row 1: missing" in msg

    def test_failure_message_truncates_at_ten(self):
        """
        GIVEN a ComparisonResult with more than 10 differences
        WHEN failure_message() is called
        THEN only first 10 are shown with a count of remaining
        """
        diffs = [f"diff_{i}" for i in range(15)]
        cr = ComparisonResult(
            match=False,
            giql_row_count=0,
            bedtools_row_count=15,
            differences=diffs,
        )
        msg = cr.failure_message()
        assert "diff_9" in msg
        assert "diff_10" not in msg
        assert "5 more" in msg

    def test_default_metadata(self):
        """
        GIVEN no comparison_metadata provided
        WHEN ComparisonResult is instantiated
        THEN metadata defaults to empty dict
        """
        cr = ComparisonResult(match=True, giql_row_count=0, bedtools_row_count=0)
        assert cr.comparison_metadata == {}
