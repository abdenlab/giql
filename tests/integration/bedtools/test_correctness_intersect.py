"""Extended correctness tests for GIQL INTERSECTS operator vs bedtools intersect.

These tests cover boundary cases, scale, and edge scenarios beyond the basic
tests in test_intersect.py, ensuring comprehensive GIQL/bedtools equivalence.
"""

import pytest

from giql import transpile

from .utils.bedtools_wrapper import intersect
from .utils.comparison import compare_results
from .utils.data_models import GenomicInterval
from .utils.duckdb_loader import load_intervals
from .utils.random_intervals import generate_random_intervals

pytestmark = pytest.mark.integration


def _run_intersect_comparison(
    duckdb_connection,
    intervals_a,
    intervals_b,
    strand_filter="",
):
    """Run GIQL INTERSECTS and bedtools intersect, return ComparisonResult."""
    load_intervals(
        duckdb_connection,
        "intervals_a",
        [i.to_tuple() for i in intervals_a],
    )
    load_intervals(
        duckdb_connection,
        "intervals_b",
        [i.to_tuple() for i in intervals_b],
    )

    strand_mode = None
    if "a.strand = b.strand" in strand_filter:
        strand_mode = "same"
    elif "a.strand != b.strand" in strand_filter:
        strand_mode = "opposite"

    bedtools_result = intersect(
        [i.to_tuple() for i in intervals_a],
        [i.to_tuple() for i in intervals_b],
        strand_mode=strand_mode,
    )

    where_clause = "WHERE a.interval INTERSECTS b.interval"
    if strand_filter:
        where_clause += f" AND {strand_filter}"

    sql = transpile(
        f"""
        SELECT DISTINCT a.*
        FROM intervals_a a, intervals_b b
        {where_clause}
        """,
        tables=["intervals_a", "intervals_b"],
    )
    giql_result = duckdb_connection.execute(sql).fetchall()

    return compare_results(giql_result, bedtools_result)


def test_intersects_should_match_bedtools_when_overlap_is_one_bp(duckdb_connection):
    """Test INTERSECTS matches bedtools for a minimal 1bp overlap.

    Given:
        Two intervals that overlap by exactly one base pair
    When:
        GIQL INTERSECTS is compared to bedtools intersect
    Then:
        It should detect the 1bp overlap identically to bedtools
    """
    # Arrange
    a = [GenomicInterval("chr1", 100, 200, "a1", 0, "+")]
    b = [GenomicInterval("chr1", 199, 300, "b1", 0, "+")]

    # Act
    comparison = _run_intersect_comparison(duckdb_connection, a, b)

    # Assert
    assert comparison.match, comparison.failure_message()
    assert comparison.giql_row_count == 1


def test_intersects_should_match_bedtools_when_a_contains_b(duckdb_connection):
    """Test INTERSECTS matches bedtools when A fully contains B.

    Given:
        Interval A that fully contains interval B
    When:
        GIQL INTERSECTS is compared to bedtools intersect
    Then:
        It should report A as intersecting B
    """
    # Arrange
    a = [GenomicInterval("chr1", 100, 500, "a1", 0, "+")]
    b = [GenomicInterval("chr1", 200, 300, "b1", 0, "+")]

    # Act
    comparison = _run_intersect_comparison(duckdb_connection, a, b)

    # Assert
    assert comparison.match, comparison.failure_message()
    assert comparison.giql_row_count == 1


def test_intersects_should_match_bedtools_when_b_contains_a(duckdb_connection):
    """Test INTERSECTS matches bedtools when B fully contains A.

    Given:
        Interval B that fully contains interval A
    When:
        GIQL INTERSECTS is compared to bedtools intersect
    Then:
        It should report A as intersecting B
    """
    # Arrange
    a = [GenomicInterval("chr1", 200, 300, "a1", 0, "+")]
    b = [GenomicInterval("chr1", 100, 500, "b1", 0, "+")]

    # Act
    comparison = _run_intersect_comparison(duckdb_connection, a, b)

    # Assert
    assert comparison.match, comparison.failure_message()
    assert comparison.giql_row_count == 1


def test_intersects_should_deduplicate_when_a_overlaps_multiple_b(duckdb_connection):
    """Test INTERSECTS with DISTINCT matches bedtools -u deduplication.

    Given:
        One interval in A that overlaps several intervals in B
    When:
        GIQL INTERSECTS with DISTINCT is compared to bedtools intersect -u
    Then:
        It should report the A interval exactly once
    """
    # Arrange
    a = [GenomicInterval("chr1", 100, 300, "a1", 0, "+")]
    b = [
        GenomicInterval("chr1", 150, 200, "b1", 0, "+"),
        GenomicInterval("chr1", 200, 250, "b2", 0, "+"),
        GenomicInterval("chr1", 250, 350, "b3", 0, "+"),
    ]

    # Act
    comparison = _run_intersect_comparison(duckdb_connection, a, b)

    # Assert
    assert comparison.match, comparison.failure_message()
    assert comparison.giql_row_count == 1


def test_intersects_should_match_bedtools_when_chromosome_names_are_non_standard(
    duckdb_connection,
):
    """Test INTERSECTS matches bedtools on non-standard chromosome names.

    Given:
        Intervals on non-standard chromosome names like chrM and chrUn
    When:
        GIQL INTERSECTS is compared to bedtools intersect
    Then:
        It should match bedtools regardless of chromosome naming
    """
    # Arrange
    a = [
        GenomicInterval("chrM", 100, 200, "a1", 0, "+"),
        GenomicInterval("chrUn", 100, 200, "a2", 0, "+"),
    ]
    b = [
        GenomicInterval("chrM", 150, 250, "b1", 0, "+"),
        GenomicInterval("chrUn", 150, 250, "b2", 0, "+"),
    ]

    # Act
    comparison = _run_intersect_comparison(duckdb_connection, a, b)

    # Assert
    assert comparison.match, comparison.failure_message()
    assert comparison.giql_row_count == 2


def test_intersects_should_match_bedtools_when_intervals_are_very_large(
    duckdb_connection,
):
    """Test INTERSECTS matches bedtools for multi-megabase intervals.

    Given:
        Very large genomic intervals spanning millions of bases
    When:
        GIQL INTERSECTS is compared to bedtools intersect
    Then:
        It should produce the same overlap result as bedtools
    """
    # Arrange
    a = [GenomicInterval("chr1", 0, 10_000_000, "a1", 0, "+")]
    b = [GenomicInterval("chr1", 5_000_000, 15_000_000, "b1", 0, "+")]

    # Act
    comparison = _run_intersect_comparison(duckdb_connection, a, b)

    # Assert
    assert comparison.match, comparison.failure_message()
    assert comparison.giql_row_count == 1


def test_intersects_should_match_bedtools_at_scale(duckdb_connection):
    """Test INTERSECTS matches bedtools on a larger generated dataset.

    Given:
        A generated dataset with 100 intervals per chromosome on 3 chromosomes
    When:
        GIQL INTERSECTS is compared to bedtools intersect
    Then:
        It should match bedtools on the full dataset
    """
    # Arrange
    intervals_a = generate_random_intervals(
        seed=42,
        prefix="a",
        count_per_chrom=100,
        n_chroms=3,
        start_max=900_000,
    )
    intervals_b = generate_random_intervals(
        seed=43,
        prefix="b",
        count_per_chrom=100,
        n_chroms=3,
        start_max=900_000,
    )

    # Act
    comparison = _run_intersect_comparison(duckdb_connection, intervals_a, intervals_b)

    # Assert
    assert comparison.match, comparison.failure_message()


def test_intersects_should_match_bedtools_when_same_strand_filter_applied(
    duckdb_connection,
):
    """Test INTERSECTS with same-strand filter matches bedtools -s.

    Given:
        Overlapping intervals with mixed strand orientations
    When:
        GIQL INTERSECTS with a same-strand filter is compared to bedtools -s
    Then:
        It should return only the same-strand overlaps
    """
    # Arrange
    a = [
        GenomicInterval("chr1", 100, 200, "a_plus", 0, "+"),
        GenomicInterval("chr1", 100, 200, "a_minus", 0, "-"),
    ]
    b = [GenomicInterval("chr1", 150, 250, "b_plus", 0, "+")]

    # Act
    comparison = _run_intersect_comparison(
        duckdb_connection,
        a,
        b,
        strand_filter="a.strand = b.strand",
    )

    # Assert
    assert comparison.match, comparison.failure_message()
    assert comparison.giql_row_count == 1


def test_intersects_should_match_bedtools_when_opposite_strand_filter_applied(
    duckdb_connection,
):
    """Test INTERSECTS with opposite-strand filter matches bedtools -S.

    Given:
        Overlapping intervals with mixed strand orientations
    When:
        GIQL INTERSECTS with an opposite-strand filter is compared to bedtools -S
    Then:
        It should return only the opposite-strand overlaps
    """
    # Arrange
    a = [
        GenomicInterval("chr1", 100, 200, "a_plus", 0, "+"),
        GenomicInterval("chr1", 100, 200, "a_minus", 0, "-"),
    ]
    b = [GenomicInterval("chr1", 150, 250, "b_plus", 0, "+")]

    # Act
    comparison = _run_intersect_comparison(
        duckdb_connection,
        a,
        b,
        strand_filter="a.strand != b.strand",
    )

    # Assert
    assert comparison.match, comparison.failure_message()
    assert comparison.giql_row_count == 1
