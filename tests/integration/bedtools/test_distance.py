"""Integration tests for GIQL DISTANCE operator.

These tests validate that GIQL's DISTANCE function computes the correct
genomic distance between intervals, cross-validated against bedtools
closest -d output.
"""

import pytest

from .utils.bedtools_wrapper import closest
from .utils.data_models import GenomicInterval


def _oracle_distance(a: GenomicInterval, b: GenomicInterval) -> int:
    """Return the bedtools ``closest -d`` distance between a single A and B."""
    return closest([a.to_tuple()], [b.to_tuple()])[0][12]


def test_distance_non_overlapping(giql_query):
    """
    GIVEN two non-overlapping intervals with a known gap
    WHEN DISTANCE is computed via GIQL
    THEN the distance matches live bedtools closest -d
    """
    a = GenomicInterval("chr1", 100, 200, "a1", 100, "+")
    b = GenomicInterval("chr1", 300, 400, "b1", 100, "+")
    result = giql_query(
        """
        SELECT DISTANCE(a.interval, b.interval) AS dist
        FROM intervals_a a, intervals_b b
        """,
        tables=["intervals_a", "intervals_b"],
        intervals_a=[a],
        intervals_b=[b],
    )

    assert len(result) == 1
    assert result[0][0] == _oracle_distance(a, b)


def test_distance_overlapping(giql_query):
    """
    GIVEN two overlapping intervals
    WHEN DISTANCE is computed via GIQL
    THEN the distance is 0, matching live bedtools closest -d
    """
    a = GenomicInterval("chr1", 100, 300, "a1", 100, "+")
    b = GenomicInterval("chr1", 200, 400, "b1", 100, "+")
    result = giql_query(
        """
        SELECT DISTANCE(a.interval, b.interval) AS dist
        FROM intervals_a a, intervals_b b
        """,
        tables=["intervals_a", "intervals_b"],
        intervals_a=[a],
        intervals_b=[b],
    )

    assert len(result) == 1
    assert result[0][0] == _oracle_distance(a, b)


def test_distance_adjacent(giql_query):
    """
    GIVEN two adjacent intervals (touching, half-open coordinates)
    WHEN DISTANCE is computed via GIQL
    THEN the distance is 1, matching live bedtools closest -d for book-ended features
    """
    a = GenomicInterval("chr1", 100, 200, "a1", 100, "+")
    b = GenomicInterval("chr1", 200, 300, "b1", 100, "+")
    result = giql_query(
        """
        SELECT DISTANCE(a.interval, b.interval) AS dist
        FROM intervals_a a, intervals_b b
        """,
        tables=["intervals_a", "intervals_b"],
        intervals_a=[a],
        intervals_b=[b],
    )

    assert len(result) == 1
    assert result[0][0] == _oracle_distance(a, b)


@pytest.mark.parametrize(
    ("b_start", "b_end"),
    [
        (150, 250),  # overlapping
        (200, 300),  # book-ended downstream
        (201, 300),  # 1 bp gap
        (202, 300),  # 2 bp gap
        (300, 400),  # 100 bp gap
        (10200, 10300),  # large gap
    ],
    ids=["overlap", "book_ended", "gap_1bp", "gap_2bp", "gap_100bp", "gap_large"],
)
def test_distance_matches_bedtools_across_gap_sizes(giql_query, b_start, b_end):
    """
    GIVEN A chr1:100-200 and B at a range of gaps from overlap to a large gap
    WHEN DISTANCE is computed via GIQL
    THEN the result equals live bedtools closest -d for every gap, so a future
        off-by-one in either direction fails automatically
    """
    a = GenomicInterval("chr1", 100, 200, "a1", 100, "+")
    b = GenomicInterval("chr1", b_start, b_end, "b1", 100, "+")
    result = giql_query(
        """
        SELECT DISTANCE(a.interval, b.interval) AS dist
        FROM intervals_a a, intervals_b b
        """,
        tables=["intervals_a", "intervals_b"],
        intervals_a=[a],
        intervals_b=[b],
    )

    assert len(result) == 1
    assert result[0][0] == _oracle_distance(a, b)


def test_distance_book_ended_upstream_matches_bedtools(giql_query):
    """
    GIVEN a book-ended pair with B upstream of A (A chr1:200-300, B chr1:100-200)
    WHEN DISTANCE is computed via GIQL
    THEN the result is 1, matching bedtools, confirming direction symmetry of
        the parity +1
    """
    a = GenomicInterval("chr1", 200, 300, "a1", 100, "+")
    b = GenomicInterval("chr1", 100, 200, "b1", 100, "+")
    result = giql_query(
        """
        SELECT DISTANCE(a.interval, b.interval) AS dist
        FROM intervals_a a, intervals_b b
        """,
        tables=["intervals_a", "intervals_b"],
        intervals_a=[a],
        intervals_b=[b],
    )

    assert len(result) == 1
    assert result[0][0] == _oracle_distance(a, b)


def test_distance_multi_row_cross_join_matches_bedtools(giql_query):
    """
    GIVEN two A intervals and two B intervals at mixed gaps in one query
    WHEN DISTANCE is computed for every (a, b) pair via a GIQL cross join
    THEN each pair's distance equals live bedtools closest -d for that pair,
        exercising the multi-row table-scan path
    """
    a_rows = [
        GenomicInterval("chr1", 100, 200, "a1", 100, "+"),
        GenomicInterval("chr1", 500, 600, "a2", 100, "+"),
    ]
    b_rows = [
        GenomicInterval("chr1", 200, 300, "b1", 100, "+"),
        GenomicInterval("chr1", 800, 900, "b2", 100, "+"),
    ]
    result = giql_query(
        """
        SELECT a.name, b.name, DISTANCE(a.interval, b.interval) AS dist
        FROM intervals_a a, intervals_b b
        """,
        tables=["intervals_a", "intervals_b"],
        intervals_a=a_rows,
        intervals_b=b_rows,
    )

    by_pair = {(row[0], row[1]): row[2] for row in result}
    assert len(by_pair) == 4
    for a in a_rows:
        for b in b_rows:
            assert by_pair[(a.name, b.name)] == _oracle_distance(a, b)


def test_distance_cross_chromosome(giql_query):
    """
    GIVEN two intervals on different chromosomes
    WHEN DISTANCE is computed via GIQL
    THEN the distance is NULL (cross-chromosome distance is undefined)
    """
    result = giql_query(
        """
        SELECT DISTANCE(a.interval, b.interval) AS dist
        FROM intervals_a a, intervals_b b
        """,
        tables=["intervals_a", "intervals_b"],
        intervals_a=[GenomicInterval("chr1", 100, 200, "a1", 100, "+")],
        intervals_b=[GenomicInterval("chr2", 100, 200, "b1", 100, "+")],
    )

    assert len(result) == 1
    assert result[0][0] is None, (
        f"Cross-chromosome distance should be NULL, got {result[0][0]}"
    )


def test_distance_signed_downstream(giql_query):
    """
    GIVEN B is downstream of A (B starts after A ends) on + strand
    WHEN DISTANCE with signed := true is computed
    THEN the distance is +101, the gap (100) plus the bedtools parity 1
    """
    result = giql_query(
        """
        SELECT DISTANCE(a.interval, b.interval, signed := true) AS dist
        FROM intervals_a a, intervals_b b
        """,
        tables=["intervals_a", "intervals_b"],
        intervals_a=[GenomicInterval("chr1", 100, 200, "a1", 100, "+")],
        intervals_b=[GenomicInterval("chr1", 300, 400, "b1", 100, "+")],
    )

    assert len(result) == 1
    assert result[0][0] == 101, f"Expected +101 (downstream), got {result[0][0]}"


def test_distance_signed_upstream(giql_query):
    """
    GIVEN B is upstream of A (B ends before A starts) on + strand
    WHEN DISTANCE with signed := true is computed
    THEN the distance is -101, the gap (100) plus parity 1, negated for upstream
    """
    result = giql_query(
        """
        SELECT DISTANCE(a.interval, b.interval, signed := true) AS dist
        FROM intervals_a a, intervals_b b
        """,
        tables=["intervals_a", "intervals_b"],
        intervals_a=[GenomicInterval("chr1", 300, 400, "a1", 100, "+")],
        intervals_b=[GenomicInterval("chr1", 100, 200, "b1", 100, "+")],
    )

    assert len(result) == 1
    assert result[0][0] == -101, f"Expected -101 (upstream), got {result[0][0]}"


def test_distance_stranded_unstranded_input(giql_query):
    """
    GIVEN one interval with strand "." (unstranded)
    WHEN DISTANCE with stranded := true is computed
    THEN the distance is NULL (stranded mode requires valid strand)
    """
    result = giql_query(
        """
        SELECT DISTANCE(a.interval, b.interval, stranded := true) AS dist
        FROM intervals_a a, intervals_b b
        """,
        tables=["intervals_a", "intervals_b"],
        intervals_a=[GenomicInterval("chr1", 100, 200, "a1", 100, ".")],
        intervals_b=[GenomicInterval("chr1", 300, 400, "b1", 100, "+")],
    )

    assert len(result) == 1
    assert result[0][0] is None, (
        f"Stranded distance with '.' strand should be NULL, got {result[0][0]}"
    )


def test_distance_stranded_same_strand(giql_query):
    """
    GIVEN two non-overlapping intervals both on + strand
    WHEN DISTANCE with stranded := true is computed
    THEN the distance matches live bedtools closest -d (same + strand, no flip)
    """
    a = GenomicInterval("chr1", 100, 200, "a1", 100, "+")
    b = GenomicInterval("chr1", 300, 400, "b1", 100, "+")
    result = giql_query(
        """
        SELECT DISTANCE(a.interval, b.interval, stranded := true) AS dist
        FROM intervals_a a, intervals_b b
        """,
        tables=["intervals_a", "intervals_b"],
        intervals_a=[a],
        intervals_b=[b],
    )

    assert len(result) == 1
    # Same + strand applies no flip, so the value matches the unsigned oracle.
    assert result[0][0] == _oracle_distance(a, b)


def test_distance_signed_stranded_minus_strand(giql_query):
    """
    GIVEN two non-overlapping intervals on - strand, B downstream genomically
    WHEN DISTANCE with signed := true, stranded := true is computed
    THEN the sign is inverted due to - strand (downstream
    genomically is upstream in transcript orientation)
    """
    result = giql_query(
        """
        SELECT DISTANCE(
            a.interval, b.interval,
            signed := true,
            stranded := true
        ) AS dist
        FROM intervals_a a, intervals_b b
        """,
        tables=["intervals_a", "intervals_b"],
        intervals_a=[GenomicInterval("chr1", 100, 200, "a1", 100, "-")],
        intervals_b=[GenomicInterval("chr1", 300, 400, "b1", 100, "-")],
    )

    assert len(result) == 1
    # On - strand with signed+stranded, genomic downstream becomes transcript
    # upstream, so the sign flips: -(gap 100 + parity 1) = -101.
    assert result[0][0] == -101, f"Expected -101 on - strand, got {result[0][0]}"
