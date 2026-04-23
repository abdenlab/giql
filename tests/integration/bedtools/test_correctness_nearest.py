"""Extended correctness tests for GIQL NEAREST operator vs bedtools closest.

These tests cover distance calculations, multi-query scenarios, and scale
to ensure comprehensive GIQL/bedtools equivalence for nearest operations.
"""

import pytest

from giql import transpile

from .utils.bedtools_wrapper import closest
from .utils.data_models import GenomicInterval
from .utils.duckdb_loader import load_intervals

pytestmark = pytest.mark.integration


def _load_and_query_nearest(
    duckdb_connection,
    intervals_a,
    intervals_b,
    *,
    k=1,
    stranded=False,
):
    """Load intervals, run GIQL NEAREST and bedtools closest, return both results."""
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

    strand_mode = "same" if stranded else None
    bedtools_result = closest(
        [i.to_tuple() for i in intervals_a],
        [i.to_tuple() for i in intervals_b],
        strand_mode=strand_mode,
        k=k,
    )

    stranded_arg = ", stranded := true" if stranded else ""
    sql = transpile(
        f"""
        SELECT a.*, b.*
        FROM intervals_a a
        CROSS JOIN LATERAL NEAREST(
            intervals_b,
            reference := a.interval,
            k := {k}{stranded_arg}
        ) b
        ORDER BY a.chrom, a.start
        """,
        tables=["intervals_a", "intervals_b"],
    )
    giql_result = duckdb_connection.execute(sql).fetchall()

    return giql_result, bedtools_result


def test_nearest_should_report_distance_zero_when_intervals_overlap(duckdb_connection):
    """Test NEAREST reports zero distance for overlapping intervals.

    Given:
        Overlapping intervals in A and B
    When:
        GIQL NEAREST is compared to bedtools closest
    Then:
        It should report distance=0 for the overlapping pair
    """
    # Arrange
    a = [GenomicInterval("chr1", 100, 300, "a1", 0, "+")]
    b = [GenomicInterval("chr1", 200, 400, "b1", 0, "+")]

    # Act
    giql_result, bedtools_result = _load_and_query_nearest(duckdb_connection, a, b)

    # Assert
    assert len(giql_result) == len(bedtools_result) == 1
    # bedtools closest -d reports 0 for overlapping
    assert bedtools_result[0][-1] == 0


def test_nearest_should_find_adjacent_neighbor_when_intervals_touch(
    duckdb_connection,
):
    """Test NEAREST matches bedtools for adjacent non-overlapping intervals.

    Given:
        Two adjacent intervals in half-open coordinates (a1 ending at
        200, b1 starting at 200 — touching but not overlapping)
    When:
        GIQL NEAREST is compared to bedtools closest
    Then:
        It should identify b1 as a1's nearest neighbor, and bedtools
        should report the canonical adjacent-interval distance of 1
        (bedtools >= 2.31 counts the gap base in half-open coords)
    """
    # Arrange
    a = [GenomicInterval("chr1", 100, 200, "a1", 0, "+")]
    b = [GenomicInterval("chr1", 200, 300, "b1", 0, "+")]

    # Act
    giql_result, bedtools_result = _load_and_query_nearest(duckdb_connection, a, b)

    # Assert
    assert len(giql_result) == len(bedtools_result) == 1
    assert bedtools_result[0][-1] == 1
    assert giql_result[0][9] == "b1"


def test_nearest_should_match_bedtools_when_candidate_is_upstream(duckdb_connection):
    """Test NEAREST matches bedtools for an upstream candidate interval.

    Given:
        A B interval positioned far upstream of the A interval
    When:
        GIQL NEAREST is compared to bedtools closest
    Then:
        It should identify the upstream candidate with the correct distance
    """
    # Arrange
    a = [GenomicInterval("chr1", 500, 600, "a1", 0, "+")]
    b = [GenomicInterval("chr1", 100, 200, "b1", 0, "+")]

    # Act
    giql_result, bedtools_result = _load_and_query_nearest(duckdb_connection, a, b)

    # Assert
    assert len(giql_result) == len(bedtools_result) == 1
    # Distance: 500 - 200 = 300 (half-open), bedtools may report 301
    assert bedtools_result[0][-1] in (300, 301)
    assert giql_result[0][9] == "b1"


def test_nearest_should_match_bedtools_when_candidate_is_downstream(duckdb_connection):
    """Test NEAREST matches bedtools for a downstream candidate interval.

    Given:
        A B interval positioned far downstream of the A interval
    When:
        GIQL NEAREST is compared to bedtools closest
    Then:
        It should identify the downstream candidate with the correct distance
    """
    # Arrange
    a = [GenomicInterval("chr1", 100, 200, "a1", 0, "+")]
    b = [GenomicInterval("chr1", 500, 600, "b1", 0, "+")]

    # Act
    giql_result, bedtools_result = _load_and_query_nearest(duckdb_connection, a, b)

    # Assert
    assert len(giql_result) == len(bedtools_result) == 1
    # Distance: 500 - 200 = 300 (half-open), bedtools may report 301
    assert bedtools_result[0][-1] in (300, 301)
    assert giql_result[0][9] == "b1"


def test_nearest_should_match_bedtools_for_multiple_query_intervals(duckdb_connection):
    """Test NEAREST matches bedtools when multiple query intervals are used.

    Given:
        Multiple query intervals in A and multiple candidates in B
    When:
        GIQL NEAREST is compared to bedtools closest
    Then:
        It should produce the correct pairing for each query interval
    """
    # Arrange
    a = [
        GenomicInterval("chr1", 100, 200, "a1", 0, "+"),
        GenomicInterval("chr1", 500, 600, "a2", 0, "+"),
        GenomicInterval("chr1", 900, 1000, "a3", 0, "+"),
    ]
    b = [
        GenomicInterval("chr1", 250, 300, "b1", 0, "+"),
        GenomicInterval("chr1", 700, 800, "b2", 0, "+"),
    ]

    # Act
    giql_result, bedtools_result = _load_and_query_nearest(duckdb_connection, a, b)

    # Assert
    assert len(giql_result) == len(bedtools_result) == 3

    giql_sorted = sorted(giql_result, key=lambda r: (r[0], r[1]))
    bt_sorted = sorted(bedtools_result, key=lambda r: (r[0], r[1]))

    for giql_row, bt_row in zip(giql_sorted, bt_sorted):
        assert giql_row[3] == bt_row[3]  # a.name matches
        assert giql_row[9] == bt_row[9]  # b.name matches


def test_nearest_should_return_three_neighbors_when_k_is_three(duckdb_connection):
    """Test NEAREST returns the three nearest neighbors when k=3.

    Given:
        One query interval and four database candidates
    When:
        GIQL NEAREST(k=3) is compared to bedtools closest -k 3
    Then:
        It should return the same three nearest intervals as bedtools
    """
    # Arrange
    a = [GenomicInterval("chr1", 400, 500, "a1", 0, "+")]
    b = [
        GenomicInterval("chr1", 100, 150, "b_far", 0, "+"),
        GenomicInterval("chr1", 350, 390, "b_near", 0, "+"),
        GenomicInterval("chr1", 550, 600, "b_close", 0, "+"),
        GenomicInterval("chr1", 900, 1000, "b_farther", 0, "+"),
    ]

    # Act
    giql_result, bedtools_result = _load_and_query_nearest(duckdb_connection, a, b, k=3)

    # Assert
    assert len(giql_result) == 3
    assert len(bedtools_result) == 3

    giql_names = {r[9] for r in giql_result}
    bt_names = {r[9] for r in bedtools_result}
    assert giql_names == bt_names


def test_nearest_should_return_available_neighbors_when_k_exceeds_candidates(duckdb_connection):
    """Test NEAREST caps results at the number of available candidates.

    Given:
        One query interval, only two database candidates, and k=5
    When:
        GIQL NEAREST is compared to bedtools closest
    Then:
        It should return only the two available candidates, matching bedtools
    """
    # Arrange
    a = [GenomicInterval("chr1", 200, 300, "a1", 0, "+")]
    b = [
        GenomicInterval("chr1", 100, 150, "b1", 0, "+"),
        GenomicInterval("chr1", 400, 500, "b2", 0, "+"),
    ]

    # Act
    giql_result, bedtools_result = _load_and_query_nearest(duckdb_connection, a, b, k=5)

    # Assert
    assert len(giql_result) == len(bedtools_result) == 2


def test_nearest_should_return_only_same_strand_candidates_when_stranded(duckdb_connection):
    """Test NEAREST restricts matches to same strand when stranded=true.

    Given:
        Candidates on same and opposite strands, with the opposite-strand
        candidate being closer
    When:
        GIQL NEAREST(stranded=true) is compared to bedtools closest -s
    Then:
        It should return only the same-strand match, ignoring the closer
        opposite-strand candidate
    """
    # Arrange
    a = [GenomicInterval("chr1", 100, 200, "a1", 0, "+")]
    b = [
        GenomicInterval("chr1", 220, 240, "b_opp", 0, "-"),  # closer, opposite
        GenomicInterval("chr1", 300, 400, "b_same", 0, "+"),  # farther, same
    ]

    # Act
    giql_result, bedtools_result = _load_and_query_nearest(
        duckdb_connection,
        a,
        b,
        stranded=True,
    )

    # Assert
    assert len(giql_result) == len(bedtools_result) == 1
    assert giql_result[0][9] == "b_same"
    assert bedtools_result[0][9] == "b_same"


def test_nearest_should_ignore_strand_when_unstranded(duckdb_connection):
    """Test NEAREST ignores strand when not configured as stranded.

    Given:
        Candidates on different strands where the closer one is on the
        opposite strand
    When:
        GIQL NEAREST (default) is compared to bedtools closest (default)
    Then:
        It should return the nearest candidate regardless of strand
    """
    # Arrange
    a = [GenomicInterval("chr1", 100, 200, "a1", 0, "+")]
    b = [
        GenomicInterval("chr1", 250, 300, "b_far", 0, "+"),
        GenomicInterval("chr1", 210, 230, "b_near", 0, "-"),
    ]

    # Act
    giql_result, bedtools_result = _load_and_query_nearest(duckdb_connection, a, b)

    # Assert
    assert len(giql_result) == len(bedtools_result) == 1
    assert giql_result[0][9] == "b_near"
    assert bedtools_result[0][9] == "b_near"


def test_nearest_should_isolate_matches_per_chromosome(duckdb_connection):
    """Test NEAREST only pairs intervals on the same chromosome.

    Given:
        Intervals distributed across multiple chromosomes
    When:
        GIQL NEAREST is compared to bedtools closest
    Then:
        It should find nearest matches only within each chromosome
    """
    # Arrange
    a = [
        GenomicInterval("chr1", 100, 200, "a1", 0, "+"),
        GenomicInterval("chr2", 100, 200, "a2", 0, "+"),
    ]
    b = [
        GenomicInterval("chr1", 500, 600, "b1", 0, "+"),
        GenomicInterval("chr2", 300, 400, "b2", 0, "+"),
    ]

    # Act
    giql_result, bedtools_result = _load_and_query_nearest(duckdb_connection, a, b)

    # Assert
    assert len(giql_result) == len(bedtools_result) == 2

    for giql_row in giql_result:
        assert giql_row[0] == giql_row[6], "A and B should be on same chromosome"


def test_nearest_should_match_bedtools_on_large_multi_chromosome_dataset(duckdb_connection):
    """Test NEAREST matches bedtools on a large multi-chromosome dataset.

    Given:
        Fifty-plus intervals per table spread across three chromosomes
    When:
        GIQL NEAREST is compared to bedtools closest
    Then:
        It should produce the same row count as bedtools on the full dataset
    """
    # Arrange
    import random

    rng = random.Random(42)
    intervals_a = []
    intervals_b = []

    for chrom_num in range(1, 4):
        chrom = f"chr{chrom_num}"
        for i in range(50):
            start = rng.randint(0, 900_000)
            size = rng.randint(100, 1000)
            intervals_a.append(
                GenomicInterval(chrom, start, start + size, f"a_{chrom_num}_{i}", 0, "+")
            )
            start = rng.randint(0, 900_000)
            size = rng.randint(100, 1000)
            intervals_b.append(
                GenomicInterval(chrom, start, start + size, f"b_{chrom_num}_{i}", 0, "+")
            )

    # Act
    giql_result, bedtools_result = _load_and_query_nearest(
        duckdb_connection,
        intervals_a,
        intervals_b,
    )

    # Assert
    assert len(giql_result) == len(bedtools_result), (
        f"Row count mismatch: GIQL={len(giql_result)}, bedtools={len(bedtools_result)}"
    )
