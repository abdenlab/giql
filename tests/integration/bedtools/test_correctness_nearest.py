"""Extended correctness tests for GIQL NEAREST operator vs bedtools closest.

These tests cover distance calculations, multi-query scenarios, and scale
to ensure comprehensive GIQL/bedtools equivalence for nearest operations.
"""

from giql import transpile

from .utils.bedtools_wrapper import closest
from .utils.data_models import GenomicInterval
from .utils.duckdb_loader import load_intervals


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


def test_nearest_overlapping_distance_zero(duckdb_connection):
    """
    GIVEN overlapping intervals in A and B
    WHEN GIQL NEAREST is compared to bedtools closest
    THEN overlapping intervals report distance=0 in bedtools
    """
    a = [GenomicInterval("chr1", 100, 300, "a1", 0, "+")]
    b = [GenomicInterval("chr1", 200, 400, "b1", 0, "+")]
    giql_result, bedtools_result = _load_and_query_nearest(duckdb_connection, a, b)

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


def test_nearest_upstream_distance(duckdb_connection):
    """
    GIVEN B interval far upstream of A
    WHEN GIQL NEAREST is compared to bedtools closest
    THEN distance calculated correctly
    """
    a = [GenomicInterval("chr1", 500, 600, "a1", 0, "+")]
    b = [GenomicInterval("chr1", 100, 200, "b1", 0, "+")]
    giql_result, bedtools_result = _load_and_query_nearest(duckdb_connection, a, b)

    assert len(giql_result) == len(bedtools_result) == 1
    # Distance: 500 - 200 = 300 (half-open), bedtools may report 301
    assert bedtools_result[0][-1] in (300, 301)
    assert giql_result[0][9] == "b1"


def test_nearest_downstream_distance(duckdb_connection):
    """
    GIVEN B interval far downstream of A
    WHEN GIQL NEAREST is compared to bedtools closest
    THEN distance calculated correctly
    """
    a = [GenomicInterval("chr1", 100, 200, "a1", 0, "+")]
    b = [GenomicInterval("chr1", 500, 600, "b1", 0, "+")]
    giql_result, bedtools_result = _load_and_query_nearest(duckdb_connection, a, b)

    assert len(giql_result) == len(bedtools_result) == 1
    # Distance: 500 - 200 = 300 (half-open), bedtools may report 301
    assert bedtools_result[0][-1] in (300, 301)
    assert giql_result[0][9] == "b1"


def test_nearest_multi_query_correctness(duckdb_connection):
    """
    GIVEN multiple query intervals and multiple candidates
    WHEN GIQL NEAREST is compared to bedtools closest
    THEN correct pairing for each query interval
    """
    a = [
        GenomicInterval("chr1", 100, 200, "a1", 0, "+"),
        GenomicInterval("chr1", 500, 600, "a2", 0, "+"),
        GenomicInterval("chr1", 900, 1000, "a3", 0, "+"),
    ]
    b = [
        GenomicInterval("chr1", 250, 300, "b1", 0, "+"),
        GenomicInterval("chr1", 700, 800, "b2", 0, "+"),
    ]
    giql_result, bedtools_result = _load_and_query_nearest(duckdb_connection, a, b)

    assert len(giql_result) == len(bedtools_result) == 3

    giql_sorted = sorted(giql_result, key=lambda r: (r[0], r[1]))
    bt_sorted = sorted(bedtools_result, key=lambda r: (r[0], r[1]))

    for giql_row, bt_row in zip(giql_sorted, bt_sorted):
        assert giql_row[3] == bt_row[3]  # a.name matches
        assert giql_row[9] == bt_row[9]  # b.name matches


def test_nearest_k3_correctness(duckdb_connection):
    """
    GIVEN one query interval and 4 database intervals
    WHEN GIQL NEAREST(k=3) is compared to bedtools closest -k 3
    THEN both return 3 nearest intervals
    """
    a = [GenomicInterval("chr1", 400, 500, "a1", 0, "+")]
    b = [
        GenomicInterval("chr1", 100, 150, "b_far", 0, "+"),
        GenomicInterval("chr1", 350, 390, "b_near", 0, "+"),
        GenomicInterval("chr1", 550, 600, "b_close", 0, "+"),
        GenomicInterval("chr1", 900, 1000, "b_farther", 0, "+"),
    ]
    giql_result, bedtools_result = _load_and_query_nearest(duckdb_connection, a, b, k=3)

    assert len(giql_result) == 3
    assert len(bedtools_result) == 3

    giql_names = {r[9] for r in giql_result}
    bt_names = {r[9] for r in bedtools_result}
    assert giql_names == bt_names


def test_nearest_k_exceeds_available_correctness(duckdb_connection):
    """
    GIVEN one query and only 2 database intervals, k=5
    WHEN GIQL NEAREST is compared to bedtools closest
    THEN both return only 2 (available) results
    """
    a = [GenomicInterval("chr1", 200, 300, "a1", 0, "+")]
    b = [
        GenomicInterval("chr1", 100, 150, "b1", 0, "+"),
        GenomicInterval("chr1", 400, 500, "b2", 0, "+"),
    ]
    giql_result, bedtools_result = _load_and_query_nearest(duckdb_connection, a, b, k=5)

    assert len(giql_result) == len(bedtools_result) == 2


def test_nearest_same_strand_correctness(duckdb_connection):
    """
    GIVEN intervals with candidates on same and opposite strands
    WHEN GIQL NEAREST(stranded=true) is compared to bedtools closest -s
    THEN only same-strand matches
    """
    a = [GenomicInterval("chr1", 100, 200, "a1", 0, "+")]
    b = [
        GenomicInterval("chr1", 220, 240, "b_opp", 0, "-"),  # closer, opposite
        GenomicInterval("chr1", 300, 400, "b_same", 0, "+"),  # farther, same
    ]
    giql_result, bedtools_result = _load_and_query_nearest(
        duckdb_connection,
        a,
        b,
        stranded=True,
    )

    assert len(giql_result) == len(bedtools_result) == 1
    assert giql_result[0][9] == "b_same"
    assert bedtools_result[0][9] == "b_same"


def test_nearest_strand_ignorant_correctness(duckdb_connection):
    """
    GIVEN intervals on different strands
    WHEN GIQL NEAREST (default) is compared to bedtools closest (default)
    THEN nearest found regardless of strand
    """
    a = [GenomicInterval("chr1", 100, 200, "a1", 0, "+")]
    b = [
        GenomicInterval("chr1", 250, 300, "b_far", 0, "+"),
        GenomicInterval("chr1", 210, 230, "b_near", 0, "-"),
    ]
    giql_result, bedtools_result = _load_and_query_nearest(duckdb_connection, a, b)

    assert len(giql_result) == len(bedtools_result) == 1
    assert giql_result[0][9] == "b_near"
    assert bedtools_result[0][9] == "b_near"


def test_nearest_cross_chromosome_isolation(duckdb_connection):
    """
    GIVEN intervals on multiple chromosomes
    WHEN GIQL NEAREST is compared to bedtools closest
    THEN nearest found per-chromosome only
    """
    a = [
        GenomicInterval("chr1", 100, 200, "a1", 0, "+"),
        GenomicInterval("chr2", 100, 200, "a2", 0, "+"),
    ]
    b = [
        GenomicInterval("chr1", 500, 600, "b1", 0, "+"),
        GenomicInterval("chr2", 300, 400, "b2", 0, "+"),
    ]
    giql_result, bedtools_result = _load_and_query_nearest(duckdb_connection, a, b)

    assert len(giql_result) == len(bedtools_result) == 2

    for giql_row in giql_result:
        assert giql_row[0] == giql_row[6], "A and B should be on same chromosome"


def test_nearest_large_scale(duckdb_connection):
    """
    GIVEN 50+ intervals per table across 3 chromosomes
    WHEN GIQL NEAREST is compared to bedtools closest
    THEN row counts match on the full dataset
    """
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

    giql_result, bedtools_result = _load_and_query_nearest(
        duckdb_connection,
        intervals_a,
        intervals_b,
    )

    assert len(giql_result) == len(bedtools_result), (
        f"Row count mismatch: GIQL={len(giql_result)}, bedtools={len(bedtools_result)}"
    )
