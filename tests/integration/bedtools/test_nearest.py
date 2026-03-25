"""Integration tests for GIQL NEAREST operator.

These tests validate that GIQL's NEAREST operator produces results
consistent with bedtools closest command.
"""

from giql import transpile

from .utils.bed_export import load_intervals
from .utils.bedtools_wrapper import closest
from .utils.data_models import GenomicInterval


def test_nearest_non_overlapping(duckdb_connection):
    """
    Given:
        Two sets of non-overlapping intervals
    When:
        NEAREST operator is applied
    Then:
        Each interval in A finds its closest neighbor in B
    """
    intervals_a = [
        GenomicInterval("chr1", 100, 200, "a1", 100, "+"),
        GenomicInterval("chr1", 500, 600, "a2", 150, "+"),
    ]
    intervals_b = [
        GenomicInterval("chr1", 250, 300, "b1", 100, "+"),
        GenomicInterval("chr1", 350, 400, "b2", 150, "+"),
        GenomicInterval("chr1", 700, 800, "b3", 200, "+"),
    ]

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

    bedtools_result = closest(
        [i.to_tuple() for i in intervals_a],
        [i.to_tuple() for i in intervals_b],
    )

    sql = transpile(
        """
        SELECT a.*, b.*
        FROM intervals_a a
        CROSS JOIN LATERAL NEAREST(
            intervals_b,
            reference := a.interval,
            k := 1
        ) b
        ORDER BY a.chrom, a.start
        """,
        tables=["intervals_a", "intervals_b"],
    )
    giql_result = duckdb_connection.execute(sql).fetchall()

    # Verify row counts match
    assert len(giql_result) == len(bedtools_result), (
        f"Row count mismatch: GIQL={len(giql_result)}, bedtools={len(bedtools_result)}"
    )

    # Verify each A interval found the correct B neighbor
    for giql_row, bt_row in zip(
        sorted(giql_result, key=lambda r: (r[0], r[1])),
        sorted(bedtools_result, key=lambda r: (r[0], r[1])),
    ):
        # Compare A interval name
        assert giql_row[3] == bt_row[3], (
            f"A name mismatch: GIQL={giql_row[3]}, bedtools={bt_row[3]}"
        )
        # Compare matched B interval name
        assert giql_row[9] == bt_row[9], (
            f"B name mismatch: GIQL={giql_row[9]}, bedtools={bt_row[9]}"
        )


def test_nearest_multiple_candidates(duckdb_connection):
    """
    Given:
        Interval in A with multiple equidistant intervals in B
    When:
        NEAREST operator is applied
    Then:
        One of the equidistant intervals is returned
    """
    intervals_a = [
        GenomicInterval("chr1", 300, 400, "a1", 100, "+"),
    ]
    intervals_b = [
        GenomicInterval("chr1", 100, 200, "b1", 100, "+"),
        GenomicInterval("chr1", 500, 600, "b2", 150, "+"),
    ]

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

    sql = transpile(
        """
        SELECT a.*, b.*
        FROM intervals_a a
        CROSS JOIN LATERAL NEAREST(
            intervals_b,
            reference := a.interval,
            k := 1
        ) b
        """,
        tables=["intervals_a", "intervals_b"],
    )
    giql_result = duckdb_connection.execute(sql).fetchall()

    assert len(giql_result) == 1
    assert giql_result[0][3] == "a1"
    # Nearest could be either b1 or b2 (both equidistant at 100bp)
    assert giql_result[0][9] in ("b1", "b2")


def test_nearest_cross_chromosome(duckdb_connection):
    """
    Given:
        Intervals on different chromosomes
    When:
        NEAREST operator is applied
    Then:
        Each interval finds nearest only on the same chromosome
    """
    intervals_a = [
        GenomicInterval("chr1", 100, 200, "a1", 100, "+"),
        GenomicInterval("chr2", 100, 200, "a2", 150, "+"),
    ]
    intervals_b = [
        GenomicInterval("chr1", 300, 400, "b1", 100, "+"),
        GenomicInterval("chr2", 300, 400, "b2", 150, "+"),
    ]

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

    bedtools_result = closest(
        [i.to_tuple() for i in intervals_a],
        [i.to_tuple() for i in intervals_b],
    )

    sql = transpile(
        """
        SELECT a.*, b.*
        FROM intervals_a a
        CROSS JOIN LATERAL NEAREST(
            intervals_b,
            reference := a.interval,
            k := 1
        ) b
        ORDER BY a.chrom, a.start
        """,
        tables=["intervals_a", "intervals_b"],
    )
    giql_result = duckdb_connection.execute(sql).fetchall()

    assert len(giql_result) == len(bedtools_result)

    for giql_row, bt_row in zip(
        sorted(giql_result, key=lambda r: (r[0], r[1])),
        sorted(bedtools_result, key=lambda r: (r[0], r[1])),
    ):
        # A and matched B should be on the same chromosome
        assert giql_row[0] == giql_row[6], (
            f"Cross-chromosome match: A on {giql_row[0]}, B on {giql_row[6]}"
        )
        assert giql_row[3] == bt_row[3]
        assert giql_row[9] == bt_row[9]


def test_nearest_boundary_cases(duckdb_connection):
    """
    Given:
        Adjacent intervals (touching but not overlapping)
    When:
        NEAREST operator is applied
    Then:
        Adjacent interval is reported as nearest (distance = 0)
    """
    intervals_a = [
        GenomicInterval("chr1", 100, 200, "a1", 100, "+"),
    ]
    intervals_b = [
        GenomicInterval("chr1", 200, 300, "b1", 150, "+"),
        GenomicInterval("chr1", 500, 600, "b2", 200, "+"),
    ]

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

    sql = transpile(
        """
        SELECT a.*, b.*
        FROM intervals_a a
        CROSS JOIN LATERAL NEAREST(
            intervals_b,
            reference := a.interval,
            k := 1
        ) b
        """,
        tables=["intervals_a", "intervals_b"],
    )
    giql_result = duckdb_connection.execute(sql).fetchall()

    assert len(giql_result) == 1
    # b1 is adjacent (distance 0), b2 is far (distance 300)
    assert giql_result[0][9] == "b1"
