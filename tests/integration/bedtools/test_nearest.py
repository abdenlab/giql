"""Integration tests for GIQL NEAREST operator.

These tests validate that GIQL's NEAREST operator produces results
consistent with bedtools closest command.
"""

from giql import transpile

from .utils.bedtools_wrapper import closest
from .utils.data_models import GenomicInterval
from .utils.duckdb_loader import load_intervals


def test_nearest_non_overlapping(duckdb_connection):
    """
    GIVEN two sets of non-overlapping intervals
    WHEN NEAREST operator is applied
    THEN each interval in A finds its closest neighbor in B, matching bedtools closest
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
        SELECT a.name, b.name, b.distance
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

    # Align both sides by A name; the bedtools closest row is BED6+BED6+distance,
    # so A name is field 3, B name field 9, and distance field 12.
    giql_by_a = {a_name: (b_name, distance) for a_name, b_name, distance in giql_result}
    bt_by_a = {row[3]: (row[9], row[12]) for row in bedtools_result}

    assert giql_by_a.keys() == bt_by_a.keys(), (
        f"A coverage mismatch: GIQL={set(giql_by_a)}, bedtools={set(bt_by_a)}"
    )

    for a_name, (bt_b_name, bt_distance) in bt_by_a.items():
        giql_b_name, giql_distance = giql_by_a[a_name]
        # Matched B neighbor must agree.
        assert giql_b_name == bt_b_name, (
            f"{a_name} B mismatch: GIQL={giql_b_name}, bedtools={bt_b_name}"
        )
        # GIQL distance is signed; its magnitude must equal bedtools closest -d.
        assert abs(giql_distance) == bt_distance, (
            f"{a_name} distance mismatch: GIQL={giql_distance}, bedtools={bt_distance}"
        )


def test_nearest_multiple_candidates(giql_query):
    """
    GIVEN an interval in A with multiple equidistant intervals in B
    WHEN NEAREST operator is applied
    THEN one of the equidistant intervals is returned
    """
    result = giql_query(
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
        intervals_a=[GenomicInterval("chr1", 300, 400, "a1", 100, "+")],
        intervals_b=[
            GenomicInterval("chr1", 100, 200, "b1", 100, "+"),
            GenomicInterval("chr1", 500, 600, "b2", 150, "+"),
        ],
    )

    assert len(result) == 1
    assert result[0][3] == "a1"
    # Nearest could be either b1 or b2 (both equidistant at 100bp)
    assert result[0][9] in ("b1", "b2")


def test_nearest_cross_chromosome(duckdb_connection):
    """
    GIVEN intervals on different chromosomes
    WHEN NEAREST operator is applied
    THEN each interval finds nearest only on same chromosome,
    matching bedtools closest
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


def test_nearest_boundary_cases(giql_query):
    """
    GIVEN adjacent intervals (touching but not overlapping)
    WHEN NEAREST operator is applied
    THEN adjacent interval is reported as nearest (distance = 1, bedtools parity)
    """
    result = giql_query(
        """
        SELECT b.name, b.distance
        FROM intervals_a a
        CROSS JOIN LATERAL NEAREST(
            intervals_b,
            reference := a.interval,
            k := 1
        ) b
        """,
        tables=["intervals_a", "intervals_b"],
        intervals_a=[GenomicInterval("chr1", 100, 200, "a1", 100, "+")],
        intervals_b=[
            GenomicInterval("chr1", 200, 300, "b1", 150, "+"),
            GenomicInterval("chr1", 500, 600, "b2", 200, "+"),
        ],
    )

    assert len(result) == 1
    # b1 is adjacent (book-ended), reported as distance 1 under bedtools parity;
    # b2 is far (301).
    assert result[0][0] == "b1"
    assert result[0][1] == 1


def test_nearest_k_greater_than_one(giql_query):
    """
    GIVEN one query interval and three database intervals at different distances
    WHEN NEAREST with k := 3 is applied
    THEN all 3 neighbors are returned, ordered by distance
    """
    result = giql_query(
        """
        SELECT a.name, b.name
        FROM intervals_a a
        CROSS JOIN LATERAL NEAREST(
            intervals_b,
            reference := a.interval,
            k := 3
        ) b
        """,
        tables=["intervals_a", "intervals_b"],
        intervals_a=[GenomicInterval("chr1", 200, 300, "a1", 100, "+")],
        intervals_b=[
            GenomicInterval("chr1", 100, 150, "b_far", 100, "+"),
            GenomicInterval("chr1", 310, 350, "b_near", 100, "+"),
            GenomicInterval("chr1", 500, 600, "b_farther", 100, "+"),
        ],
    )

    b_names = [row[1] for row in result]
    assert len(b_names) == 3, f"Expected 3 results for k=3, got {len(b_names)}"
    assert set(b_names) == {"b_far", "b_near", "b_farther"}


def test_nearest_k3_reported_distance_matches_bedtools(duckdb_connection):
    """
    GIVEN one query interval and three candidates at distinct (untied) gaps
    WHEN NEAREST with k := 3 is applied and the distance column is read
    THEN each neighbor's reported distance magnitude equals live bedtools
        closest -d -k 3 for the same neighbor
    """
    intervals_a = [GenomicInterval("chr1", 200, 300, "a1", 100, "+")]
    intervals_b = [
        GenomicInterval("chr1", 310, 350, "b_near", 100, "+"),  # gap 10 -> 11
        GenomicInterval("chr1", 400, 450, "b_mid", 100, "+"),  # gap 100 -> 101
        GenomicInterval("chr1", 600, 700, "b_far", 100, "+"),  # gap 300 -> 301
    ]

    load_intervals(duckdb_connection, "intervals_a", [i.to_tuple() for i in intervals_a])
    load_intervals(duckdb_connection, "intervals_b", [i.to_tuple() for i in intervals_b])

    bedtools_result = closest(
        [i.to_tuple() for i in intervals_a],
        [i.to_tuple() for i in intervals_b],
        k=3,
    )
    oracle_by_name = {row[9]: row[12] for row in bedtools_result}

    sql = transpile(
        """
        SELECT b.name, b.distance
        FROM intervals_a a
        CROSS JOIN LATERAL NEAREST(
            intervals_b,
            reference := a.interval,
            k := 3
        ) b
        """,
        tables=["intervals_a", "intervals_b"],
    )
    giql_result = duckdb_connection.execute(sql).fetchall()

    assert len(giql_result) == 3
    for name, distance in giql_result:
        assert abs(distance) == oracle_by_name[name], (
            f"{name}: GIQL={distance}, bedtools={oracle_by_name[name]}"
        )


def test_nearest_k_exceeds_available(giql_query):
    """
    GIVEN one query interval and only two database intervals
    WHEN NEAREST with k := 5 is applied
    THEN only 2 rows are returned (fewer than k available)
    """
    result = giql_query(
        """
        SELECT a.name, b.name
        FROM intervals_a a
        CROSS JOIN LATERAL NEAREST(
            intervals_b,
            reference := a.interval,
            k := 5
        ) b
        """,
        tables=["intervals_a", "intervals_b"],
        intervals_a=[GenomicInterval("chr1", 200, 300, "a1", 100, "+")],
        intervals_b=[
            GenomicInterval("chr1", 100, 150, "b1", 100, "+"),
            GenomicInterval("chr1", 400, 500, "b2", 100, "+"),
        ],
    )

    assert len(result) == 2, f"Expected 2 results (fewer than k=5), got {len(result)}"


def test_nearest_max_distance(giql_query):
    """
    GIVEN candidates whose reported distances straddle max_distance := 50: one
        below, one exactly at 50, and one at 51
    WHEN NEAREST with max_distance := 50 is applied
    THEN the candidates at 11 and exactly 50 are kept and the one at 51 is
        excluded, pinning the inclusive `<= max_distance` boundary on both sides
    """
    result = giql_query(
        """
        SELECT b.name, b.distance
        FROM intervals_a a
        CROSS JOIN LATERAL NEAREST(
            intervals_b,
            reference := a.interval,
            k := 5,
            max_distance := 50
        ) b
        """,
        tables=["intervals_a", "intervals_b"],
        intervals_a=[GenomicInterval("chr1", 200, 300, "a1", 100, "+")],
        intervals_b=[
            # gap 10 -> reported 11 (well within)
            GenomicInterval("chr1", 310, 350, "b_near", 100, "+"),
            # gap 49 -> reported 50 == max_distance (inclusive, kept)
            GenomicInterval("chr1", 349, 400, "b_edge", 100, "+"),
            # gap 50 -> reported 51 == max_distance + 1 (excluded)
            GenomicInterval("chr1", 350, 450, "b_over", 100, "+"),
        ],
    )

    by_name = {row[0]: row[1] for row in result}
    assert sorted(by_name) == ["b_edge", "b_near"], (
        f"Expected b_near and b_edge within 50bp, got {sorted(by_name)}"
    )
    assert by_name["b_near"] == 11
    assert by_name["b_edge"] == 50  # exactly at the boundary, kept


def test_nearest_standalone_literal_reference(giql_query):
    """
    GIVEN a table with intervals
    WHEN NEAREST is used in standalone mode with a literal reference
    THEN the nearest intervals to the literal position are returned
    """
    result = giql_query(
        """
        SELECT *
        FROM NEAREST(
            intervals,
            reference := 'chr1:350-360',
            k := 2
        )
        """,
        tables=["intervals"],
        intervals=[
            GenomicInterval("chr1", 100, 200, "near", 100, "+"),
            GenomicInterval("chr1", 400, 500, "mid", 100, "+"),
            GenomicInterval("chr1", 800, 900, "far", 100, "+"),
        ],
    )

    names = [row[3] for row in result]
    assert len(names) == 2, f"Expected 2 results for k=2, got {len(names)}"
    # near is 150bp away, mid is 40bp away, far is 440bp away
    assert set(names) == {"near", "mid"}
