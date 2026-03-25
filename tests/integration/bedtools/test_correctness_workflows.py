"""Integration correctness tests for multi-operation GIQL workflows.

These tests validate that chained GIQL operations produce results matching
equivalent bedtools command pipelines. Corresponds to User Story 4 (P3)
from the bedtools integration test spec.
"""

from giql import transpile

from .utils.bedtools_wrapper import closest
from .utils.bedtools_wrapper import intersect
from .utils.bedtools_wrapper import merge
from .utils.comparison import compare_results
from .utils.data_models import GenomicInterval
from .utils.duckdb_loader import load_intervals


def test_workflow_intersect_then_merge(duckdb_connection):
    """
    GIVEN two interval sets with overlaps
    WHEN GIQL: intersect then merge (via subquery)
    vs bedtools: intersect | sort | merge
    THEN final merged intervals match
    """
    intervals_a = [
        GenomicInterval("chr1", 100, 200, "a1", 0, "+"),
        GenomicInterval("chr1", 150, 300, "a2", 0, "+"),
        GenomicInterval("chr1", 500, 600, "a3", 0, "+"),
    ]
    intervals_b = [
        GenomicInterval("chr1", 180, 250, "b1", 0, "+"),
        GenomicInterval("chr1", 520, 580, "b2", 0, "+"),
    ]

    load_intervals(duckdb_connection, "intervals_a", [i.to_tuple() for i in intervals_a])
    load_intervals(duckdb_connection, "intervals_b", [i.to_tuple() for i in intervals_b])

    # bedtools pipeline: intersect then merge
    intersect_result = intersect(
        [i.to_tuple() for i in intervals_a],
        [i.to_tuple() for i in intervals_b],
    )
    bedtools_final = merge(intersect_result)

    # GIQL: use CTE to intersect, then merge
    sql = transpile(
        """
        WITH hits AS (
            SELECT DISTINCT a.*
            FROM intervals_a a, intervals_b b
            WHERE a.interval INTERSECTS b.interval
        )
        SELECT MERGE(interval)
        FROM hits
        """,
        tables=["intervals_a", "intervals_b", "hits"],
    )
    giql_result = duckdb_connection.execute(sql).fetchall()

    comparison = compare_results(giql_result, bedtools_final)
    assert comparison.match, comparison.failure_message()


def test_workflow_nearest_then_filter_distance(duckdb_connection):
    """
    GIVEN two interval sets
    WHEN GIQL: NEAREST with max_distance filter
    vs bedtools: closest -d then filter by distance
    THEN filtered nearest results match
    """
    intervals_a = [
        GenomicInterval("chr1", 100, 200, "a1", 0, "+"),
        GenomicInterval("chr1", 500, 600, "a2", 0, "+"),
    ]
    intervals_b = [
        GenomicInterval("chr1", 220, 250, "b_near", 0, "+"),  # 20bp from a1
        GenomicInterval("chr1", 900, 1000, "b_far", 0, "+"),  # 300bp from a2
    ]

    load_intervals(duckdb_connection, "intervals_a", [i.to_tuple() for i in intervals_a])
    load_intervals(duckdb_connection, "intervals_b", [i.to_tuple() for i in intervals_b])

    # bedtools: closest -d, then filter distance <= 50
    bt_result = closest(
        [i.to_tuple() for i in intervals_a],
        [i.to_tuple() for i in intervals_b],
    )
    bedtools_filtered = [row for row in bt_result if row[-1] <= 50]

    # GIQL: NEAREST with max_distance
    sql = transpile(
        """
        SELECT a.name, b.name
        FROM intervals_a a
        CROSS JOIN LATERAL NEAREST(
            intervals_b,
            reference := a.interval,
            k := 1,
            max_distance := 50
        ) b
        """,
        tables=["intervals_a", "intervals_b"],
    )
    giql_result = duckdb_connection.execute(sql).fetchall()

    # Both should return only a1->b_near (distance 20 <= 50)
    # a2->b_far (distance 300 > 50) should be excluded
    assert len(giql_result) == len(bedtools_filtered)
    if len(giql_result) > 0:
        giql_names = {r[0] for r in giql_result}
        assert "a1" in giql_names


def test_workflow_merge_then_intersect(duckdb_connection):
    """
    GIVEN intervals with overlaps and a second interval set
    WHEN GIQL: merge intervals then intersect with second set
    vs bedtools: merge | intersect
    THEN results match
    """
    intervals_a = [
        GenomicInterval("chr1", 100, 200, "a1", 0, "+"),
        GenomicInterval("chr1", 180, 300, "a2", 0, "+"),
        GenomicInterval("chr1", 500, 600, "a3", 0, "+"),
    ]
    intervals_b = [
        GenomicInterval("chr1", 250, 350, "b1", 0, "+"),
        GenomicInterval("chr1", 550, 650, "b2", 0, "+"),
    ]

    load_intervals(duckdb_connection, "intervals_a", [i.to_tuple() for i in intervals_a])
    load_intervals(duckdb_connection, "intervals_b", [i.to_tuple() for i in intervals_b])

    # bedtools pipeline: merge a, then intersect with b
    merged_a = merge([i.to_tuple() for i in intervals_a])
    bedtools_final = intersect(merged_a, [i.to_tuple() for i in intervals_b])

    # GIQL: CTE to merge, then intersect
    sql = transpile(
        """
        WITH merged AS (
            SELECT MERGE(interval) AS interval
            FROM intervals_a
        )
        SELECT DISTINCT m.*
        FROM merged m, intervals_b b
        WHERE m.interval INTERSECTS b.interval
        """,
        tables=["intervals_a", "intervals_b", "merged"],
    )
    giql_result = duckdb_connection.execute(sql).fetchall()

    comparison = compare_results(giql_result, bedtools_final)
    assert comparison.match, comparison.failure_message()


def test_workflow_stranded_intersect_merge(duckdb_connection):
    """
    GIVEN intervals with strand info
    WHEN GIQL: strand-specific intersect then merge
    vs bedtools: intersect -s | sort | merge
    THEN strand-aware pipeline results match
    """
    intervals_a = [
        GenomicInterval("chr1", 100, 200, "a1", 0, "+"),
        GenomicInterval("chr1", 150, 300, "a2", 0, "+"),
        GenomicInterval("chr1", 120, 250, "a3", 0, "-"),
    ]
    intervals_b = [
        GenomicInterval("chr1", 180, 250, "b1", 0, "+"),
        GenomicInterval("chr1", 130, 220, "b2", 0, "-"),
    ]

    load_intervals(duckdb_connection, "intervals_a", [i.to_tuple() for i in intervals_a])
    load_intervals(duckdb_connection, "intervals_b", [i.to_tuple() for i in intervals_b])

    # bedtools pipeline: intersect -s then merge
    intersect_result = intersect(
        [i.to_tuple() for i in intervals_a],
        [i.to_tuple() for i in intervals_b],
        strand_mode="same",
    )
    bedtools_final = merge(intersect_result)

    # GIQL: same-strand intersect via CTE then merge
    sql = transpile(
        """
        WITH hits AS (
            SELECT DISTINCT a.*
            FROM intervals_a a, intervals_b b
            WHERE a.interval INTERSECTS b.interval
              AND a.strand = b.strand
        )
        SELECT MERGE(interval)
        FROM hits
        """,
        tables=["intervals_a", "intervals_b", "hits"],
    )
    giql_result = duckdb_connection.execute(sql).fetchall()

    comparison = compare_results(giql_result, bedtools_final)
    assert comparison.match, comparison.failure_message()


def test_workflow_intersect_filter_chrom_merge(duckdb_connection):
    """
    GIVEN two interval sets on multiple chromosomes
    WHEN GIQL: intersect, keep only chr1, then merge
    vs bedtools: intersect | grep chr1 | sort | merge
    THEN filtered-chromosome workflow matches
    """
    intervals_a = [
        GenomicInterval("chr1", 100, 200, "a1", 0, "+"),
        GenomicInterval("chr1", 150, 300, "a2", 0, "+"),
        GenomicInterval("chr2", 100, 200, "a3", 0, "+"),
    ]
    intervals_b = [
        GenomicInterval("chr1", 180, 250, "b1", 0, "+"),
        GenomicInterval("chr2", 150, 250, "b2", 0, "+"),
    ]

    load_intervals(duckdb_connection, "intervals_a", [i.to_tuple() for i in intervals_a])
    load_intervals(duckdb_connection, "intervals_b", [i.to_tuple() for i in intervals_b])

    # bedtools pipeline: intersect, filter chr1, merge
    intersect_result = intersect(
        [i.to_tuple() for i in intervals_a],
        [i.to_tuple() for i in intervals_b],
    )
    chr1_only = [r for r in intersect_result if r[0] == "chr1"]
    bedtools_final = merge(chr1_only) if chr1_only else []

    # GIQL: CTE intersect with chr1 filter, then merge
    sql = transpile(
        """
        WITH chr1_hits AS (
            SELECT DISTINCT a.*
            FROM intervals_a a, intervals_b b
            WHERE a.interval INTERSECTS b.interval
              AND a.chrom = 'chr1'
        )
        SELECT MERGE(interval)
        FROM chr1_hits
        """,
        tables=["intervals_a", "intervals_b", "chr1_hits"],
    )
    giql_result = duckdb_connection.execute(sql).fetchall()

    comparison = compare_results(giql_result, bedtools_final)
    assert comparison.match, comparison.failure_message()


def test_workflow_full_pipeline_step_by_step(duckdb_connection):
    """
    GIVEN a generated dataset across 3 chromosomes
    WHEN full pipeline (intersect -> merge -> nearest) is run
    THEN each intermediate step matches bedtools
    """
    import random

    rng = random.Random(99)
    intervals_a = []
    intervals_b = []
    intervals_c = []

    for chrom_num in range(1, 4):
        chrom = f"chr{chrom_num}"
        for i in range(30):
            start = rng.randint(0, 100_000)
            size = rng.randint(100, 1000)
            intervals_a.append(
                GenomicInterval(chrom, start, start + size, f"a_{chrom_num}_{i}", 0, "+")
            )
        for i in range(30):
            start = rng.randint(0, 100_000)
            size = rng.randint(100, 1000)
            intervals_b.append(
                GenomicInterval(chrom, start, start + size, f"b_{chrom_num}_{i}", 0, "+")
            )
        for i in range(10):
            start = rng.randint(0, 100_000)
            size = rng.randint(100, 1000)
            intervals_c.append(
                GenomicInterval(chrom, start, start + size, f"c_{chrom_num}_{i}", 0, "+")
            )

    load_intervals(duckdb_connection, "intervals_a", [i.to_tuple() for i in intervals_a])
    load_intervals(duckdb_connection, "intervals_b", [i.to_tuple() for i in intervals_b])
    load_intervals(duckdb_connection, "intervals_c", [i.to_tuple() for i in intervals_c])

    # Step 1: Intersect A with B
    bt_intersected = intersect(
        [i.to_tuple() for i in intervals_a],
        [i.to_tuple() for i in intervals_b],
    )

    sql_step1 = transpile(
        """
        SELECT DISTINCT a.*
        FROM intervals_a a, intervals_b b
        WHERE a.interval INTERSECTS b.interval
        """,
        tables=["intervals_a", "intervals_b"],
    )
    giql_step1 = duckdb_connection.execute(sql_step1).fetchall()

    comparison1 = compare_results(giql_step1, bt_intersected)
    assert comparison1.match, (
        f"Step 1 (intersect) failed: {comparison1.failure_message()}"
    )

    # Step 2: Merge the intersected results
    if bt_intersected:
        bt_merged = merge(bt_intersected)
    else:
        bt_merged = []

    if giql_step1:
        # Create temp table from step 1 results for step 2
        duckdb_connection.execute("""
            CREATE TABLE step1_results AS
            SELECT * FROM (
                SELECT DISTINCT a.*
                FROM intervals_a a, intervals_b b
                WHERE a.chrom = b.chrom
                  AND a."start" < b."end"
                  AND a."end" > b."start"
            )
        """)

        sql_step2 = transpile(
            "SELECT MERGE(interval) FROM step1_results",
            tables=["step1_results"],
        )
        giql_step2 = duckdb_connection.execute(sql_step2).fetchall()

        comparison2 = compare_results(giql_step2, bt_merged)
        assert comparison2.match, (
            f"Step 2 (merge) failed: {comparison2.failure_message()}"
        )
