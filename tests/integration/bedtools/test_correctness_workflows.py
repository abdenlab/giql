"""Integration correctness tests for multi-operation GIQL workflows.

These tests validate that chained GIQL operations produce results matching
equivalent bedtools command pipelines. Corresponds to User Story 4 (P3)
from the bedtools integration test spec.
"""

import pytest

from giql import transpile

from .utils.bedtools_wrapper import closest
from .utils.bedtools_wrapper import intersect
from .utils.bedtools_wrapper import merge
from .utils.comparison import compare_results
from .utils.data_models import GenomicInterval
from .utils.duckdb_loader import load_intervals

pytestmark = pytest.mark.integration


def test_pipeline_should_match_bedtools_when_intersect_chained_into_merge(duckdb_connection):
    """Test that chaining intersect into merge in GIQL matches the bedtools pipeline.

    Given:
        Two interval sets with overlaps on chr1
    When:
        GIQL intersects via CTE and then merges, compared against
        bedtools intersect piped into bedtools merge
    Then:
        It should produce identical merged intervals
    """
    # Arrange
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

    # Act
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

    # Assert
    comparison = compare_results(giql_result, bedtools_final)
    assert comparison.match, comparison.failure_message()


def test_pipeline_should_filter_by_distance_when_nearest_max_distance_applied(duckdb_connection):
    """Test that NEAREST with max_distance matches bedtools closest filtered by distance.

    Given:
        Two interval sets where one B interval is within 50bp of an A
        interval and another is far beyond that threshold
    When:
        GIQL runs NEAREST with max_distance := 50, compared against
        bedtools closest -d post-filtered to distance <= 50
    Then:
        It should return only the close neighbor pair and drop the far one
    """
    # Arrange
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

    # Act
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

    # Assert
    # Both should return only a1->b_near (distance 20 <= 50)
    # a2->b_far (distance 300 > 50) should be excluded
    assert len(giql_result) == len(bedtools_filtered)
    if len(giql_result) > 0:
        giql_names = {r[0] for r in giql_result}
        assert "a1" in giql_names


def test_pipeline_should_match_bedtools_when_merge_chained_into_intersect(duckdb_connection):
    """Test that merging then intersecting in GIQL matches the bedtools pipeline.

    Given:
        An A interval set with overlapping intervals and a disjoint B set
    When:
        GIQL merges A via CTE then intersects with B, compared against
        bedtools merge of A piped into bedtools intersect against B
    Then:
        It should produce matching interval coordinates
    """
    # Arrange
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

    # Act
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

    # Assert
    # MERGE outputs BED3 (chrom, start, end); compare only coordinates
    bedtools_coords = [row[:3] for row in bedtools_final]
    comparison = compare_results(giql_result, bedtools_coords)
    assert comparison.match, comparison.failure_message()


def test_pipeline_should_preserve_strand_when_intersect_then_merge_stranded(duckdb_connection):
    """Test that strand-aware intersect chained into merge matches the bedtools pipeline.

    Given:
        Two interval sets carrying strand information, with mixed plus
        and minus strand overlaps
    When:
        GIQL performs a same-strand intersect via CTE then merges,
        compared against bedtools intersect -s piped into bedtools merge
    Then:
        It should produce matching merged intervals honoring strand
    """
    # Arrange
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

    # Act
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

    # Assert
    comparison = compare_results(giql_result, bedtools_final)
    assert comparison.match, comparison.failure_message()


def test_pipeline_should_match_bedtools_when_intersect_chrom_filter_then_merge(duckdb_connection):
    """Test that intersect followed by a chr1 filter and a merge matches the bedtools pipeline.

    Given:
        Two interval sets spanning chr1 and chr2 with overlaps on both
        chromosomes
    When:
        GIQL intersects with a chrom = 'chr1' predicate inside a CTE
        and then merges, compared against bedtools intersect, then a
        Python-side chr1 filter, then bedtools merge
    Then:
        It should produce matching merged intervals restricted to chr1
    """
    # Arrange
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

    # Act
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

    # Assert
    comparison = compare_results(giql_result, bedtools_final)
    assert comparison.match, comparison.failure_message()


def test_pipeline_should_match_bedtools_when_chained_step_by_step(
    duckdb_connection,
):
    """Test chained GIQL pipeline matches bedtools at each step.

    Given:
        Three interval sets across two chromosomes — A and B as inputs
        for intersect + merge, C as reference for nearest — hand-crafted
        so the pipeline output is unambiguous (no tie-breaking)
    When:
        Each GIQL step's output is materialized as a table and fed
        back into the next GIQL step, and each bedtools equivalent
        operates on its own prior step's output
    Then:
        GIQL and bedtools outputs should match at each of the three
        stages: full row equality for intersect and merge, and
        equal (a_name, b_name) neighbor pairs for nearest (distance
        values are compared in the dedicated nearest tests because
        bedtools 2.31+ uses the N+1 half-open gap convention)
    """
    # Arrange
    intervals_a = [
        GenomicInterval("chr1", 100, 300, "a1", 0, "+"),
        GenomicInterval("chr1", 500, 700, "a2", 0, "+"),
        GenomicInterval("chr2", 100, 300, "a3", 0, "+"),
    ]
    intervals_b = [
        GenomicInterval("chr1", 200, 400, "b1", 0, "+"),
        GenomicInterval("chr1", 600, 800, "b2", 0, "+"),
        GenomicInterval("chr2", 200, 400, "b3", 0, "+"),
    ]
    intervals_c = [
        GenomicInterval("chr1", 5000, 5100, "c1", 0, "+"),
        GenomicInterval("chr2", 5000, 5100, "c2", 0, "+"),
    ]
    load_intervals(
        duckdb_connection, "intervals_a", [i.to_tuple() for i in intervals_a]
    )
    load_intervals(
        duckdb_connection, "intervals_b", [i.to_tuple() for i in intervals_b]
    )
    load_intervals(
        duckdb_connection, "intervals_c", [i.to_tuple() for i in intervals_c]
    )

    # Act & Assert — Step 1: GIQL intersect vs bedtools intersect
    bt_step1 = intersect(
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
    c1 = compare_results(giql_step1, bt_step1)
    assert c1.match, f"Step 1 (intersect): {c1.failure_message()}"

    # Act & Assert — Step 2: materialize GIQL step-1 output, GIQL MERGE
    assert giql_step1, "fixture should produce at least one intersecting row"
    load_intervals(duckdb_connection, "step1_results", giql_step1)
    bt_step2 = merge(bt_step1)
    sql_step2 = transpile(
        "SELECT MERGE(interval) FROM step1_results",
        tables=["step1_results"],
    )
    giql_step2 = duckdb_connection.execute(sql_step2).fetchall()
    c2 = compare_results(giql_step2, bt_step2)
    assert c2.match, f"Step 2 (merge): {c2.failure_message()}"

    # Act & Assert — Step 3: pad BED3 step-2 output to BED6, GIQL NEAREST
    assert giql_step2, "step 2 should produce at least one merged interval"
    giql_step2_bed6 = [
        (row[0], row[1], row[2], f"step2_{i}", 0, "+")
        for i, row in enumerate(giql_step2)
    ]
    load_intervals(duckdb_connection, "step2_results", giql_step2_bed6)
    bt_step3 = closest(
        giql_step2_bed6, [i.to_tuple() for i in intervals_c]
    )
    sql_step3 = transpile(
        """
        SELECT a.*, b.*
        FROM step2_results a
        CROSS JOIN LATERAL NEAREST(intervals_c, reference := a.interval) b
        ORDER BY a.chrom, a.start
        """,
        tables=["step2_results", "intervals_c"],
    )
    giql_step3 = duckdb_connection.execute(sql_step3).fetchall()
    giql_pairs = {(row[3], row[9]) for row in giql_step3}
    bt_pairs = {(row[3], row[9]) for row in bt_step3}
    assert giql_pairs == bt_pairs, (
        f"Step 3 (nearest) neighbor pairs differ\n"
        f"  GIQL: {sorted(giql_pairs)}\n"
        f"  bedtools: {sorted(bt_pairs)}"
    )
