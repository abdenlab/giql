"""Integration tests for GIQL CLUSTER operator.

These tests validate that GIQL's CLUSTER operator correctly assigns
cluster IDs to overlapping intervals. Since bedtools has no direct
CLUSTER equivalent, we cross-validate against bedtools merge: the
number of distinct clusters should equal the number of merged intervals.
"""

import pytest

from giql import transpile

from .utils.bedtools_wrapper import merge
from .utils.data_models import GenomicInterval
from .utils.duckdb_loader import load_intervals

pytestmark = pytest.mark.integration


def test_cluster_basic(duckdb_connection):
    """
    GIVEN a set of intervals with two overlapping groups
    WHEN CLUSTER operator is applied via GIQL
    THEN overlapping intervals share cluster IDs and distinct
    cluster count matches bedtools merge
    """
    intervals = [
        GenomicInterval("chr1", 100, 200, "i1", 100, "+"),
        GenomicInterval("chr1", 150, 250, "i2", 150, "+"),
        GenomicInterval("chr1", 400, 500, "i3", 200, "+"),
    ]

    load_intervals(
        duckdb_connection,
        "intervals",
        [i.to_tuple() for i in intervals],
    )

    bedtools_merged = merge([i.to_tuple() for i in intervals])

    sql = transpile(
        """
        SELECT *, CLUSTER(interval) AS cluster_id
        FROM intervals
        """,
        tables=["intervals"],
    )
    giql_result = duckdb_connection.execute(sql).fetchall()

    assert len(giql_result) == len(intervals)

    cluster_ids = {row[-1] for row in giql_result}
    assert len(cluster_ids) == len(bedtools_merged), (
        f"Expected {len(bedtools_merged)} clusters (matching merge "
        f"count), got {len(cluster_ids)}"
    )

    # i1 and i2 overlap, so they should share a cluster_id
    i1_cluster = giql_result[0][-1]
    i2_cluster = giql_result[1][-1]
    i3_cluster = giql_result[2][-1]
    assert i1_cluster == i2_cluster, (
        "Overlapping intervals i1 and i2 should share cluster_id"
    )
    assert i3_cluster != i1_cluster, (
        "Separated interval i3 should have a different cluster_id"
    )


def test_cluster_separated(duckdb_connection):
    """
    GIVEN non-overlapping intervals with gaps
    WHEN CLUSTER operator is applied
    THEN each interval gets its own cluster_id, matching bedtools merge count
    """
    intervals = [
        GenomicInterval("chr1", 100, 200, "i1", 100, "+"),
        GenomicInterval("chr1", 300, 400, "i2", 150, "+"),
        GenomicInterval("chr1", 500, 600, "i3", 200, "+"),
    ]

    load_intervals(
        duckdb_connection,
        "intervals",
        [i.to_tuple() for i in intervals],
    )

    bedtools_merged = merge([i.to_tuple() for i in intervals])

    sql = transpile(
        """
        SELECT *, CLUSTER(interval) AS cluster_id
        FROM intervals
        """,
        tables=["intervals"],
    )
    giql_result = duckdb_connection.execute(sql).fetchall()

    cluster_ids = {row[-1] for row in giql_result}
    assert len(cluster_ids) == len(bedtools_merged) == 3


def test_cluster_multiple_chromosomes(duckdb_connection):
    """
    GIVEN overlapping intervals on different chromosomes
    WHEN CLUSTER operator is applied
    THEN clustering occurs per chromosome independently, matching bedtools merge count
    """
    intervals = [
        GenomicInterval("chr1", 100, 200, "i1", 100, "+"),
        GenomicInterval("chr1", 150, 250, "i2", 150, "+"),
        GenomicInterval("chr2", 100, 200, "i3", 100, "+"),
        GenomicInterval("chr2", 150, 250, "i4", 150, "+"),
    ]

    load_intervals(
        duckdb_connection,
        "intervals",
        [i.to_tuple() for i in intervals],
    )

    bedtools_merged = merge([i.to_tuple() for i in intervals])

    sql = transpile(
        """
        SELECT *, CLUSTER(interval) AS cluster_id
        FROM intervals
        """,
        tables=["intervals"],
    )
    giql_result = duckdb_connection.execute(sql).fetchall()

    # Cluster IDs are per-partition (per-chrom), so we count
    # distinct (chrom, cluster_id) pairs
    chrom_clusters = {(row[0], row[-1]) for row in giql_result}
    assert len(chrom_clusters) == len(bedtools_merged) == 2


def test_cluster_stranded(duckdb_connection):
    """
    GIVEN overlapping intervals on different strands
    WHEN CLUSTER with stranded := true is applied
    THEN clustering occurs per strand independently, matching
    bedtools strand-aware merge count
    """
    intervals = [
        GenomicInterval("chr1", 100, 200, "i1", 100, "+"),
        GenomicInterval("chr1", 150, 250, "i2", 150, "+"),
        GenomicInterval("chr1", 120, 220, "i3", 200, "-"),
        GenomicInterval("chr1", 180, 280, "i4", 100, "-"),
    ]

    load_intervals(
        duckdb_connection,
        "intervals",
        [i.to_tuple() for i in intervals],
    )

    bedtools_strand_merged = merge(
        [i.to_tuple() for i in intervals],
        strand_mode="same",
    )

    sql = transpile(
        """
        SELECT *, CLUSTER(interval, stranded := true) AS cluster_id
        FROM intervals
        """,
        tables=["intervals"],
    )
    giql_result = duckdb_connection.execute(sql).fetchall()

    # Cluster IDs are per-partition (per-chrom-strand), so count
    # distinct (chrom, strand, cluster_id) pairs
    strand_clusters = {(row[0], row[5], row[-1]) for row in giql_result}
    assert len(strand_clusters) == len(bedtools_strand_merged)


def test_cluster_with_distance(giql_query):
    """
    GIVEN intervals with gaps of 50bp and 150bp
    WHEN CLUSTER with distance=100 is applied
    THEN gaps <= 100bp are in the same cluster and gaps > 100bp are separate
    """
    result = giql_query(
        """
        SELECT *, CLUSTER(interval, 100) AS cluster_id
        FROM intervals
        """,
        tables=["intervals"],
        intervals=[
            GenomicInterval("chr1", 100, 200, "i1", 100, "+"),
            GenomicInterval("chr1", 250, 350, "i2", 150, "+"),  # 50bp gap
            GenomicInterval("chr1", 500, 600, "i3", 200, "+"),  # 150bp gap
        ],
    )

    assert len(result) == 3

    # i1 and i2 should share a cluster (50bp gap <= 100)
    i1_cluster = result[0][-1]
    i2_cluster = result[1][-1]
    i3_cluster = result[2][-1]

    assert i1_cluster == i2_cluster, (
        "i1 and i2 (50bp gap) should be in the same cluster with distance=100"
    )
    assert i3_cluster != i1_cluster, (
        "i3 (150bp gap) should be in a different cluster with distance=100"
    )
