"""Extended correctness tests for GIQL INTERSECTS operator vs bedtools intersect.

These tests cover boundary cases, scale, and edge scenarios beyond the basic
tests in test_intersect.py, ensuring comprehensive GIQL/bedtools equivalence.
"""

from giql import transpile

from .utils.bedtools_wrapper import intersect
from .utils.comparison import compare_results
from .utils.data_models import GenomicInterval
from .utils.duckdb_loader import load_intervals


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


def test_intersect_single_bp_overlap(duckdb_connection):
    """
    GIVEN two intervals overlapping by exactly 1bp
    WHEN GIQL INTERSECTS is compared to bedtools intersect
    THEN both detect the 1bp overlap
    """
    a = [GenomicInterval("chr1", 100, 200, "a1", 0, "+")]
    b = [GenomicInterval("chr1", 199, 300, "b1", 0, "+")]
    comparison = _run_intersect_comparison(duckdb_connection, a, b)
    assert comparison.match, comparison.failure_message()
    assert comparison.giql_row_count == 1


def test_intersect_containment_a_contains_b(duckdb_connection):
    """
    GIVEN interval A fully contains interval B
    WHEN GIQL INTERSECTS is compared to bedtools intersect
    THEN A is reported as intersecting
    """
    a = [GenomicInterval("chr1", 100, 500, "a1", 0, "+")]
    b = [GenomicInterval("chr1", 200, 300, "b1", 0, "+")]
    comparison = _run_intersect_comparison(duckdb_connection, a, b)
    assert comparison.match, comparison.failure_message()
    assert comparison.giql_row_count == 1


def test_intersect_containment_b_contains_a(duckdb_connection):
    """
    GIVEN interval B fully contains interval A
    WHEN GIQL INTERSECTS is compared to bedtools intersect
    THEN A is reported as intersecting
    """
    a = [GenomicInterval("chr1", 200, 300, "a1", 0, "+")]
    b = [GenomicInterval("chr1", 100, 500, "b1", 0, "+")]
    comparison = _run_intersect_comparison(duckdb_connection, a, b)
    assert comparison.match, comparison.failure_message()
    assert comparison.giql_row_count == 1


def test_intersect_deduplication(duckdb_connection):
    """
    GIVEN one interval in A overlapping multiple intervals in B
    WHEN GIQL INTERSECTS with DISTINCT is compared to bedtools intersect -u
    THEN A interval reported once
    """
    a = [GenomicInterval("chr1", 100, 300, "a1", 0, "+")]
    b = [
        GenomicInterval("chr1", 150, 200, "b1", 0, "+"),
        GenomicInterval("chr1", 200, 250, "b2", 0, "+"),
        GenomicInterval("chr1", 250, 350, "b3", 0, "+"),
    ]
    comparison = _run_intersect_comparison(duckdb_connection, a, b)
    assert comparison.match, comparison.failure_message()
    assert comparison.giql_row_count == 1


def test_intersect_non_standard_chroms(duckdb_connection):
    """
    GIVEN intervals on non-standard chromosome names (chrM, chrUn)
    WHEN GIQL INTERSECTS is compared to bedtools intersect
    THEN results match regardless of chromosome naming
    """
    a = [
        GenomicInterval("chrM", 100, 200, "a1", 0, "+"),
        GenomicInterval("chrUn", 100, 200, "a2", 0, "+"),
    ]
    b = [
        GenomicInterval("chrM", 150, 250, "b1", 0, "+"),
        GenomicInterval("chrUn", 150, 250, "b2", 0, "+"),
    ]
    comparison = _run_intersect_comparison(duckdb_connection, a, b)
    assert comparison.match, comparison.failure_message()
    assert comparison.giql_row_count == 2


def test_intersect_large_intervals(duckdb_connection):
    """
    GIVEN very large genomic intervals (spanning millions of bases)
    WHEN GIQL INTERSECTS is compared to bedtools intersect
    THEN results match correctly
    """
    a = [GenomicInterval("chr1", 0, 10_000_000, "a1", 0, "+")]
    b = [GenomicInterval("chr1", 5_000_000, 15_000_000, "b1", 0, "+")]
    comparison = _run_intersect_comparison(duckdb_connection, a, b)
    assert comparison.match, comparison.failure_message()
    assert comparison.giql_row_count == 1


def test_intersect_many_intervals_scale(duckdb_connection):
    """
    GIVEN a generated dataset with 100 intervals per chromosome on 3 chromosomes
    WHEN GIQL INTERSECTS is compared to bedtools intersect
    THEN results match on the full dataset
    """
    import random

    rng = random.Random(42)
    intervals_a = []
    intervals_b = []

    for chrom_num in range(1, 4):
        chrom = f"chr{chrom_num}"
        for i in range(100):
            start = rng.randint(0, 900_000)
            size = rng.randint(100, 1000)
            strand = rng.choice(["+", "-"])
            intervals_a.append(
                GenomicInterval(
                    chrom,
                    start,
                    start + size,
                    f"a_{chrom_num}_{i}",
                    0,
                    strand,
                )
            )
            start = rng.randint(0, 900_000)
            size = rng.randint(100, 1000)
            strand = rng.choice(["+", "-"])
            intervals_b.append(
                GenomicInterval(
                    chrom,
                    start,
                    start + size,
                    f"b_{chrom_num}_{i}",
                    0,
                    strand,
                )
            )

    comparison = _run_intersect_comparison(duckdb_connection, intervals_a, intervals_b)
    assert comparison.match, comparison.failure_message()


def test_intersect_same_strand_correctness(duckdb_connection):
    """
    GIVEN overlapping intervals with mixed strands
    WHEN GIQL INTERSECTS with same-strand filter is compared to bedtools -s
    THEN only same-strand overlaps match
    """
    a = [
        GenomicInterval("chr1", 100, 200, "a_plus", 0, "+"),
        GenomicInterval("chr1", 100, 200, "a_minus", 0, "-"),
    ]
    b = [GenomicInterval("chr1", 150, 250, "b_plus", 0, "+")]
    comparison = _run_intersect_comparison(
        duckdb_connection,
        a,
        b,
        strand_filter="a.strand = b.strand",
    )
    assert comparison.match, comparison.failure_message()
    assert comparison.giql_row_count == 1


def test_intersect_opposite_strand_correctness(duckdb_connection):
    """
    GIVEN overlapping intervals with mixed strands
    WHEN GIQL INTERSECTS with opposite-strand filter is compared to bedtools -S
    THEN only opposite-strand overlaps match
    """
    a = [
        GenomicInterval("chr1", 100, 200, "a_plus", 0, "+"),
        GenomicInterval("chr1", 100, 200, "a_minus", 0, "-"),
    ]
    b = [GenomicInterval("chr1", 150, 250, "b_plus", 0, "+")]
    comparison = _run_intersect_comparison(
        duckdb_connection,
        a,
        b,
        strand_filter="a.strand != b.strand",
    )
    assert comparison.match, comparison.failure_message()
    assert comparison.giql_row_count == 1
