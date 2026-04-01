"""Property-based correctness tests for INTERSECTS binned equi-join.

These tests use hypothesis to generate random genomic intervals of
varying sizes — including intervals that span multiple bins — and
verify that GIQL's binned equi-join produces identical results to
bedtools intersect.
"""

from hypothesis import HealthCheck
from hypothesis import given
from hypothesis import settings
from hypothesis import strategies as st

from giql import transpile

from .utils.bedtools_wrapper import intersect
from .utils.comparison import compare_results
from .utils.data_models import GenomicInterval
from .utils.duckdb_loader import load_intervals

duckdb = __import__("pytest").importorskip("duckdb")


# ---------------------------------------------------------------------------
# Strategies
# ---------------------------------------------------------------------------

CHROMS = ["chr1", "chr2", "chr3"]


@st.composite
def genomic_interval_st(draw, idx=None):
    """Generate a random GenomicInterval that can span multiple 10k bins."""
    chrom = draw(st.sampled_from(CHROMS))
    start = draw(st.integers(min_value=0, max_value=1_000_000))
    length = draw(st.integers(min_value=1, max_value=200_000))
    score = draw(st.integers(min_value=0, max_value=1000))
    strand = draw(st.sampled_from(["+", "-"]))
    # Name is set by the list strategy to guarantee uniqueness, avoiding
    # the known DISTINCT duplicate-collapse limitation.
    name = (
        f"r{idx}"
        if idx is not None
        else draw(st.from_regex(r"r[0-9]{1,6}", fullmatch=True))
    )
    return GenomicInterval(chrom, start, start + length, name, score, strand)


@st.composite
def unique_interval_list_st(draw, max_size=60):
    """Generate a list of intervals with unique names."""
    n = draw(st.integers(min_value=1, max_value=max_size))
    intervals = []
    for i in range(n):
        iv = draw(genomic_interval_st(idx=i))
        intervals.append(iv)
    return intervals


interval_list_st = unique_interval_list_st()


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _run_giql(intervals_a, intervals_b):
    """Run the binned-join INTERSECTS query via DuckDB and return result rows."""
    conn = duckdb.connect(":memory:")
    try:
        load_intervals(conn, "intervals_a", [i.to_tuple() for i in intervals_a])
        load_intervals(conn, "intervals_b", [i.to_tuple() for i in intervals_b])

        sql = transpile(
            """
            SELECT DISTINCT a.*
            FROM intervals_a a, intervals_b b
            WHERE a.interval INTERSECTS b.interval
            """,
            tables=["intervals_a", "intervals_b"],
        )
        return conn.execute(sql).fetchall()
    finally:
        conn.close()


def _run_bedtools(intervals_a, intervals_b):
    """Run bedtools intersect -u and return result tuples."""
    return intersect(
        [i.to_tuple() for i in intervals_a],
        [i.to_tuple() for i in intervals_b],
    )


# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------


@given(intervals_a=interval_list_st, intervals_b=interval_list_st)
@settings(
    max_examples=50,
    deadline=None,
    suppress_health_check=[HealthCheck.too_slow],
)
def test_binned_join_matches_bedtools(intervals_a, intervals_b):
    """
    GIVEN two randomly generated sets of genomic intervals
    WHEN GIQL INTERSECTS binned equi-join is executed
    THEN results match bedtools intersect -u exactly
    """
    giql_result = _run_giql(intervals_a, intervals_b)
    bedtools_result = _run_bedtools(intervals_a, intervals_b)

    comparison = compare_results(giql_result, bedtools_result)
    assert comparison.match, comparison.failure_message()


@given(intervals=interval_list_st)
@settings(
    max_examples=30,
    deadline=None,
    suppress_health_check=[HealthCheck.too_slow],
)
def test_self_join_matches_bedtools(intervals):
    """
    GIVEN a randomly generated set of genomic intervals
    WHEN GIQL INTERSECTS self-join is executed
    THEN results match bedtools intersect -u with the same file as A and B
    """
    giql_result = _run_giql(intervals, intervals)
    bedtools_result = _run_bedtools(intervals, intervals)

    comparison = compare_results(giql_result, bedtools_result)
    assert comparison.match, comparison.failure_message()


@given(
    intervals_a=unique_interval_list_st(max_size=30),
    intervals_b=unique_interval_list_st(max_size=30),
    bin_size=st.sampled_from([100, 1_000, 10_000, 100_000]),
)
@settings(
    max_examples=40,
    deadline=None,
    suppress_health_check=[HealthCheck.too_slow],
)
def test_bin_size_does_not_affect_correctness(intervals_a, intervals_b, bin_size):
    """
    GIVEN two randomly generated sets of genomic intervals and a bin size
    WHEN GIQL INTERSECTS is executed with that bin size
    THEN results match bedtools intersect -u regardless of bin size
    """
    conn = duckdb.connect(":memory:")
    try:
        load_intervals(conn, "intervals_a", [i.to_tuple() for i in intervals_a])
        load_intervals(conn, "intervals_b", [i.to_tuple() for i in intervals_b])

        sql = transpile(
            """
            SELECT DISTINCT a.*
            FROM intervals_a a, intervals_b b
            WHERE a.interval INTERSECTS b.interval
            """,
            tables=["intervals_a", "intervals_b"],
            bin_size=bin_size,
        )
        giql_result = conn.execute(sql).fetchall()
    finally:
        conn.close()

    bedtools_result = _run_bedtools(intervals_a, intervals_b)

    comparison = compare_results(giql_result, bedtools_result)
    assert comparison.match, f"bin_size={bin_size}: {comparison.failure_message()}"


@given(
    intervals_a=unique_interval_list_st(max_size=8),
    intervals_b=unique_interval_list_st(max_size=8),
    intervals_c=unique_interval_list_st(max_size=8),
)
@settings(
    max_examples=40,
    deadline=None,
    suppress_health_check=[HealthCheck.too_slow],
)
def test_multi_table_join_matches_bedtools(intervals_a, intervals_b, intervals_c):
    """
    GIVEN three randomly generated sets of genomic intervals
    WHEN GIQL three-way INTERSECTS join is executed
    THEN the A-side rows match bedtools intersect chained A->B then ->C
    """
    conn = duckdb.connect(":memory:")
    try:
        load_intervals(conn, "intervals_a", [i.to_tuple() for i in intervals_a])
        load_intervals(conn, "intervals_b", [i.to_tuple() for i in intervals_b])
        load_intervals(conn, "intervals_c", [i.to_tuple() for i in intervals_c])

        sql = transpile(
            """
            SELECT DISTINCT a.*
            FROM intervals_a a
            JOIN intervals_b b ON a.interval INTERSECTS b.interval
            JOIN intervals_c c ON a.interval INTERSECTS c.interval
            """,
            tables=["intervals_a", "intervals_b", "intervals_c"],
        )
        giql_result = conn.execute(sql).fetchall()
    finally:
        conn.close()

    # bedtools equivalent: chain A∩B then filter against C
    tuples_a = [i.to_tuple() for i in intervals_a]
    tuples_b = [i.to_tuple() for i in intervals_b]
    tuples_c = [i.to_tuple() for i in intervals_c]
    ab_result = intersect(tuples_a, tuples_b)
    if ab_result:
        bedtools_result = intersect(ab_result, tuples_c)
    else:
        bedtools_result = []

    comparison = compare_results(giql_result, bedtools_result)
    assert comparison.match, comparison.failure_message()


@given(
    intervals_a=unique_interval_list_st(max_size=30),
    intervals_b=unique_interval_list_st(max_size=30),
)
@settings(
    max_examples=40,
    deadline=None,
    suppress_health_check=[HealthCheck.too_slow],
)
def test_left_join_matches_bedtools_loj(intervals_a, intervals_b):
    """
    GIVEN two randomly generated sets of genomic intervals
    WHEN GIQL LEFT JOIN INTERSECTS is executed
    THEN results match bedtools intersect -loj exactly
    """
    conn = duckdb.connect(":memory:")
    try:
        load_intervals(conn, "intervals_a", [i.to_tuple() for i in intervals_a])
        load_intervals(conn, "intervals_b", [i.to_tuple() for i in intervals_b])

        sql = transpile(
            """
            SELECT DISTINCT
                a.chrom, a.start, a.end, a.name, a.score, a.strand,
                b.chrom AS b_chrom, b.start AS b_start, b.end AS b_end,
                b.name AS b_name, b.score AS b_score, b.strand AS b_strand
            FROM intervals_a a
            LEFT JOIN intervals_b b ON a.interval INTERSECTS b.interval
            """,
            tables=["intervals_a", "intervals_b"],
        )
        giql_result = conn.execute(sql).fetchall()
    finally:
        conn.close()

    bedtools_result = intersect(
        [i.to_tuple() for i in intervals_a],
        [i.to_tuple() for i in intervals_b],
        loj=True,
    )

    comparison = compare_results(giql_result, bedtools_result)
    assert comparison.match, comparison.failure_message()


# ---------------------------------------------------------------------------
# -v (inverse / anti-join)
# ---------------------------------------------------------------------------


@given(
    intervals_a=unique_interval_list_st(max_size=30),
    intervals_b=unique_interval_list_st(max_size=30),
)
@settings(
    max_examples=40,
    deadline=None,
    suppress_health_check=[HealthCheck.too_slow],
)
def test_inverse_matches_bedtools_v(intervals_a, intervals_b):
    """
    GIVEN two randomly generated sets of genomic intervals
    WHEN GIQL anti-join (LEFT JOIN WHERE b IS NULL) is executed
    THEN results match bedtools intersect -v exactly
    """
    conn = duckdb.connect(":memory:")
    try:
        load_intervals(conn, "intervals_a", [i.to_tuple() for i in intervals_a])
        load_intervals(conn, "intervals_b", [i.to_tuple() for i in intervals_b])

        sql = transpile(
            """
            SELECT DISTINCT a.*
            FROM intervals_a a
            LEFT JOIN intervals_b b ON a.interval INTERSECTS b.interval
            WHERE b.chrom IS NULL
            """,
            tables=["intervals_a", "intervals_b"],
        )
        giql_result = conn.execute(sql).fetchall()
    finally:
        conn.close()

    bedtools_result = intersect(
        [i.to_tuple() for i in intervals_a],
        [i.to_tuple() for i in intervals_b],
        inverse=True,
    )

    comparison = compare_results(giql_result, bedtools_result)
    assert comparison.match, comparison.failure_message()


# ---------------------------------------------------------------------------
# -wa -wb (write both A and B entries)
# ---------------------------------------------------------------------------


@given(
    intervals_a=unique_interval_list_st(max_size=20),
    intervals_b=unique_interval_list_st(max_size=20),
)
@settings(
    max_examples=40,
    deadline=None,
    suppress_health_check=[HealthCheck.too_slow],
)
def test_write_both_matches_bedtools_wa_wb(intervals_a, intervals_b):
    """
    GIVEN two randomly generated sets of genomic intervals
    WHEN GIQL INTERSECTS join selecting both sides is executed
    THEN results match bedtools intersect -wa -wb exactly
    """
    conn = duckdb.connect(":memory:")
    try:
        load_intervals(conn, "intervals_a", [i.to_tuple() for i in intervals_a])
        load_intervals(conn, "intervals_b", [i.to_tuple() for i in intervals_b])

        sql = transpile(
            """
            SELECT
                a.chrom, a.start, a.end, a.name, a.score, a.strand,
                b.chrom AS b_chrom, b.start AS b_start, b.end AS b_end,
                b.name AS b_name, b.score AS b_score, b.strand AS b_strand
            FROM intervals_a a
            JOIN intervals_b b ON a.interval INTERSECTS b.interval
            """,
            tables=["intervals_a", "intervals_b"],
        )
        giql_result = conn.execute(sql).fetchall()
    finally:
        conn.close()

    bedtools_result = intersect(
        [i.to_tuple() for i in intervals_a],
        [i.to_tuple() for i in intervals_b],
        write_both=True,
    )

    comparison = compare_results(giql_result, bedtools_result)
    assert comparison.match, comparison.failure_message()


# ---------------------------------------------------------------------------
# -c (count overlaps)
# ---------------------------------------------------------------------------


@given(
    intervals_a=unique_interval_list_st(max_size=20),
    intervals_b=unique_interval_list_st(max_size=20),
)
@settings(
    max_examples=40,
    deadline=None,
    suppress_health_check=[HealthCheck.too_slow],
)
def test_count_matches_bedtools_c(intervals_a, intervals_b):
    """
    GIVEN two randomly generated sets of unique genomic intervals
    WHEN GIQL COUNT of overlapping B per A is computed
    THEN results match bedtools intersect -c exactly
    """
    conn = duckdb.connect(":memory:")
    try:
        load_intervals(conn, "intervals_a", [i.to_tuple() for i in intervals_a])
        load_intervals(conn, "intervals_b", [i.to_tuple() for i in intervals_b])

        # Use a naive overlap join for counting — the binned join's
        # DISTINCT would collapse duplicate B matches.
        count_sql = """
            SELECT
                a.chrom, a."start", a."end", a.name, a.score, a.strand,
                COUNT(b.chrom) AS cnt
            FROM intervals_a a
            LEFT JOIN intervals_b b
                ON a.chrom = b.chrom
               AND a."start" < b."end"
               AND a."end" > b."start"
            GROUP BY a.chrom, a."start", a."end", a.name, a.score, a.strand
        """
        giql_result = conn.execute(count_sql).fetchall()
    finally:
        conn.close()

    bedtools_result = intersect(
        [i.to_tuple() for i in intervals_a],
        [i.to_tuple() for i in intervals_b],
        count=True,
    )

    comparison = compare_results(giql_result, bedtools_result)
    assert comparison.match, comparison.failure_message()


# ---------------------------------------------------------------------------
# -s (same strand)
# ---------------------------------------------------------------------------


@given(
    intervals_a=unique_interval_list_st(max_size=30),
    intervals_b=unique_interval_list_st(max_size=30),
)
@settings(
    max_examples=40,
    deadline=None,
    suppress_health_check=[HealthCheck.too_slow],
)
def test_same_strand_matches_bedtools_s(intervals_a, intervals_b):
    """
    GIVEN two randomly generated sets of genomic intervals with strands
    WHEN GIQL INTERSECTS with same-strand filter is executed
    THEN results match bedtools intersect -s exactly
    """
    conn = duckdb.connect(":memory:")
    try:
        load_intervals(conn, "intervals_a", [i.to_tuple() for i in intervals_a])
        load_intervals(conn, "intervals_b", [i.to_tuple() for i in intervals_b])

        sql = transpile(
            """
            SELECT DISTINCT a.*
            FROM intervals_a a, intervals_b b
            WHERE a.interval INTERSECTS b.interval
              AND a.strand = b.strand
            """,
            tables=["intervals_a", "intervals_b"],
        )
        giql_result = conn.execute(sql).fetchall()
    finally:
        conn.close()

    bedtools_result = intersect(
        [i.to_tuple() for i in intervals_a],
        [i.to_tuple() for i in intervals_b],
        strand_mode="same",
    )

    comparison = compare_results(giql_result, bedtools_result)
    assert comparison.match, comparison.failure_message()


# ---------------------------------------------------------------------------
# -S (opposite strand)
# ---------------------------------------------------------------------------


@given(
    intervals_a=unique_interval_list_st(max_size=30),
    intervals_b=unique_interval_list_st(max_size=30),
)
@settings(
    max_examples=40,
    deadline=None,
    suppress_health_check=[HealthCheck.too_slow],
)
def test_opposite_strand_matches_bedtools_S(intervals_a, intervals_b):
    """
    GIVEN two randomly generated sets of genomic intervals with strands
    WHEN GIQL INTERSECTS with opposite-strand filter is executed
    THEN results match bedtools intersect -S exactly
    """
    conn = duckdb.connect(":memory:")
    try:
        load_intervals(conn, "intervals_a", [i.to_tuple() for i in intervals_a])
        load_intervals(conn, "intervals_b", [i.to_tuple() for i in intervals_b])

        sql = transpile(
            """
            SELECT DISTINCT a.*
            FROM intervals_a a, intervals_b b
            WHERE a.interval INTERSECTS b.interval
              AND a.strand != b.strand
            """,
            tables=["intervals_a", "intervals_b"],
        )
        giql_result = conn.execute(sql).fetchall()
    finally:
        conn.close()

    bedtools_result = intersect(
        [i.to_tuple() for i in intervals_a],
        [i.to_tuple() for i in intervals_b],
        strand_mode="opposite",
    )

    comparison = compare_results(giql_result, bedtools_result)
    assert comparison.match, comparison.failure_message()


# ---------------------------------------------------------------------------
# -f (minimum overlap fraction of A)
# ---------------------------------------------------------------------------


@given(
    intervals_a=unique_interval_list_st(max_size=20),
    intervals_b=unique_interval_list_st(max_size=20),
    fraction=st.sampled_from([0.1, 0.25, 0.5, 0.75, 0.9]),
)
@settings(
    max_examples=40,
    deadline=None,
    suppress_health_check=[HealthCheck.too_slow],
)
def test_fraction_a_matches_bedtools_f(intervals_a, intervals_b, fraction):
    """
    GIVEN two randomly generated sets of genomic intervals and a fraction
    WHEN GIQL INTERSECTS with minimum overlap fraction of A is executed
    THEN results match bedtools intersect -f exactly
    """
    conn = duckdb.connect(":memory:")
    try:
        load_intervals(conn, "intervals_a", [i.to_tuple() for i in intervals_a])
        load_intervals(conn, "intervals_b", [i.to_tuple() for i in intervals_b])

        inner_sql = transpile(
            """
            SELECT DISTINCT
                a.chrom, a.start, a.end, a.name, a.score, a.strand,
                b.start AS b_start, b.end AS b_end
            FROM intervals_a a
            JOIN intervals_b b ON a.interval INTERSECTS b.interval
            """,
            tables=["intervals_a", "intervals_b"],
        )
        sql = f"""
            SELECT DISTINCT chrom, "start", "end", name, score, strand
            FROM ({inner_sql})
            WHERE (LEAST("end", b_end) - GREATEST("start", b_start))
                  >= {fraction} * ("end" - "start")
        """
        giql_result = conn.execute(sql).fetchall()
    finally:
        conn.close()

    bedtools_result = intersect(
        [i.to_tuple() for i in intervals_a],
        [i.to_tuple() for i in intervals_b],
        fraction_a=fraction,
    )

    comparison = compare_results(giql_result, bedtools_result)
    assert comparison.match, f"fraction_a={fraction}: {comparison.failure_message()}"


# ---------------------------------------------------------------------------
# -F (minimum overlap fraction of B)
# ---------------------------------------------------------------------------


@given(
    intervals_a=unique_interval_list_st(max_size=20),
    intervals_b=unique_interval_list_st(max_size=20),
    fraction=st.sampled_from([0.1, 0.25, 0.5, 0.75, 0.9]),
)
@settings(
    max_examples=40,
    deadline=None,
    suppress_health_check=[HealthCheck.too_slow],
)
def test_fraction_b_matches_bedtools_F(intervals_a, intervals_b, fraction):
    """
    GIVEN two randomly generated sets of genomic intervals and a fraction
    WHEN GIQL INTERSECTS with minimum overlap fraction of B is executed
    THEN results match bedtools intersect -F exactly
    """
    conn = duckdb.connect(":memory:")
    try:
        load_intervals(conn, "intervals_a", [i.to_tuple() for i in intervals_a])
        load_intervals(conn, "intervals_b", [i.to_tuple() for i in intervals_b])

        inner_sql = transpile(
            """
            SELECT DISTINCT
                a.chrom, a.start, a.end, a.name, a.score, a.strand,
                b.start AS b_start, b.end AS b_end
            FROM intervals_a a
            JOIN intervals_b b ON a.interval INTERSECTS b.interval
            """,
            tables=["intervals_a", "intervals_b"],
        )
        sql = f"""
            SELECT DISTINCT chrom, "start", "end", name, score, strand
            FROM ({inner_sql})
            WHERE (LEAST("end", b_end) - GREATEST("start", b_start))
                  >= {fraction} * (b_end - b_start)
        """
        giql_result = conn.execute(sql).fetchall()
    finally:
        conn.close()

    bedtools_result = intersect(
        [i.to_tuple() for i in intervals_a],
        [i.to_tuple() for i in intervals_b],
        fraction_b=fraction,
    )

    comparison = compare_results(giql_result, bedtools_result)
    assert comparison.match, f"fraction_b={fraction}: {comparison.failure_message()}"


# ---------------------------------------------------------------------------
# -r (reciprocal overlap fraction)
# ---------------------------------------------------------------------------


@given(
    intervals_a=unique_interval_list_st(max_size=20),
    intervals_b=unique_interval_list_st(max_size=20),
    fraction=st.sampled_from([0.1, 0.25, 0.5, 0.75]),
)
@settings(
    max_examples=40,
    deadline=None,
    suppress_health_check=[HealthCheck.too_slow],
)
def test_reciprocal_fraction_matches_bedtools_r(intervals_a, intervals_b, fraction):
    """
    GIVEN two randomly generated sets of genomic intervals and a fraction
    WHEN GIQL INTERSECTS with reciprocal overlap fraction is executed
    THEN results match bedtools intersect -f -F -r exactly
    """
    conn = duckdb.connect(":memory:")
    try:
        load_intervals(conn, "intervals_a", [i.to_tuple() for i in intervals_a])
        load_intervals(conn, "intervals_b", [i.to_tuple() for i in intervals_b])

        inner_sql = transpile(
            """
            SELECT DISTINCT
                a.chrom, a.start, a.end, a.name, a.score, a.strand,
                b.start AS b_start, b.end AS b_end
            FROM intervals_a a
            JOIN intervals_b b ON a.interval INTERSECTS b.interval
            """,
            tables=["intervals_a", "intervals_b"],
        )
        sql = f"""
            SELECT DISTINCT chrom, "start", "end", name, score, strand
            FROM ({inner_sql})
            WHERE (LEAST("end", b_end) - GREATEST("start", b_start))
                  >= {fraction} * ("end" - "start")
              AND (LEAST("end", b_end) - GREATEST("start", b_start))
                  >= {fraction} * (b_end - b_start)
        """
        giql_result = conn.execute(sql).fetchall()
    finally:
        conn.close()

    # -r applies -f reciprocally to both sides and requires -wa output.
    # Deduplicate to match GIQL's SELECT DISTINCT.
    bedtools_raw = intersect(
        [i.to_tuple() for i in intervals_a],
        [i.to_tuple() for i in intervals_b],
        fraction_a=fraction,
        reciprocal=True,
    )
    bedtools_result = list(set(bedtools_raw))

    comparison = compare_results(giql_result, bedtools_result)
    assert comparison.match, (
        f"reciprocal fraction={fraction}: {comparison.failure_message()}"
    )
