"""DataFusion execution smoke tests for the ``dialect="datafusion"`` target.

Issue #132 makes ``"datafusion"`` a routing key in :func:`giql.transpile`. The
existing ``TestBinnedJoinDataFusion`` suite proves the binned plan executes on
DataFusion, but always with ``dialect`` omitted (``None``); nothing executes the
``dialect="datafusion"`` entry point end-to-end. These tests close that seam.
"""

import pytest

from giql import transpile

pytestmark = pytest.mark.integration


class TestDataFusionTargetExecution:
    """End-to-end execution of dialect="datafusion" output on DataFusion."""

    def test_intersects_literal_returns_overlapping_row(self, datafusion_ctx):
        """Test a literal INTERSECTS query through the datafusion target.

        Given:
            A peaks table in a real DataFusion context with one overlapping
            interval, one non-overlapping interval, and one on another
            chromosome.
        When:
            A literal-range INTERSECTS query is transpiled with
            dialect="datafusion" and executed.
        Then:
            It should parse, execute, and return exactly the overlapping
            chr1 row.
        """
        # Arrange
        ctx = datafusion_ctx(
            peaks=[
                ("chr1", 1500, 1800),  # overlaps chr1:1000-2000
                ("chr1", 5000, 6000),  # same chrom, no overlap
                ("chr2", 1500, 1800),  # overlapping position, wrong chrom
            ]
        )
        sql = transpile(
            'SELECT chrom, start, "end" FROM peaks '
            "WHERE interval INTERSECTS 'chr1:1000-2000'",
            tables=["peaks"],
            dialect="datafusion",
        )

        # Act
        df = ctx.sql(sql).to_pandas()

        # Assert
        assert len(df) == 1
        assert df.iloc[0]["chrom"] == "chr1"
        assert df.iloc[0]["start"] == 1500

    def test_intersects_literal_no_overlap_returns_zero_rows(self, datafusion_ctx):
        """Test that a non-overlapping literal INTERSECTS yields no rows.

        Given:
            A peaks table in a real DataFusion context with no interval
            overlapping the query range.
        When:
            The literal INTERSECTS query is transpiled with
            dialect="datafusion" and executed.
        Then:
            It should execute and return zero rows.
        """
        # Arrange
        ctx = datafusion_ctx(
            peaks=[
                ("chr1", 5000, 6000),
                ("chr2", 1500, 1800),
            ]
        )
        sql = transpile(
            'SELECT chrom, start, "end" FROM peaks '
            "WHERE interval INTERSECTS 'chr1:1000-2000'",
            tables=["peaks"],
            dialect="datafusion",
        )

        # Act
        df = ctx.sql(sql).to_pandas()

        # Assert
        assert len(df) == 0

    def test_intersects_join_with_bin_size_dedups_across_bins(self, datafusion_ctx):
        """Test the datafusion-routed binned join honours an explicit bin size.

        Given:
            Two interval tables in a real DataFusion context whose single
            intervals overlap and span several bins under a small bin size.
        When:
            A column-to-column INTERSECTS join is transpiled with
            dialect="datafusion" and an explicit intersects_bin_size, then
            executed.
        Then:
            It should be accepted (not rejected as it would be for duckdb),
            the explicit bin size should reach emission (changing the SQL vs
            the default), and execution should return exactly one row — the
            overlap, de-duplicated by the pairs CTE across the shared bins.
        """
        # Arrange
        ctx = datafusion_ctx(
            peaks=[("chr1", 0, 5000)],
            genes=[("chr1", 500, 4500)],
        )
        query = (
            'SELECT a.chrom, a.start, a."end", '
            'b.start AS b_start, b."end" AS b_end '
            "FROM peaks a JOIN genes b ON a.interval INTERSECTS b.interval"
        )
        default_sql = transpile(query, tables=["peaks", "genes"], dialect="datafusion")
        sized_sql = transpile(
            query,
            tables=["peaks", "genes"],
            dialect="datafusion",
            intersects_bin_size=1000,
        )

        # Act
        df = ctx.sql(sized_sql).to_pandas()

        # Assert
        assert sized_sql != default_sql  # the explicit bin size reaches emission
        assert len(df) == 1  # de-duplicated across the shared bins
        assert df.iloc[0]["start"] == 0
        assert df.iloc[0]["b_start"] == 500
