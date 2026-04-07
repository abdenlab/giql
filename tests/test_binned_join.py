"""Tests for the INTERSECTS binned equi-join transpilation."""

import math

import pytest

from giql import Table
from giql import transpile


def _is_null(value) -> bool:
    """Check if a value is null/NaN (DataFusion returns NaN for nullable int64)."""
    if value is None:
        return True
    try:
        return math.isnan(value)
    except (TypeError, ValueError):
        return False


class TestTranspileBinnedJoin:
    """Unit tests for binned join SQL structure."""

    def test_basic_binned_join_rewrite(self):
        """
        GIVEN a GIQL query joining two tables with column-to-column INTERSECTS
        WHEN transpiling with default settings
        THEN should produce CTEs with UNNEST/range, equi-join and overlap in ON,
             and DISTINCT
        """
        sql = transpile(
            """
            SELECT a.*, b.*
            FROM peaks a
            JOIN genes b ON a.interval INTERSECTS b.interval
            """,
            tables=["peaks", "genes"],
        )

        sql_upper = sql.upper()

        # CTEs with UNNEST and range
        assert "WITH" in sql_upper
        assert "UNNEST" in sql_upper
        assert "range" in sql or "RANGE" in sql_upper
        assert "__giql_bin" in sql

        # Equi-join on chrom and bin
        assert '"chrom"' in sql
        assert "__giql_bin" in sql

        # Overlap filter in ON (not WHERE) for correct outer-join semantics
        assert "ON" in sql_upper
        assert '"start"' in sql or '"START"' in sql_upper
        assert '"end"' in sql or '"END"' in sql_upper

        # DISTINCT to deduplicate across bins
        assert "DISTINCT" in sql_upper

    def test_custom_bin_size(self):
        """
        GIVEN a GIQL query with column-to-column INTERSECTS join
        WHEN transpiling with intersects_bin_size=100000
        THEN should use 100000 in the range expressions
        """
        sql = transpile(
            """
            SELECT a.*, b.*
            FROM peaks a
            JOIN genes b ON a.interval INTERSECTS b.interval
            """,
            tables=["peaks", "genes"],
            intersects_bin_size=100000,
        )

        assert "100000" in sql

    def test_custom_column_mappings(self):
        """
        GIVEN two tables with different custom column schemas
        WHEN transpiling a binned join query
        THEN should use each table's custom column names in CTEs, ON, and WHERE
        """
        sql = transpile(
            """
            SELECT a.*, b.*
            FROM peaks a
            JOIN features b ON a.interval INTERSECTS b.location
            """,
            tables=[
                Table(
                    "peaks",
                    genomic_col="interval",
                    chrom_col="chromosome",
                    start_col="start_pos",
                    end_col="end_pos",
                ),
                Table(
                    "features",
                    genomic_col="location",
                    chrom_col="seqname",
                    start_col="begin",
                    end_col="terminus",
                ),
            ],
        )

        # Custom column names for peaks
        assert '"chromosome"' in sql
        assert '"start_pos"' in sql
        assert '"end_pos"' in sql

        # Custom column names for features
        assert '"seqname"' in sql
        assert '"begin"' in sql
        assert '"terminus"' in sql

        # Default column names should NOT appear
        assert '"chrom"' not in sql
        assert '"start"' not in sql
        assert '"end"' not in sql

    def test_literal_intersects_no_binned_ctes(self):
        """
        GIVEN a GIQL query with a literal-range INTERSECTS in WHERE (not a join)
        WHEN transpiling
        THEN should NOT produce binned CTEs
        """
        sql = transpile(
            "SELECT * FROM peaks WHERE interval INTERSECTS 'chr1:1000-2000'",
            tables=["peaks"],
        )

        assert "__giql_bin" not in sql
        assert "UNNEST" not in sql.upper()

    def test_no_join_passthrough(self):
        """
        GIVEN a simple SELECT query with no JOIN
        WHEN transpiling
        THEN should NOT produce binned CTEs
        """
        sql = transpile(
            "SELECT * FROM peaks",
            tables=["peaks"],
        )

        assert "__giql_bin" not in sql
        assert "UNNEST" not in sql.upper()
        assert "__giql_" not in sql

    def test_existing_where_preserved(self):
        """
        GIVEN a GIQL join query that already has a WHERE clause
        WHEN transpiling a binned join
        THEN should preserve the original WHERE condition alongside the overlap filter
        """
        sql = transpile(
            """
            SELECT a.*, b.*
            FROM peaks a
            JOIN genes b ON a.interval INTERSECTS b.interval
            WHERE a.score > 100
            """,
            tables=["peaks", "genes"],
        )

        sql_upper = sql.upper()

        # Original WHERE condition preserved
        assert "100" in sql
        assert "score" in sql.lower()

        # Overlap filter also present
        assert "WHERE" in sql_upper
        # Both conditions combined with AND
        assert "AND" in sql_upper

    def test_bin_size_none_defaults_to_10000(self):
        """
        GIVEN a GIQL join query
        WHEN transpiling with intersects_bin_size=None (explicit)
        THEN should produce the same output as default (10000)
        """
        sql_default = transpile(
            """
            SELECT a.*, b.*
            FROM peaks a
            JOIN genes b ON a.interval INTERSECTS b.interval
            """,
            tables=["peaks", "genes"],
        )

        sql_none = transpile(
            """
            SELECT a.*, b.*
            FROM peaks a
            JOIN genes b ON a.interval INTERSECTS b.interval
            """,
            tables=["peaks", "genes"],
            intersects_bin_size=None,
        )

        assert sql_default == sql_none
        assert "10000" in sql_default

    def test_implicit_cross_join_uses_binned_optimization(self):
        """
        GIVEN a GIQL query with implicit cross-join (FROM a, b WHERE INTERSECTS)
        WHEN transpiling
        THEN should use the binned equi-join optimization without leaking
             __giql_bin into SELECT * output columns
        """
        sql = transpile(
            """
            SELECT DISTINCT a.*
            FROM peaks a, genes b
            WHERE a.interval INTERSECTS b.interval
            """,
            tables=["peaks", "genes"],
        )

        # Binned CTEs are present
        assert "WITH" in sql.upper()
        assert "__giql_bin" in sql
        assert "UNNEST" in sql.upper()

        # Original table references preserved — no CTE leak into SELECT *
        assert "peaks" in sql
        assert '"chrom"' in sql
        assert '"start"' in sql
        assert '"end"' in sql

    def test_self_join_single_shared_cte(self):
        """
        GIVEN a self-join query where the same table appears with two aliases
        WHEN transpiling a binned join
        THEN should produce one shared key-only CTE for the underlying table,
             joined twice through distinct connector aliases
        """
        sql = transpile(
            """
            SELECT a.*, b.*
            FROM peaks a
            JOIN peaks b ON a.interval INTERSECTS b.interval
            """,
            tables=["peaks"],
        )

        sql_upper = sql.upper()

        # One shared CTE keyed on the table name
        assert "__giql_peaks_bins" in sql

        # Original table preserved in FROM
        assert "peaks" in sql

        # Should still have DISTINCT
        assert "DISTINCT" in sql_upper

    def test_invalid_bin_size_raises(self):
        """
        GIVEN intersects_bin_size=0 or a negative value
        WHEN calling transpile
        THEN should raise ValueError
        """
        with pytest.raises(ValueError, match="positive"):
            transpile(
                "SELECT * FROM a JOIN b ON a.interval INTERSECTS b.interval",
                tables=["a", "b"],
                intersects_bin_size=0,
            )

        with pytest.raises(ValueError, match="positive"):
            transpile(
                "SELECT * FROM a JOIN b ON a.interval INTERSECTS b.interval",
                tables=["a", "b"],
                intersects_bin_size=-1,
            )

    def test_multi_join_all_intersects_rewritten(self):
        """
        GIVEN a three-way join with two INTERSECTS conditions
        WHEN transpiling
        THEN should create one key-only CTE per underlying table and rewrite
             each INTERSECTS join as a three-join bridge through those CTEs
        """
        sql = transpile(
            """
            SELECT a.*, b.*, c.*
            FROM peaks a
            JOIN genes b ON a.interval INTERSECTS b.interval
            JOIN exons c ON a.interval INTERSECTS c.interval
            """,
            tables=["peaks", "genes", "exons"],
        )

        # One CTE per underlying table
        assert "__giql_peaks_bins" in sql
        assert "__giql_genes_bins" in sql
        assert "__giql_exons_bins" in sql

        # __giql_bin appears in CTE definitions and ON conditions
        sql_upper = sql.upper()
        assert sql_upper.count("__GIQL_BIN") >= 4  # at least 2 per INTERSECTS join

    def test_explicit_columns_uses_pairs_cte(self):
        """
        GIVEN a join query with only explicit columns in SELECT
        WHEN transpiling
        THEN should use pairs-CTE approach with key-only bin CTEs
        """
        sql = transpile(
            """
            SELECT a.chrom, a.start, b.start AS b_start
            FROM peaks a
            JOIN genes b ON a.interval INTERSECTS b.interval
            """,
            tables=["peaks", "genes"],
        )

        assert "__giql_peaks_bins" in sql
        assert "__giql_genes_bins" in sql
        assert "__giql_pairs_0" in sql

    def test_wildcard_select_uses_pairs_cte(self):
        """
        GIVEN a join query with wildcard expressions in SELECT
        WHEN transpiling
        THEN should use pairs-CTE approach, no __giql_bin in output
        """
        sql = transpile(
            """
            SELECT a.*, b.*
            FROM peaks a
            JOIN genes b ON a.interval INTERSECTS b.interval
            """,
            tables=["peaks", "genes"],
        )

        assert "__giql_peaks_bins" in sql
        assert "__giql_genes_bins" in sql
        assert "__giql_pairs_0" in sql


class TestBinnedJoinDataFusion:
    """End-to-end DataFusion correctness tests for binned INTERSECTS joins."""

    @staticmethod
    def _make_ctx(peaks_data, genes_data):
        """Create a DataFusion context with two interval tables."""
        import pyarrow as pa
        from datafusion import SessionContext

        schema = pa.schema(
            [
                ("chrom", pa.utf8()),
                ("start", pa.int64()),
                ("end", pa.int64()),
            ]
        )

        ctx = SessionContext()

        peaks_arrays = {
            "chrom": [r[0] for r in peaks_data],
            "start": [r[1] for r in peaks_data],
            "end": [r[2] for r in peaks_data],
        }
        genes_arrays = {
            "chrom": [r[0] for r in genes_data],
            "start": [r[1] for r in genes_data],
            "end": [r[2] for r in genes_data],
        }

        ctx.register_record_batches(
            "peaks", [pa.table(peaks_arrays, schema=schema).to_batches()]
        )
        ctx.register_record_batches(
            "genes", [pa.table(genes_arrays, schema=schema).to_batches()]
        )
        return ctx

    def test_overlapping_intervals_correct_rows_no_duplicates(self):
        """
        GIVEN two tables with overlapping intervals
        WHEN executing a binned INTERSECTS join via DataFusion
        THEN should return the correct matching rows with no duplicates
        """
        ctx = self._make_ctx(
            peaks_data=[("chr1", 100, 500), ("chr1", 1000, 2000)],
            genes_data=[("chr1", 300, 600), ("chr1", 5000, 6000)],
        )

        sql = transpile(
            """
            SELECT a.chrom, a.start, a."end", b.start AS b_start, b."end" AS b_end
            FROM peaks a
            JOIN genes b ON a.interval INTERSECTS b.interval
            """,
            tables=[
                Table("peaks", chrom_col="chrom", start_col="start", end_col="end"),
                Table("genes", chrom_col="chrom", start_col="start", end_col="end"),
            ],
        )

        df = ctx.sql(sql).to_pandas()

        # Only the first peak (100-500) overlaps the first gene (300-600)
        assert len(df) == 1
        assert df.iloc[0]["start"] == 100
        assert df.iloc[0]["b_start"] == 300

    def test_non_overlapping_intervals_zero_rows(self):
        """
        GIVEN two tables with no overlapping intervals
        WHEN executing a binned INTERSECTS join via DataFusion
        THEN should return zero rows
        """
        ctx = self._make_ctx(
            peaks_data=[("chr1", 100, 200), ("chr1", 300, 400)],
            genes_data=[("chr1", 500, 600), ("chr1", 700, 800)],
        )

        sql = transpile(
            """
            SELECT a.chrom, a.start, a."end"
            FROM peaks a
            JOIN genes b ON a.interval INTERSECTS b.interval
            """,
            tables=[
                Table("peaks", chrom_col="chrom", start_col="start", end_col="end"),
                Table("genes", chrom_col="chrom", start_col="start", end_col="end"),
            ],
        )

        df = ctx.sql(sql).to_pandas()
        assert len(df) == 0

    def test_adjacent_intervals_zero_rows_half_open(self):
        """
        GIVEN two tables with adjacent (touching) intervals under half-open coordinates
        WHEN executing a binned INTERSECTS join via DataFusion
        THEN should return zero rows because [100, 200) and [200, 300) do not overlap
        """
        ctx = self._make_ctx(
            peaks_data=[("chr1", 100, 200)],
            genes_data=[("chr1", 200, 300)],
        )

        sql = transpile(
            """
            SELECT a.chrom, a.start, a."end"
            FROM peaks a
            JOIN genes b ON a.interval INTERSECTS b.interval
            """,
            tables=[
                Table("peaks", chrom_col="chrom", start_col="start", end_col="end"),
                Table("genes", chrom_col="chrom", start_col="start", end_col="end"),
            ],
        )

        df = ctx.sql(sql).to_pandas()
        assert len(df) == 0

    def test_different_chromosomes_only_same_chrom(self):
        """
        GIVEN two tables with intervals on different chromosomes that would overlap positionally
        WHEN executing a binned INTERSECTS join via DataFusion
        THEN should only return overlaps on the same chromosome
        """
        ctx = self._make_ctx(
            peaks_data=[("chr1", 100, 500), ("chr2", 100, 500)],
            genes_data=[("chr1", 200, 400), ("chr3", 200, 400)],
        )

        sql = transpile(
            """
            SELECT a.chrom, a.start, a."end",
                   b.chrom AS b_chrom, b.start AS b_start, b."end" AS b_end
            FROM peaks a
            JOIN genes b ON a.interval INTERSECTS b.interval
            """,
            tables=[
                Table("peaks", chrom_col="chrom", start_col="start", end_col="end"),
                Table("genes", chrom_col="chrom", start_col="start", end_col="end"),
            ],
        )

        df = ctx.sql(sql).to_pandas()

        assert len(df) == 1
        assert df.iloc[0]["chrom"] == "chr1"
        assert df.iloc[0]["b_chrom"] == "chr1"

    def test_intervals_spanning_multiple_bins_no_duplicates(self):
        """
        GIVEN intervals that span multiple bins
        WHEN executing a binned INTERSECTS join via DataFusion
        THEN overlapping pairs should be returned exactly once (DISTINCT dedup)
        """
        ctx = self._make_ctx(
            peaks_data=[("chr1", 0, 50000)],
            genes_data=[("chr1", 25000, 75000)],
        )

        sql = transpile(
            """
            SELECT a.chrom, a.start, a."end", b.start AS b_start, b."end" AS b_end
            FROM peaks a
            JOIN genes b ON a.interval INTERSECTS b.interval
            """,
            tables=[
                Table("peaks", chrom_col="chrom", start_col="start", end_col="end"),
                Table("genes", chrom_col="chrom", start_col="start", end_col="end"),
            ],
            intersects_bin_size=10000,
        )

        df = ctx.sql(sql).to_pandas()

        # Despite sharing multiple bins (2, 3, 4), should appear exactly once
        assert len(df) == 1
        assert df.iloc[0]["start"] == 0
        assert df.iloc[0]["end"] == 50000
        assert df.iloc[0]["b_start"] == 25000
        assert df.iloc[0]["b_end"] == 75000

    def test_equivalence_with_naive_cross_join(self):
        """
        GIVEN two tables with a mix of overlapping and non-overlapping intervals
        WHEN executing a binned INTERSECTS join via DataFusion
        THEN results should match a naive cross-join with overlap filter
        """
        ctx = self._make_ctx(
            peaks_data=[
                ("chr1", 0, 100),
                ("chr1", 150, 300),
                ("chr1", 500, 1000),
                ("chr2", 0, 500),
            ],
            genes_data=[
                ("chr1", 50, 200),
                ("chr1", 250, 600),
                ("chr1", 900, 1100),
                ("chr2", 400, 800),
            ],
        )

        naive_sql = """
            SELECT a.chrom AS a_chrom, a.start AS a_start, a."end" AS a_end,
                   b.chrom AS b_chrom, b.start AS b_start, b."end" AS b_end
            FROM peaks a, genes b
            WHERE a.chrom = b.chrom
              AND a.start < b."end"
              AND a."end" > b.start
            ORDER BY a.chrom, a.start, b.start
        """
        naive_df = ctx.sql(naive_sql).to_pandas()

        binned_sql = transpile(
            """
            SELECT a.chrom AS a_chrom, a.start AS a_start, a."end" AS a_end,
                   b.chrom AS b_chrom, b.start AS b_start, b."end" AS b_end
            FROM peaks a
            JOIN genes b ON a.interval INTERSECTS b.interval
            """,
            tables=[
                Table("peaks", chrom_col="chrom", start_col="start", end_col="end"),
                Table("genes", chrom_col="chrom", start_col="start", end_col="end"),
            ],
        )
        binned_df = (
            ctx.sql(binned_sql)
            .to_pandas()
            .sort_values(by=["a_chrom", "a_start", "b_start"])
            .reset_index(drop=True)
        )
        naive_df = naive_df.reset_index(drop=True)

        assert len(binned_df) == len(naive_df)
        assert binned_df.values.tolist() == naive_df.values.tolist()

    def test_implicit_cross_join_correct_rows_no_bin_leak(self):
        """
        GIVEN two tables with overlapping intervals queried via implicit cross-join syntax
        WHEN executing a binned INTERSECTS join via DataFusion
        THEN results should be correct and SELECT a.* should not include __giql_bin
        """
        ctx = self._make_ctx(
            peaks_data=[("chr1", 100, 500), ("chr1", 1000, 2000)],
            genes_data=[("chr1", 300, 600), ("chr1", 5000, 6000)],
        )

        sql = transpile(
            """
            SELECT DISTINCT a.*
            FROM peaks a, genes b
            WHERE a.interval INTERSECTS b.interval
            """,
            tables=[
                Table("peaks", chrom_col="chrom", start_col="start", end_col="end"),
                Table("genes", chrom_col="chrom", start_col="start", end_col="end"),
            ],
        )

        df = ctx.sql(sql).to_pandas()

        # Only the first peak (100-500) overlaps the first gene (300-600)
        assert len(df) == 1
        assert df.iloc[0]["start"] == 100

        # SELECT a.* must return exactly the original table columns — no __giql_bin
        assert list(df.columns) == ["chrom", "start", "end"]


class TestBinnedJoinOuterJoinSemantics:
    """Regression tests: outer join kinds must be preserved after rewrite.

    Bug: the bridge path only applied the join kind (LEFT, RIGHT, FULL) to
    join3, while join1 and join2 were always INNER — silently converting
    outer joins into inner joins.
    """

    @staticmethod
    def _make_ctx(peaks_data, genes_data):
        """Create a DataFusion context with peaks and genes tables."""
        import pyarrow as pa
        from datafusion import SessionContext

        schema = pa.schema(
            [
                ("chrom", pa.utf8()),
                ("start", pa.int64()),
                ("end", pa.int64()),
            ]
        )
        ctx = SessionContext()
        ctx.register_record_batches(
            "peaks",
            [
                pa.table(
                    {
                        "chrom": [r[0] for r in peaks_data],
                        "start": [r[1] for r in peaks_data],
                        "end": [r[2] for r in peaks_data],
                    },
                    schema=schema,
                ).to_batches()
            ],
        )
        ctx.register_record_batches(
            "genes",
            [
                pa.table(
                    {
                        "chrom": [r[0] for r in genes_data],
                        "start": [r[1] for r in genes_data],
                        "end": [r[2] for r in genes_data],
                    },
                    schema=schema,
                ).to_batches()
            ],
        )
        return ctx

    def test_left_join_preserves_unmatched_left_rows_full_cte(self):
        """
        GIVEN peaks with one matching and one non-matching interval
        WHEN a LEFT JOIN with INTERSECTS is transpiled (no wildcards, full-CTE path)
        THEN the SQL must contain LEFT keyword and execution must return all
             left rows including unmatched ones with NULLs on the right
        """
        sql = transpile(
            """
            SELECT a.chrom, a.start, a."end", b.start AS b_start
            FROM peaks a
            LEFT JOIN genes b ON a.interval INTERSECTS b.interval
            """,
            tables=[
                Table("peaks", chrom_col="chrom", start_col="start", end_col="end"),
                Table("genes", chrom_col="chrom", start_col="start", end_col="end"),
            ],
        )

        assert "LEFT" in sql.upper()

        ctx = self._make_ctx(
            peaks_data=[("chr1", 100, 500), ("chr1", 1000, 2000)],
            genes_data=[("chr1", 300, 600)],
        )
        df = ctx.sql(sql).to_pandas().sort_values("start").reset_index(drop=True)

        assert len(df) == 2
        assert df.iloc[0]["start"] == 100
        assert df.iloc[0]["b_start"] == 300
        assert df.iloc[1]["start"] == 1000
        assert _is_null(df.iloc[1]["b_start"])

    def test_left_join_preserves_unmatched_left_rows_bridge(self):
        """
        GIVEN peaks with one matching and one non-matching interval
        WHEN a LEFT JOIN with INTERSECTS is transpiled (wildcards, bridge path)
        THEN the SQL must contain LEFT keyword and execution must return all
             left rows including unmatched ones with NULLs on the right
        """
        sql = transpile(
            """
            SELECT a.*, b.start AS b_start
            FROM peaks a
            LEFT JOIN genes b ON a.interval INTERSECTS b.interval
            """,
            tables=[
                Table("peaks", chrom_col="chrom", start_col="start", end_col="end"),
                Table("genes", chrom_col="chrom", start_col="start", end_col="end"),
            ],
        )

        assert "LEFT" in sql.upper()

        ctx = self._make_ctx(
            peaks_data=[("chr1", 100, 500), ("chr1", 1000, 2000)],
            genes_data=[("chr1", 300, 600)],
        )
        df = ctx.sql(sql).to_pandas().sort_values("start").reset_index(drop=True)

        assert len(df) == 2
        assert df.iloc[0]["start"] == 100
        assert df.iloc[0]["b_start"] == 300
        assert df.iloc[1]["start"] == 1000
        assert _is_null(df.iloc[1]["b_start"])

    def test_right_join_preserves_unmatched_right_rows_full_cte(self):
        """
        GIVEN genes with one matching and one non-matching interval
        WHEN a RIGHT JOIN with INTERSECTS is transpiled (no wildcards, full-CTE path)
        THEN the SQL must contain RIGHT keyword and execution must return all
             right rows including unmatched ones with NULLs on the left
        """
        sql = transpile(
            """
            SELECT a.start AS a_start, b.chrom, b.start, b."end"
            FROM peaks a
            RIGHT JOIN genes b ON a.interval INTERSECTS b.interval
            """,
            tables=[
                Table("peaks", chrom_col="chrom", start_col="start", end_col="end"),
                Table("genes", chrom_col="chrom", start_col="start", end_col="end"),
            ],
        )

        assert "RIGHT" in sql.upper()

        ctx = self._make_ctx(
            peaks_data=[("chr1", 100, 500)],
            genes_data=[("chr1", 300, 600), ("chr1", 5000, 6000)],
        )
        df = ctx.sql(sql).to_pandas().sort_values("start").reset_index(drop=True)

        assert len(df) == 2
        matched = df[df["a_start"].notna()]
        unmatched = df[df["a_start"].isna()]
        assert len(matched) == 1
        assert matched.iloc[0]["start"] == 300
        assert len(unmatched) == 1
        assert unmatched.iloc[0]["start"] == 5000

    def test_right_join_preserves_unmatched_right_rows_bridge(self):
        """
        GIVEN genes with one matching and one non-matching interval
        WHEN a RIGHT JOIN with INTERSECTS is transpiled (wildcards, bridge path)
        THEN the SQL must contain RIGHT keyword and execution must return all
             right rows including unmatched ones with NULLs on the left
        """
        sql = transpile(
            """
            SELECT a.start AS a_start, b.*
            FROM peaks a
            RIGHT JOIN genes b ON a.interval INTERSECTS b.interval
            """,
            tables=[
                Table("peaks", chrom_col="chrom", start_col="start", end_col="end"),
                Table("genes", chrom_col="chrom", start_col="start", end_col="end"),
            ],
        )

        assert "RIGHT" in sql.upper()

        ctx = self._make_ctx(
            peaks_data=[("chr1", 100, 500)],
            genes_data=[("chr1", 300, 600), ("chr1", 5000, 6000)],
        )
        df = ctx.sql(sql).to_pandas().sort_values("start").reset_index(drop=True)

        assert len(df) == 2
        matched = df[df["a_start"].notna()]
        unmatched = df[df["a_start"].isna()]
        assert len(matched) == 1
        assert matched.iloc[0]["start"] == 300
        assert len(unmatched) == 1
        assert unmatched.iloc[0]["start"] == 5000

    def test_full_outer_join_preserves_both_unmatched_full_cte(self):
        """
        GIVEN peaks and genes each with one matching and one non-matching interval
        WHEN a FULL OUTER JOIN with INTERSECTS is transpiled (no wildcards, full-CTE)
        THEN the SQL must contain FULL keyword and execution must return three
             rows: one matched pair plus one unmatched from each side
        """
        sql = transpile(
            """
            SELECT a.start AS a_start, a."end" AS a_end,
                   b.start AS b_start, b."end" AS b_end
            FROM peaks a
            FULL OUTER JOIN genes b ON a.interval INTERSECTS b.interval
            """,
            tables=[
                Table("peaks", chrom_col="chrom", start_col="start", end_col="end"),
                Table("genes", chrom_col="chrom", start_col="start", end_col="end"),
            ],
        )

        assert "FULL" in sql.upper()

        ctx = self._make_ctx(
            peaks_data=[("chr1", 100, 500), ("chr1", 8000, 9000)],
            genes_data=[("chr1", 300, 600), ("chr1", 5000, 6000)],
        )
        df = ctx.sql(sql).to_pandas()

        assert len(df) == 3
        matched = df[df["a_start"].notna() & df["b_start"].notna()]
        left_only = df[df["a_start"].notna() & df["b_start"].isna()]
        right_only = df[df["a_start"].isna() & df["b_start"].notna()]
        assert len(matched) == 1
        assert len(left_only) == 1
        assert len(right_only) == 1

    def test_full_outer_join_preserves_both_unmatched_bridge(self):
        """
        GIVEN peaks and genes each with one matching and one non-matching interval
        WHEN a FULL OUTER JOIN with INTERSECTS is transpiled (wildcards, bridge path)
        THEN the SQL must contain FULL keyword and execution must return three
             rows: one matched pair plus one unmatched from each side
        """
        sql = transpile(
            """
            SELECT a.*, b.start AS b_start, b."end" AS b_end
            FROM peaks a
            FULL OUTER JOIN genes b ON a.interval INTERSECTS b.interval
            """,
            tables=[
                Table("peaks", chrom_col="chrom", start_col="start", end_col="end"),
                Table("genes", chrom_col="chrom", start_col="start", end_col="end"),
            ],
        )

        assert "FULL" in sql.upper()

        ctx = self._make_ctx(
            peaks_data=[("chr1", 100, 500), ("chr1", 8000, 9000)],
            genes_data=[("chr1", 300, 600), ("chr1", 5000, 6000)],
        )
        df = ctx.sql(sql).to_pandas()

        assert len(df) == 3
        matched = df[df["start"].notna() & df["b_start"].notna()]
        left_only = df[df["start"].notna() & df["b_start"].isna()]
        right_only = df[df["start"].isna() & df["b_start"].notna()]
        assert len(matched) == 1
        assert len(left_only) == 1
        assert len(right_only) == 1

    def test_left_join_all_unmatched_returns_all_left_rows(self):
        """
        GIVEN peaks where no intervals overlap any gene
        WHEN a LEFT JOIN with INTERSECTS is transpiled
        THEN all left rows must still appear with NULLs on the right
        """
        sql = transpile(
            """
            SELECT a.chrom, a.start, a."end", b.start AS b_start
            FROM peaks a
            LEFT JOIN genes b ON a.interval INTERSECTS b.interval
            """,
            tables=[
                Table("peaks", chrom_col="chrom", start_col="start", end_col="end"),
                Table("genes", chrom_col="chrom", start_col="start", end_col="end"),
            ],
        )

        ctx = self._make_ctx(
            peaks_data=[("chr1", 100, 200), ("chr1", 300, 400)],
            genes_data=[("chr1", 500, 600)],
        )
        df = ctx.sql(sql).to_pandas()

        assert len(df) == 2
        assert df["b_start"].isna().all()


class TestBinnedJoinAdditionalOnConditions:
    """Regression tests: non-INTERSECTS conditions in ON must be preserved.

    Bug: the rewrite replaces the entire ON clause with the binned equi-join
    and overlap predicate, silently dropping any additional user conditions
    like ``AND a.score > b.score``.
    """

    @staticmethod
    def _make_ctx_with_score():
        """Create a DataFusion context with peaks and genes tables that include a score column."""
        import pyarrow as pa
        from datafusion import SessionContext

        schema = pa.schema(
            [
                ("chrom", pa.utf8()),
                ("start", pa.int64()),
                ("end", pa.int64()),
                ("score", pa.int64()),
            ]
        )
        ctx = SessionContext()
        ctx.register_record_batches(
            "peaks",
            [
                pa.table(
                    {
                        "chrom": ["chr1", "chr1"],
                        "start": [100, 100],
                        "end": [500, 500],
                        "score": [10, 50],
                    },
                    schema=schema,
                ).to_batches()
            ],
        )
        ctx.register_record_batches(
            "genes",
            [
                pa.table(
                    {
                        "chrom": ["chr1", "chr1"],
                        "start": [200, 200],
                        "end": [600, 600],
                        "score": [30, 60],
                    },
                    schema=schema,
                ).to_batches()
            ],
        )
        return ctx

    def test_additional_on_condition_preserved_full_cte(self):
        """
        GIVEN two overlapping intervals where only one pair satisfies score filter
        WHEN INTERSECTS is combined with a.score > b.score in ON (no wildcards)
        THEN the additional condition must survive the rewrite and filter results
        """
        sql = transpile(
            """
            SELECT a.chrom, a.start, a."end", a.score AS a_score, b.score AS b_score
            FROM peaks a
            JOIN genes b ON a.interval INTERSECTS b.interval AND a.score > b.score
            """,
            tables=[
                Table("peaks", chrom_col="chrom", start_col="start", end_col="end"),
                Table("genes", chrom_col="chrom", start_col="start", end_col="end"),
            ],
        )

        assert "score" in sql.lower()

        ctx = self._make_ctx_with_score()
        df = ctx.sql(sql).to_pandas()

        assert len(df) == 1
        assert df.iloc[0]["a_score"] == 50

    def test_additional_on_condition_preserved_bridge(self):
        """
        GIVEN two overlapping intervals where only one pair satisfies score filter
        WHEN INTERSECTS is combined with a.score > b.score in ON (wildcards)
        THEN the additional condition must survive the rewrite and filter results
        """
        sql = transpile(
            """
            SELECT a.*, b.score AS b_score
            FROM peaks a
            JOIN genes b ON a.interval INTERSECTS b.interval AND a.score > b.score
            """,
            tables=[
                Table("peaks", chrom_col="chrom", start_col="start", end_col="end"),
                Table("genes", chrom_col="chrom", start_col="start", end_col="end"),
            ],
        )

        assert "score" in sql.lower()

        ctx = self._make_ctx_with_score()
        df = ctx.sql(sql).to_pandas()

        assert len(df) == 1
        assert df.iloc[0]["score"] == 50

    def test_additional_on_condition_with_left_join(self):
        """
        GIVEN overlapping intervals with an extra ON condition that filters all matches
        WHEN LEFT JOIN with INTERSECTS AND a.score > b.score is used
        THEN unmatched left rows appear with NULL right columns
        """
        sql = transpile(
            """
            SELECT a.chrom, a.start, a.score AS a_score, b.score AS b_score
            FROM peaks a
            LEFT JOIN genes b
                ON a.interval INTERSECTS b.interval AND a.score > b.score
            """,
            tables=[
                Table("peaks", chrom_col="chrom", start_col="start", end_col="end"),
                Table("genes", chrom_col="chrom", start_col="start", end_col="end"),
            ],
        )

        ctx = self._make_ctx_with_score()
        df = ctx.sql(sql).to_pandas().sort_values("a_score").reset_index(drop=True)

        assert len(df) == 2
        row_low = df[df["a_score"] == 10].iloc[0]
        row_high = df[df["a_score"] == 50].iloc[0]
        assert _is_null(row_low["b_score"])
        assert row_high["b_score"] == 30

    def test_multiple_additional_conditions_preserved(self):
        """
        GIVEN overlapping intervals with two extra ON conditions
        WHEN INTERSECTS is combined with a.score > 20 AND b.score < 40 in ON
        THEN both conditions must survive the rewrite
        """
        sql = transpile(
            """
            SELECT a.chrom, a.start, a.score AS a_score, b.score AS b_score
            FROM peaks a
            JOIN genes b
                ON a.interval INTERSECTS b.interval
                AND a.score > 20
                AND b.score < 40
            """,
            tables=[
                Table("peaks", chrom_col="chrom", start_col="start", end_col="end"),
                Table("genes", chrom_col="chrom", start_col="start", end_col="end"),
            ],
        )

        sql_lower = sql.lower()
        assert "score" in sql_lower
        assert "20" in sql
        assert "40" in sql

    def test_additional_on_condition_implicit_cross_join(self):
        """
        GIVEN overlapping intervals queried via implicit cross-join with extra WHERE
        WHEN INTERSECTS is in WHERE alongside a.score > b.score
        THEN the score condition must be preserved in the output SQL
        """
        sql = transpile(
            """
            SELECT a.chrom, a.start, a.score AS a_score, b.score AS b_score
            FROM peaks a, genes b
            WHERE a.interval INTERSECTS b.interval AND a.score > b.score
            """,
            tables=[
                Table("peaks", chrom_col="chrom", start_col="start", end_col="end"),
                Table("genes", chrom_col="chrom", start_col="start", end_col="end"),
            ],
        )

        assert "score" in sql.lower()

        ctx = self._make_ctx_with_score()
        df = ctx.sql(sql).to_pandas()

        assert len(df) == 1
        assert df.iloc[0]["a_score"] == 50


class TestBinnedJoinDistinctSemantics:
    """Tests that the pairs-CTE approach preserves standard SQL bag semantics.

    The pairs CTE deduplicates on key columns internally, so the output
    query does not need SELECT DISTINCT.  This preserves legitimately
    duplicated source rows.
    """

    @staticmethod
    def _make_ctx_with_duplicates():
        """Create a DataFusion context where peaks has duplicate rows."""
        import pyarrow as pa
        from datafusion import SessionContext

        schema = pa.schema(
            [
                ("chrom", pa.utf8()),
                ("start", pa.int64()),
                ("end", pa.int64()),
            ]
        )
        ctx = SessionContext()
        ctx.register_record_batches(
            "peaks",
            [
                pa.table(
                    {
                        "chrom": ["chr1", "chr1"],
                        "start": [100, 100],
                        "end": [500, 500],
                    },
                    schema=schema,
                ).to_batches()
            ],
        )
        ctx.register_record_batches(
            "genes",
            [
                pa.table(
                    {
                        "chrom": ["chr1"],
                        "start": [200],
                        "end": [600],
                    },
                    schema=schema,
                ).to_batches()
            ],
        )
        return ctx

    def test_duplicate_rows_preserved_full_cte(self):
        """
        GIVEN peaks with two identical rows that both overlap one gene
        WHEN an inner join with INTERSECTS is transpiled (no wildcards)
        THEN both rows should be returned, matching naive cross-join behavior
        """
        ctx = self._make_ctx_with_duplicates()

        naive_sql = """
            SELECT a.chrom, a.start, a."end", b.start AS b_start
            FROM peaks a, genes b
            WHERE a.chrom = b.chrom AND a.start < b."end" AND a."end" > b.start
        """
        naive_df = ctx.sql(naive_sql).to_pandas()
        assert len(naive_df) == 2

        binned_sql = transpile(
            """
            SELECT a.chrom, a.start, a."end", b.start AS b_start
            FROM peaks a
            JOIN genes b ON a.interval INTERSECTS b.interval
            """,
            tables=[
                Table("peaks", chrom_col="chrom", start_col="start", end_col="end"),
                Table("genes", chrom_col="chrom", start_col="start", end_col="end"),
            ],
        )

        binned_df = ctx.sql(binned_sql).to_pandas()
        assert len(binned_df) == len(naive_df)

    def test_duplicate_rows_preserved_bridge(self):
        """
        GIVEN peaks with two identical rows that both overlap one gene
        WHEN an inner join with INTERSECTS is transpiled (wildcards)
        THEN both rows should be returned, matching naive cross-join behavior
        """
        ctx = self._make_ctx_with_duplicates()

        naive_sql = """
            SELECT a.chrom, a.start, a."end", b.start AS b_start
            FROM peaks a, genes b
            WHERE a.chrom = b.chrom AND a.start < b."end" AND a."end" > b.start
        """
        naive_df = ctx.sql(naive_sql).to_pandas()
        assert len(naive_df) == 2

        binned_sql = transpile(
            """
            SELECT a.*, b.start AS b_start
            FROM peaks a
            JOIN genes b ON a.interval INTERSECTS b.interval
            """,
            tables=[
                Table("peaks", chrom_col="chrom", start_col="start", end_col="end"),
                Table("genes", chrom_col="chrom", start_col="start", end_col="end"),
            ],
        )

        binned_df = ctx.sql(binned_sql).to_pandas()
        assert len(binned_df) == len(naive_df)

    def test_non_duplicate_rows_unaffected(self):
        """
        GIVEN peaks with two distinct rows that both overlap one gene
        WHEN an inner join with INTERSECTS is transpiled
        THEN DISTINCT does not collapse them because they differ
        """
        import pyarrow as pa
        from datafusion import SessionContext

        schema = pa.schema(
            [
                ("chrom", pa.utf8()),
                ("start", pa.int64()),
                ("end", pa.int64()),
            ]
        )
        ctx = SessionContext()
        ctx.register_record_batches(
            "peaks",
            [
                pa.table(
                    {
                        "chrom": ["chr1", "chr1"],
                        "start": [100, 150],
                        "end": [500, 550],
                    },
                    schema=schema,
                ).to_batches()
            ],
        )
        ctx.register_record_batches(
            "genes",
            [
                pa.table(
                    {
                        "chrom": ["chr1"],
                        "start": [200],
                        "end": [600],
                    },
                    schema=schema,
                ).to_batches()
            ],
        )

        binned_sql = transpile(
            """
            SELECT a.chrom, a.start, a."end", b.start AS b_start
            FROM peaks a
            JOIN genes b ON a.interval INTERSECTS b.interval
            """,
            tables=[
                Table("peaks", chrom_col="chrom", start_col="start", end_col="end"),
                Table("genes", chrom_col="chrom", start_col="start", end_col="end"),
            ],
        )

        df = ctx.sql(binned_sql).to_pandas()
        assert len(df) == 2

    def test_user_distinct_already_present_still_works(self):
        """
        GIVEN a query that already has SELECT DISTINCT
        WHEN the binned join rewrite also adds DISTINCT
        THEN the query must still execute correctly (no double-DISTINCT error)
        """
        import pyarrow as pa
        from datafusion import SessionContext

        schema = pa.schema(
            [
                ("chrom", pa.utf8()),
                ("start", pa.int64()),
                ("end", pa.int64()),
            ]
        )
        ctx = SessionContext()
        ctx.register_record_batches(
            "peaks",
            [
                pa.table(
                    {"chrom": ["chr1"], "start": [100], "end": [500]},
                    schema=schema,
                ).to_batches()
            ],
        )
        ctx.register_record_batches(
            "genes",
            [
                pa.table(
                    {"chrom": ["chr1"], "start": [200], "end": [600]},
                    schema=schema,
                ).to_batches()
            ],
        )

        binned_sql = transpile(
            """
            SELECT DISTINCT a.chrom, a.start, b.start AS b_start
            FROM peaks a
            JOIN genes b ON a.interval INTERSECTS b.interval
            """,
            tables=[
                Table("peaks", chrom_col="chrom", start_col="start", end_col="end"),
                Table("genes", chrom_col="chrom", start_col="start", end_col="end"),
            ],
        )

        df = ctx.sql(binned_sql).to_pandas()
        assert len(df) == 1

    def test_count_requires_distinguishing_columns(self):
        """
        GIVEN one interval A overlapping three intervals B
        WHEN the SELECT includes columns that distinguish each B match
        THEN DISTINCT preserves all three matches and COUNT is correct
        """
        import duckdb

        conn = duckdb.connect(":memory:")
        conn.execute(
            "CREATE TABLE regions "
            '(chrom VARCHAR, "start" INTEGER, "end" INTEGER, '
            "name VARCHAR)"
        )
        conn.execute(
            "CREATE TABLE features "
            '(chrom VARCHAR, "start" INTEGER, "end" INTEGER, '
            "name VARCHAR)"
        )
        # One large region
        conn.execute("INSERT INTO regions VALUES ('chr1', 100, 500, 'r1')")
        # Three features overlapping it
        conn.execute(
            "INSERT INTO features VALUES "
            "('chr1', 150, 200, 'f1'), "
            "('chr1', 250, 300, 'f2'), "
            "('chr1', 350, 400, 'f3')"
        )

        # With distinguishing columns: count should be 3
        inner_sql = transpile(
            """
            SELECT r.chrom, r.start, r.end, r.name,
                   f.start AS f_start, f.end AS f_end
            FROM regions r
            JOIN features f ON r.interval INTERSECTS f.interval
            """,
            tables=["regions", "features"],
        )
        count_sql = f"""
            SELECT name, COUNT(*) AS cnt
            FROM ({inner_sql})
            GROUP BY chrom, "start", "end", name
        """
        result = conn.execute(count_sql).fetchall()
        assert result[0][1] == 3, (
            f"Expected count=3 with distinguishing columns, got {result[0][1]}"
        )

        # Without distinguishing columns: count is still correct
        # because the pairs-CTE approach does not add SELECT DISTINCT
        inner_sql_no_dist = transpile(
            """
            SELECT r.chrom, r.start, r.end, r.name, f.chrom AS f_chrom
            FROM regions r
            JOIN features f ON r.interval INTERSECTS f.interval
            """,
            tables=["regions", "features"],
        )
        count_sql_no_dist = f"""
            SELECT name, COUNT(*) AS cnt
            FROM ({inner_sql_no_dist})
            GROUP BY chrom, "start", "end", name
        """
        result = conn.execute(count_sql_no_dist).fetchall()
        assert result[0][1] == 3, (
            f"Expected count=3 even without distinguishing columns, got {result[0][1]}"
        )
        conn.close()


class TestBinnedJoinBinBoundaryRounding:
    """Regression tests for bin-index calculation rounding errors.

    The original formula CAST(start / B AS BIGINT) uses float division
    followed by a cast.  When the division lands on x.5 the cast rounds
    to nearest-even instead of flooring, producing the wrong bin index
    and causing missed matches.
    """

    @staticmethod
    def _make_ctx(table_a_data, table_b_data):
        import pyarrow as pa
        from datafusion import SessionContext

        schema = pa.schema(
            [
                ("chrom", pa.utf8()),
                ("start", pa.int64()),
                ("end", pa.int64()),
                ("name", pa.utf8()),
            ]
        )
        ctx = SessionContext()
        ctx.register_record_batches(
            "intervals_a",
            [pa.table(table_a_data, schema=schema).to_batches()],
        )
        ctx.register_record_batches(
            "intervals_b",
            [pa.table(table_b_data, schema=schema).to_batches()],
        )
        return ctx

    def test_half_bin_boundary_overlap_not_missed(self):
        """
        GIVEN interval A spanning many bins and interval B whose start
              falls exactly on a .5 division boundary (e.g., 621950/100)
        WHEN INTERSECTS is evaluated with intersects_bin_size=100 on DuckDB
        THEN the overlap must be found, not missed due to rounding
        """
        import duckdb

        conn = duckdb.connect(":memory:")
        conn.execute(
            "CREATE TABLE intervals_a "
            '(chrom VARCHAR, "start" INTEGER, "end" INTEGER, '
            "name VARCHAR)"
        )
        conn.execute(
            "CREATE TABLE intervals_b "
            '(chrom VARCHAR, "start" INTEGER, "end" INTEGER, '
            "name VARCHAR)"
        )
        conn.execute("INSERT INTO intervals_a VALUES ('chr1', 421951, 621951, 'a0')")
        conn.execute("INSERT INTO intervals_b VALUES ('chr1', 621950, 621951, 'b0')")

        sql = transpile(
            """
            SELECT DISTINCT a.name, b.name AS b_name
            FROM intervals_a a
            JOIN intervals_b b ON a.interval INTERSECTS b.interval
            """,
            tables=["intervals_a", "intervals_b"],
            intersects_bin_size=100,
        )
        result = conn.execute(sql).fetchall()
        conn.close()
        assert len(result) == 1, (
            f"Expected 1 match, got {len(result)} — "
            f"bin boundary rounding likely dropped the overlap"
        )

    def test_exact_bin_boundary_start(self):
        """
        GIVEN interval B starting at an exact multiple of bin_size
        WHEN INTERSECTS is evaluated on DuckDB
        THEN the correct bin index is assigned (no off-by-one from rounding)
        """
        import duckdb

        conn = duckdb.connect(":memory:")
        conn.execute(
            "CREATE TABLE intervals_a "
            '(chrom VARCHAR, "start" INTEGER, "end" INTEGER, '
            "name VARCHAR)"
        )
        conn.execute(
            "CREATE TABLE intervals_b "
            '(chrom VARCHAR, "start" INTEGER, "end" INTEGER, '
            "name VARCHAR)"
        )
        conn.execute("INSERT INTO intervals_a VALUES ('chr1', 999, 1001, 'a0')")
        conn.execute("INSERT INTO intervals_b VALUES ('chr1', 1000, 1001, 'b0')")

        sql = transpile(
            """
            SELECT DISTINCT a.name, b.name AS b_name
            FROM intervals_a a
            JOIN intervals_b b ON a.interval INTERSECTS b.interval
            """,
            tables=["intervals_a", "intervals_b"],
            intersects_bin_size=1000,
        )
        result = conn.execute(sql).fetchall()
        conn.close()
        assert len(result) == 1, f"Expected 1 match at bin boundary, got {len(result)}"


class TestBinnedJoinOuterJoinMultiBin:
    """Regression tests for outer join with multi-bin intervals.

    When an interval spans multiple bins, the outer join produces one
    row per bin.  Bins that don't match the other side create spurious
    NULL rows.  DISTINCT can't collapse a NULL row with a matched row
    because they differ in the non-NULL columns.
    """

    @staticmethod
    def _make_ctx(table_a_data, table_b_data):
        import pyarrow as pa
        from datafusion import SessionContext

        schema = pa.schema(
            [
                ("chrom", pa.utf8()),
                ("start", pa.int64()),
                ("end", pa.int64()),
                ("name", pa.utf8()),
            ]
        )
        ctx = SessionContext()
        ctx.register_record_batches(
            "intervals_a",
            [pa.table(table_a_data, schema=schema).to_batches()],
        )
        ctx.register_record_batches(
            "intervals_b",
            [pa.table(table_b_data, schema=schema).to_batches()],
        )
        return ctx

    def test_left_join_no_spurious_null_row(self):
        """
        GIVEN interval A spanning bins 0 and 1 and interval B only in bin 1
        WHEN LEFT JOIN INTERSECTS is evaluated
        THEN only 1 matched row is returned, not a matched row plus a
             spurious NULL row from the unmatched bin-0 copy
        """
        ctx = self._make_ctx(
            {
                "chrom": ["chr1"],
                "start": [9000],
                "end": [11000],
                "name": ["a0"],
            },
            {
                "chrom": ["chr1"],
                "start": [10500],
                "end": [10600],
                "name": ["b0"],
            },
        )

        sql = transpile(
            """
            SELECT a.name, b.name AS b_name
            FROM intervals_a a
            LEFT JOIN intervals_b b ON a.interval INTERSECTS b.interval
            """,
            tables=["intervals_a", "intervals_b"],
        )
        result = ctx.sql(sql).to_pandas()
        assert len(result) == 1, (
            f"Expected 1 matched row, got {len(result)} — "
            f"spurious NULL row from unmatched bin"
        )
        assert result.iloc[0]["b_name"] == "b0"

    def test_left_join_unmatched_still_returns_null(self):
        """
        GIVEN interval A with no overlap in B
        WHEN LEFT JOIN INTERSECTS is evaluated
        THEN one row with NULL B columns is returned
        """
        ctx = self._make_ctx(
            {
                "chrom": ["chr1"],
                "start": [9000],
                "end": [11000],
                "name": ["a0"],
            },
            {
                "chrom": ["chr2"],
                "start": [9500],
                "end": [10500],
                "name": ["b0"],
            },
        )

        sql = transpile(
            """
            SELECT a.name, b.name AS b_name
            FROM intervals_a a
            LEFT JOIN intervals_b b ON a.interval INTERSECTS b.interval
            """,
            tables=["intervals_a", "intervals_b"],
        )
        result = ctx.sql(sql).to_pandas()
        assert len(result) == 1, f"Expected 1 unmatched row, got {len(result)}"
        assert _is_null(result.iloc[0]["b_name"])

    def test_right_join_no_spurious_null_row(self):
        """
        GIVEN interval B spanning bins 0 and 1 and interval A only in bin 0
        WHEN RIGHT JOIN INTERSECTS is evaluated
        THEN only 1 matched row is returned, not a matched row plus a
             spurious NULL row from the unmatched bin-1 copy of B
        """
        ctx = self._make_ctx(
            {
                "chrom": ["chr1"],
                "start": [9500],
                "end": [9600],
                "name": ["a0"],
            },
            {
                "chrom": ["chr1"],
                "start": [9000],
                "end": [11000],
                "name": ["b0"],
            },
        )

        sql = transpile(
            """
            SELECT a.name, b.name AS b_name
            FROM intervals_a a
            RIGHT JOIN intervals_b b ON a.interval INTERSECTS b.interval
            """,
            tables=["intervals_a", "intervals_b"],
        )
        result = ctx.sql(sql).to_pandas()
        assert len(result) == 1, (
            f"Expected 1 matched row, got {len(result)} — "
            f"spurious NULL row from unmatched bin"
        )
        assert result.iloc[0]["name"] == "a0"

    def test_full_outer_join_no_spurious_null_row(self):
        """
        GIVEN interval A spanning bins 0 and 1, interval B only in bin 1
        WHEN FULL OUTER JOIN INTERSECTS is evaluated
        THEN only 1 matched row is returned, not a matched row plus a
             spurious NULL row
        """
        ctx = self._make_ctx(
            {
                "chrom": ["chr1"],
                "start": [9000],
                "end": [11000],
                "name": ["a0"],
            },
            {
                "chrom": ["chr1"],
                "start": [10500],
                "end": [10600],
                "name": ["b0"],
            },
        )

        sql = transpile(
            """
            SELECT a.name, b.name AS b_name
            FROM intervals_a a
            FULL OUTER JOIN intervals_b b ON a.interval INTERSECTS b.interval
            """,
            tables=["intervals_a", "intervals_b"],
        )
        result = ctx.sql(sql).to_pandas()
        assert len(result) == 1, (
            f"Expected 1 matched row, got {len(result)} — "
            f"spurious NULL row from unmatched bin"
        )
