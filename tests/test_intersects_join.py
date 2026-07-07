"""Tests for the INTERSECTS naive overlap predicate join transpilation."""

import math

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


class TestTranspileIntersectsJoin:
    """Unit tests for the naive overlap predicate SQL structure."""

    def test_inner_join_emits_naive_overlap_predicate(self):
        """
        GIVEN a GIQL query joining two tables with column-to-column INTERSECTS
        WHEN transpiling with default settings
        THEN the INTERSECTS expands in place to the naive overlap predicate
             (chrom equality plus half-open start/end comparisons) with no
             CTEs, no UNNEST, no bin columns, and no added DISTINCT
        """
        sql = transpile(
            """
            SELECT a.*, b.*
            FROM peaks a
            JOIN genes b ON a.interval INTERSECTS b.interval
            """,
            tables=["peaks", "genes"],
        )

        # Naive overlap predicate: chrom equality and half-open comparisons
        assert '"chrom" = ' in sql
        assert '"start" <' in sql
        assert '"end" >' in sql

        # No CTE rewrite machinery
        assert "WITH" not in sql.upper()
        assert "UNNEST" not in sql.upper()
        assert "__giql_bin" not in sql
        assert "__giql_pairs" not in sql

        # No DISTINCT added by the transpiler
        assert "DISTINCT" not in sql.upper()

    def test_custom_column_mappings(self):
        """
        GIVEN two tables with different custom column schemas
        WHEN transpiling an INTERSECTS join query
        THEN should use each table's custom column names in the overlap predicate
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

    def test_literal_intersects_no_ctes(self):
        """
        GIVEN a GIQL query with a literal-range INTERSECTS in WHERE (not a join)
        WHEN transpiling
        THEN should NOT produce any CTEs
        """
        sql = transpile(
            "SELECT * FROM peaks WHERE interval INTERSECTS 'chr1:1000-2000'",
            tables=["peaks"],
        )

        assert "WITH" not in sql.upper()
        assert "__giql_bin" not in sql
        assert "UNNEST" not in sql.upper()

    def test_no_join_passthrough(self):
        """
        GIVEN a simple SELECT query with no JOIN
        WHEN transpiling
        THEN should NOT produce any CTEs or helper columns
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
        WHEN transpiling an INTERSECTS join
        THEN should preserve the original WHERE condition alongside the overlap
             predicate in the ON clause
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

        # Overlap predicate also present
        assert "WHERE" in sql_upper
        assert "AND" in sql_upper

    def test_implicit_cross_join_emits_naive_predicate(self):
        """
        GIVEN a GIQL query with implicit cross-join (FROM a, b WHERE INTERSECTS)
        WHEN transpiling
        THEN the INTERSECTS expands in place to the naive overlap predicate in
             WHERE, with no CTEs and no bin columns leaking into SELECT *
        """
        sql = transpile(
            """
            SELECT DISTINCT a.*
            FROM peaks a, genes b
            WHERE a.interval INTERSECTS b.interval
            """,
            tables=["peaks", "genes"],
        )

        # Naive overlap predicate present
        assert '"chrom" = ' in sql
        assert '"start" <' in sql
        assert '"end" >' in sql

        # No CTE rewrite machinery
        assert "WITH" not in sql.upper()
        assert "UNNEST" not in sql.upper()
        assert "__giql_bin" not in sql
        assert "__giql_pairs" not in sql

        # Original table references preserved
        assert "peaks" in sql

    def test_self_join_emits_naive_predicate(self):
        """
        GIVEN a self-join query where the same table appears with two aliases
        WHEN transpiling an INTERSECTS join
        THEN the INTERSECTS expands in place to the naive overlap predicate with
             no CTEs, keeping the original table in FROM
        """
        sql = transpile(
            """
            SELECT a.*, b.*
            FROM peaks a
            JOIN peaks b ON a.interval INTERSECTS b.interval
            """,
            tables=["peaks"],
        )

        # Naive overlap predicate present
        assert '"chrom" = ' in sql
        assert '"start" <' in sql
        assert '"end" >' in sql

        # No CTE rewrite machinery
        assert "WITH" not in sql.upper()
        assert "UNNEST" not in sql.upper()
        assert "__giql_bin" not in sql

        # Original table preserved in FROM
        assert "peaks" in sql

        # No DISTINCT added by the transpiler
        assert "DISTINCT" not in sql.upper()

    def test_multi_join_expands_each_intersects(self):
        """
        GIVEN a three-way join with two INTERSECTS conditions
        WHEN transpiling
        THEN each INTERSECTS expands in place to its own naive overlap predicate,
             with no CTEs
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

        # Both INTERSECTS rendered as naive overlap predicates
        assert sql.count('"chrom" = ') == 2
        assert sql.count('"start" <') == 2
        assert sql.count('"end" >') == 2

        # No CTE rewrite machinery
        assert "WITH" not in sql.upper()
        assert "UNNEST" not in sql.upper()
        assert "__giql_bin" not in sql

    def test_explicit_columns_emit_naive_predicate(self):
        """
        GIVEN a join query with only explicit columns in SELECT
        WHEN transpiling
        THEN the INTERSECTS expands in place to the naive overlap predicate with
             no CTEs
        """
        sql = transpile(
            """
            SELECT a.chrom, a.start, b.start AS b_start
            FROM peaks a
            JOIN genes b ON a.interval INTERSECTS b.interval
            """,
            tables=["peaks", "genes"],
        )

        # Naive overlap predicate present
        assert '"chrom" = ' in sql
        assert '"start" <' in sql
        assert '"end" >' in sql

        # No CTE rewrite machinery
        assert "WITH" not in sql.upper()
        assert "UNNEST" not in sql.upper()
        assert "__giql_bin" not in sql
        assert "__giql_pairs" not in sql


class TestIntersectsJoinDataFusion:
    """End-to-end DataFusion correctness tests for INTERSECTS joins."""

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
        WHEN executing an INTERSECTS join via DataFusion
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
        WHEN executing an INTERSECTS join via DataFusion
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
        WHEN executing an INTERSECTS join via DataFusion
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
        GIVEN two tables with intervals on different chromosomes that would
              overlap positionally
        WHEN executing an INTERSECTS join via DataFusion
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

    def test_large_span_overlap_returned_once(self):
        """
        GIVEN two large intervals that overlap over a wide span
        WHEN executing an INTERSECTS join via DataFusion
        THEN the overlapping pair should be returned exactly once
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
        )

        df = ctx.sql(sql).to_pandas()

        # The single overlapping pair should appear exactly once
        assert len(df) == 1
        assert df.iloc[0]["start"] == 0
        assert df.iloc[0]["end"] == 50000
        assert df.iloc[0]["b_start"] == 25000
        assert df.iloc[0]["b_end"] == 75000

    def test_equivalence_with_naive_cross_join(self):
        """
        GIVEN two tables with a mix of overlapping and non-overlapping intervals
        WHEN executing an INTERSECTS join via DataFusion
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

        intersects_sql = transpile(
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
        intersects_df = (
            ctx.sql(intersects_sql)
            .to_pandas()
            .sort_values(by=["a_chrom", "a_start", "b_start"])
            .reset_index(drop=True)
        )
        naive_df = naive_df.reset_index(drop=True)

        assert len(intersects_df) == len(naive_df)
        assert intersects_df.values.tolist() == naive_df.values.tolist()

    def test_implicit_cross_join_correct_rows_no_column_leak(self):
        """
        GIVEN two tables with overlapping intervals queried via implicit
              cross-join syntax
        WHEN executing an INTERSECTS join via DataFusion
        THEN results should be correct and SELECT a.* should return only the
             original table columns
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

        # SELECT a.* must return exactly the original table columns
        assert list(df.columns) == ["chrom", "start", "end"]


class TestIntersectsJoinOuterJoinSemantics:
    """Outer join kinds must be preserved when INTERSECTS expands in place.

    The naive overlap predicate is emitted directly in the ON clause, so
    LEFT, RIGHT, and FULL OUTER joins keep their outer semantics and return
    unmatched rows with NULLs on the opposite side.
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

    def test_left_join_preserves_unmatched_left_rows_explicit_cols(self):
        """
        GIVEN peaks with one matching and one non-matching interval
        WHEN a LEFT JOIN with INTERSECTS is transpiled (explicit columns)
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

    def test_left_join_preserves_unmatched_left_rows_wildcard(self):
        """
        GIVEN peaks with one matching and one non-matching interval
        WHEN a LEFT JOIN with INTERSECTS is transpiled (wildcard SELECT)
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

    def test_right_join_preserves_unmatched_right_rows_explicit_cols(self):
        """
        GIVEN genes with one matching and one non-matching interval
        WHEN a RIGHT JOIN with INTERSECTS is transpiled (explicit columns)
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

    def test_right_join_preserves_unmatched_right_rows_wildcard(self):
        """
        GIVEN genes with one matching and one non-matching interval
        WHEN a RIGHT JOIN with INTERSECTS is transpiled (wildcard SELECT)
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

    def test_full_outer_join_preserves_both_unmatched_explicit_cols(self):
        """
        GIVEN peaks and genes each with one matching and one non-matching interval
        WHEN a FULL OUTER JOIN with INTERSECTS is transpiled (explicit columns)
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

    def test_full_outer_join_preserves_both_unmatched_wildcard(self):
        """
        GIVEN peaks and genes each with one matching and one non-matching interval
        WHEN a FULL OUTER JOIN with INTERSECTS is transpiled (wildcard SELECT)
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


class TestIntersectsJoinAdditionalOnConditions:
    """Non-INTERSECTS conditions in ON must be preserved.

    INTERSECTS expands in place within the ON clause, so any additional user
    conditions like ``AND a.score > b.score`` remain alongside the naive
    overlap predicate and continue to filter results.
    """

    @staticmethod
    def _make_ctx_with_score():
        """Create a DataFusion context with peaks and genes tables including a score."""
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

    def test_additional_on_condition_preserved_explicit_cols(self):
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

    def test_additional_on_condition_preserved_wildcard(self):
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


class TestIntersectsJoinBagSemantics:
    """The naive overlap predicate preserves standard SQL bag semantics.

    INTERSECTS expands in place without adding SELECT DISTINCT, so the join
    behaves like an ordinary range join and legitimately duplicated source
    rows are preserved in the output.
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

    def test_duplicate_rows_preserved_explicit_cols(self):
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

        intersects_sql = transpile(
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

        intersects_df = ctx.sql(intersects_sql).to_pandas()
        assert len(intersects_df) == len(naive_df)

    def test_duplicate_rows_preserved_wildcard(self):
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

        intersects_sql = transpile(
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

        intersects_df = ctx.sql(intersects_sql).to_pandas()
        assert len(intersects_df) == len(naive_df)

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

        intersects_sql = transpile(
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

        df = ctx.sql(intersects_sql).to_pandas()
        assert len(df) == 2

    def test_user_distinct_already_present_still_works(self):
        """
        GIVEN a query that already has SELECT DISTINCT
        WHEN the INTERSECTS join is transpiled
        THEN the user's DISTINCT is preserved and the query executes correctly
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

        intersects_sql = transpile(
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

        df = ctx.sql(intersects_sql).to_pandas()
        assert len(df) == 1

    def test_count_requires_distinguishing_columns(self):
        """
        GIVEN one interval A overlapping three intervals B
        WHEN the SELECT includes columns that distinguish each B match
        THEN all three matches are preserved and COUNT is correct
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
        # because the naive overlap predicate does not add SELECT DISTINCT
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


class TestIntersectsJoinPositionBoundaries:
    """Overlap correctness at exact and near round-number positions.

    These cases probe intervals whose start lands exactly on or adjacent to
    round-number coordinates, confirming the naive overlap predicate finds
    the overlap without any off-by-one error.
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

    def test_overlap_at_position_boundary_found(self):
        """
        GIVEN a large interval A and interval B whose start falls on a
              round-number position boundary (621950)
        WHEN INTERSECTS is evaluated on DuckDB
        THEN the overlap must be found
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
        )
        result = conn.execute(sql).fetchall()
        conn.close()
        assert len(result) == 1, f"Expected 1 match, got {len(result)}"

    def test_overlap_at_exact_multiple_found(self):
        """
        GIVEN interval B starting at an exact round-number position (1000)
        WHEN INTERSECTS is evaluated on DuckDB
        THEN the overlap is found with no off-by-one
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
        )
        result = conn.execute(sql).fetchall()
        conn.close()
        assert len(result) == 1, f"Expected 1 match, got {len(result)}"


class TestIntersectsJoinOuterJoinNoSpuriousRows:
    """Outer joins must not emit spurious NULL rows for matched intervals.

    Because INTERSECTS expands to a single naive overlap predicate in the ON
    clause, an interval that overlaps produces exactly one matched row and no
    extra NULL-filled row appears on the outer side.
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
        GIVEN a large interval A overlapping a single interval B
        WHEN LEFT JOIN INTERSECTS is evaluated
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
            LEFT JOIN intervals_b b ON a.interval INTERSECTS b.interval
            """,
            tables=["intervals_a", "intervals_b"],
        )
        result = ctx.sql(sql).to_pandas()
        assert len(result) == 1, (
            f"Expected 1 matched row, got {len(result)} — spurious NULL row"
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
        GIVEN a large interval B overlapping a single interval A
        WHEN RIGHT JOIN INTERSECTS is evaluated
        THEN only 1 matched row is returned, not a matched row plus a
             spurious NULL row
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
            f"Expected 1 matched row, got {len(result)} — spurious NULL row"
        )
        assert result.iloc[0]["name"] == "a0"

    def test_full_outer_join_no_spurious_null_row(self):
        """
        GIVEN a large interval A overlapping a single interval B
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
            f"Expected 1 matched row, got {len(result)} — spurious NULL row"
        )
