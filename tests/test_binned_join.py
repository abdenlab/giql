"""Tests for the INTERSECTS binned equi-join transpilation."""

from giql import Table
from giql import transpile


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
        WHEN transpiling with bin_size=100000
        THEN should use 100000 in the range expressions
        """
        sql = transpile(
            """
            SELECT a.*, b.*
            FROM peaks a
            JOIN genes b ON a.interval INTERSECTS b.interval
            """,
            tables=["peaks", "genes"],
            bin_size=100000,
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
        WHEN transpiling with bin_size=None (explicit)
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
            bin_size=None,
        )

        assert sql_default == sql_none
        assert "10000" in sql_default

    def test_implicit_cross_join_uses_naive_predicate(self):
        """
        GIVEN a GIQL query with implicit cross-join (FROM a, b WHERE INTERSECTS)
        WHEN transpiling
        THEN should use the naive overlap predicate, not binned CTEs, because
             the CTE approach leaks __giql_bin into SELECT * results
        """
        sql = transpile(
            """
            SELECT DISTINCT a.*
            FROM peaks a, genes b
            WHERE a.interval INTERSECTS b.interval
            """,
            tables=["peaks", "genes"],
        )

        # No binned CTEs — falls through to generator's naive predicate
        assert "__giql_bin" not in sql
        assert "UNNEST" not in sql.upper()

        # Naive overlap predicate present
        assert '"chrom"' in sql
        assert '"start"' in sql
        assert '"end"' in sql

    def test_self_join_distinct_ctes(self):
        """
        GIVEN a self-join query where the same table appears with two aliases
        WHEN transpiling a binned join
        THEN should produce two distinct CTEs both referencing the same underlying table
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

        # Two distinct CTEs
        assert "__giql_a_binned" in sql
        assert "__giql_b_binned" in sql

        # Both reference the same underlying table
        assert "FROM peaks" in sql or "FROM PEAKS" in sql_upper

        # Should still have DISTINCT
        assert "DISTINCT" in sql_upper

    def test_invalid_bin_size_raises(self):
        """
        GIVEN bin_size=0 or a negative value
        WHEN calling transpile
        THEN should raise ValueError
        """
        import pytest

        with pytest.raises(ValueError, match="positive"):
            transpile(
                "SELECT * FROM a JOIN b ON a.interval INTERSECTS b.interval",
                tables=["a", "b"],
                bin_size=0,
            )

        with pytest.raises(ValueError, match="positive"):
            transpile(
                "SELECT * FROM a JOIN b ON a.interval INTERSECTS b.interval",
                tables=["a", "b"],
                bin_size=-1,
            )

    def test_multi_join_all_intersects_rewritten(self):
        """
        GIVEN a three-way join with two INTERSECTS conditions
        WHEN transpiling
        THEN should create binned CTEs for all three tables, reusing the
             FROM table CTE, and place equi-join + overlap in each JOIN ON
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

        # Three distinct CTEs
        assert "__giql_a_binned" in sql
        assert "__giql_b_binned" in sql
        assert "__giql_c_binned" in sql

        # FROM table CTE created only once
        assert sql.count("FROM peaks") == 1 or sql.upper().count("FROM PEAKS") == 1

        # Both JOINs have equi-join + overlap in ON
        sql_upper = sql.upper()
        assert sql_upper.count("__GIQL_BIN") >= 4  # at least 2 per JOIN ON


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
            bin_size=10000,
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
