"""Tests for the transpile() function."""

import pytest

import giql
from giql import Table
from giql import transpile


class TestTranspileBasic:
    """Tests for basic transpilation with string table names."""

    def test_transpile_intersects_literal(self):
        """
        GIVEN a GIQL query with INTERSECTS and literal range
        WHEN transpiling with string table name
        THEN should return valid SQL with default column names
        """
        sql = transpile(
            "SELECT * FROM peaks WHERE interval INTERSECTS 'chr1:1000-2000'",
            tables=["peaks"],
        )

        assert "SELECT" in sql
        assert "peaks" in sql
        assert "chrom" in sql
        assert "start" in sql
        assert "end" in sql
        assert "chr1" in sql

    def test_transpile_contains_literal(self):
        """
        GIVEN a GIQL query with CONTAINS and literal point
        WHEN transpiling with string table name
        THEN should return valid SQL for point containment
        """
        sql = transpile(
            "SELECT * FROM peaks WHERE interval CONTAINS 'chr1:1500'",
            tables=["peaks"],
        )

        assert "SELECT" in sql
        assert "peaks" in sql
        assert "1500" in sql

    def test_transpile_within_literal(self):
        """
        GIVEN a GIQL query with WITHIN and literal range
        WHEN transpiling with string table name
        THEN should return valid SQL for interval within range
        """
        sql = transpile(
            "SELECT * FROM peaks WHERE interval WITHIN 'chr1:1000-2000'",
            tables=["peaks"],
        )

        assert "SELECT" in sql
        assert "peaks" in sql

    def test_transpile_no_tables(self):
        """
        GIVEN a GIQL query with INTERSECTS
        WHEN transpiling with no tables parameter
        THEN should return valid SQL with default column names
        """
        sql = transpile(
            "SELECT * FROM peaks WHERE interval INTERSECTS 'chr1:1000-2000'",
        )

        assert "SELECT" in sql
        assert "peaks" in sql


class TestTranspileWithTableObjects:
    """Tests for transpilation with Table objects."""

    def test_transpile_custom_columns(self):
        """
        GIVEN a GIQL query with INTERSECTS
        WHEN transpiling with custom column mappings
        THEN should use custom column names in generated SQL
        """
        sql = transpile(
            "SELECT * FROM peaks WHERE interval INTERSECTS 'chr1:1000-2000'",
            tables=[
                Table(
                    "peaks",
                    genomic_col="interval",
                    chrom_col="chromosome",
                    start_col="start_pos",
                    end_col="end_pos",
                )
            ],
        )

        assert "SELECT" in sql
        assert "peaks" in sql
        assert '"chromosome"' in sql
        assert '"start_pos"' in sql
        assert '"end_pos"' in sql
        # Should NOT contain default column names
        assert '"chrom"' not in sql
        assert '"start"' not in sql
        assert '"end"' not in sql

    def test_transpile_no_strand_column(self):
        """
        GIVEN a Table with strand_col=None
        WHEN transpiling a query
        THEN should not require strand column
        """
        sql = transpile(
            "SELECT * FROM peaks WHERE interval INTERSECTS 'chr1:1000-2000'",
            tables=[Table("peaks", strand_col=None)],
        )

        assert "SELECT" in sql


class TestTranspileMultipleTables:
    """Tests for transpilation with multiple tables."""

    def test_transpile_join_intersects(self):
        """
        GIVEN a GIQL query joining two tables with INTERSECTS
        WHEN transpiling with both tables configured
        THEN should generate correct join conditions
        """
        sql = transpile(
            """
            SELECT a.*, b.*
            FROM peaks a
            JOIN genes b ON a.interval INTERSECTS b.region
            """,
            tables=[
                Table("peaks", genomic_col="interval"),
                Table("genes", genomic_col="region"),
            ],
        )

        assert "SELECT" in sql
        assert "peaks" in sql
        assert "genes" in sql
        assert "JOIN" in sql.upper()

    def test_transpile_different_schemas(self):
        """
        GIVEN two tables with different column schemas
        WHEN transpiling a join query
        THEN should use correct columns for each table
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

        assert "SELECT" in sql
        # Both table's column names should appear
        assert "chromosome" in sql or "start_pos" in sql
        assert "seqname" in sql or "begin" in sql or "terminus" in sql


class TestTranspileSpatialOperators:
    """Tests for all spatial operators."""

    def test_intersects_any(self):
        """
        GIVEN a GIQL query with INTERSECTS ANY
        WHEN transpiling
        THEN should generate OR conditions for multiple ranges
        """
        sql = transpile(
            """
            SELECT * FROM peaks
            WHERE interval INTERSECTS ANY('chr1:1000-2000', 'chr2:500-1000')
            """,
            tables=["peaks"],
        )

        assert "SELECT" in sql
        assert "chr1" in sql
        assert "chr2" in sql
        assert " OR " in sql

    def test_intersects_all(self):
        """
        GIVEN a GIQL query with INTERSECTS ALL
        WHEN transpiling
        THEN should generate AND conditions for multiple ranges
        """
        sql = transpile(
            """
            SELECT * FROM peaks
            WHERE interval INTERSECTS ALL('chr1:1000-2000', 'chr1:1500-2500')
            """,
            tables=["peaks"],
        )

        assert "SELECT" in sql
        assert " AND " in sql


class TestTranspileCluster:
    """Tests for CLUSTER operation."""

    def test_cluster_basic(self):
        """
        GIVEN a GIQL query with CLUSTER
        WHEN transpiling
        THEN should generate window function for clustering
        """
        sql = transpile(
            """
            SELECT *, CLUSTER(interval) AS cluster_id
            FROM peaks
            """,
            tables=["peaks"],
        )

        assert "SELECT" in sql
        assert "SUM" in sql.upper() or "LAG" in sql.upper()

    def test_cluster_with_distance(self):
        """
        GIVEN a GIQL query with CLUSTER and distance parameter
        WHEN transpiling
        THEN should include distance in clustering logic
        """
        sql = transpile(
            """
            SELECT *, CLUSTER(interval, 100) AS cluster_id
            FROM peaks
            """,
            tables=["peaks"],
        )

        assert "SELECT" in sql
        assert "100" in sql

    def test_cluster_stranded(self):
        """
        GIVEN a GIQL query with stranded CLUSTER
        WHEN transpiling
        THEN should partition by strand
        """
        sql = transpile(
            """
            SELECT *, CLUSTER(interval, stranded := true) AS cluster_id
            FROM peaks
            """,
            tables=["peaks"],
        )

        assert "SELECT" in sql
        assert "strand" in sql.lower()


class TestTranspileMerge:
    """Tests for MERGE operation."""

    def test_merge_basic(self):
        """
        GIVEN a GIQL query with MERGE
        WHEN transpiling
        THEN should generate GROUP BY with MIN/MAX aggregation
        """
        sql = transpile(
            "SELECT MERGE(interval) FROM peaks",
            tables=["peaks"],
        )

        assert "SELECT" in sql
        assert "MIN" in sql.upper()
        assert "MAX" in sql.upper()
        assert "GROUP BY" in sql.upper()

    def test_merge_with_distance(self):
        """
        GIVEN a GIQL query with MERGE and distance parameter
        WHEN transpiling
        THEN should include distance in merge logic
        """
        sql = transpile(
            "SELECT MERGE(interval, 100) FROM peaks",
            tables=["peaks"],
        )

        assert "SELECT" in sql
        assert "100" in sql

    def test_merge_with_aggregation(self):
        """
        GIVEN a GIQL query with MERGE and additional aggregation
        WHEN transpiling
        THEN should include both merge and custom aggregation
        """
        sql = transpile(
            "SELECT MERGE(interval), COUNT(*) as count FROM peaks",
            tables=["peaks"],
        )

        assert "SELECT" in sql
        assert "COUNT" in sql.upper()


class TestTranspileNearest:
    """Tests for NEAREST operation."""

    def test_nearest_standalone(self):
        """
        GIVEN a GIQL query with standalone NEAREST
        WHEN transpiling
        THEN should generate subquery with ORDER BY and LIMIT
        """
        sql = transpile(
            "SELECT * FROM NEAREST(genes, reference := 'chr1:1000-2000', k := 3)",
            tables=["genes"],
        )

        assert "SELECT" in sql
        assert "ORDER BY" in sql.upper()
        assert "LIMIT 3" in sql

    def test_nearest_with_max_distance(self):
        """
        GIVEN a GIQL query with NEAREST and max_distance
        WHEN transpiling
        THEN should include distance filter
        """
        sql = transpile(
            """
            SELECT * FROM NEAREST(genes, reference := 'chr1:1000-2000', k := 5, max_distance := 100000)
            """,
            tables=["genes"],
        )

        assert "SELECT" in sql
        assert "100000" in sql
        assert "LIMIT 5" in sql

    def test_nearest_lateral(self):
        """
        GIVEN a GIQL query with NEAREST in LATERAL join
        WHEN transpiling
        THEN should generate LATERAL subquery
        """
        sql = transpile(
            """
            SELECT *
            FROM peaks
            CROSS JOIN LATERAL NEAREST(genes, reference := peaks.interval, k := 3)
            """,
            tables=["peaks", "genes"],
        )

        assert "SELECT" in sql
        assert "LATERAL" in sql.upper()
        assert "LIMIT 3" in sql


class TestTranspileErrors:
    """Tests for error handling."""

    def test_invalid_syntax(self):
        """
        GIVEN an invalid GIQL query
        WHEN transpiling
        THEN should raise ValueError with parse error
        """
        with pytest.raises(ValueError, match="Parse error"):
            transpile("SELECT * FORM peaks")  # typo: FORM instead of FROM

    def test_unknown_dialect(self):
        """Test transpilation with an unrecognized dialect.

        Given:
            A GIQL query and a dialect that is not a supported target.
        When:
            Transpiling.
        Then:
            It should raise ValueError naming the supported targets.
        """
        # Act & assert
        with pytest.raises(ValueError, match="Unknown dialect: 'postgres'"):
            transpile(
                "SELECT * FROM peaks WHERE interval INTERSECTS 'chr1:1000-2000'",
                tables=["peaks"],
                dialect="postgres",  # type: ignore[call-overload]
            )


class TestTranspileDialects:
    """Tests for dialect/target selection."""

    def test_transpile_datafusion_accepted(self):
        """Test that the datafusion dialect is a valid target.

        Given:
            A GIQL query and dialect="datafusion".
        When:
            Transpiling.
        Then:
            It should return SQL referencing the table without raising.
        """
        # Act
        sql = transpile(
            "SELECT * FROM peaks WHERE interval INTERSECTS 'chr1:1000-2000'",
            tables=["peaks"],
            dialect="datafusion",
        )

        # Assert
        assert "SELECT" in sql
        assert "peaks" in sql

    @pytest.mark.parametrize(
        "query, tables",
        [
            (
                "SELECT * FROM peaks WHERE interval INTERSECTS 'chr1:1000-2000'",
                ["peaks"],
            ),
            ("SELECT * FROM peaks WHERE interval CONTAINS 'chr1:1500'", ["peaks"]),
            (
                "SELECT * FROM peaks WHERE interval WITHIN 'chr1:1000-2000'",
                ["peaks"],
            ),
            (
                "SELECT * FROM NEAREST(genes, reference := 'chr1:1000-2000', k := 3)",
                ["genes"],
            ),
            (
                "SELECT a.chrom, b.chrom FROM peaks a "
                "JOIN genes b ON a.interval INTERSECTS b.interval",
                ["peaks", "genes"],
            ),
            (
                "SELECT * FROM peaks WHERE interval INTERSECTS "
                "ANY('chr1:1000-2000', 'chr1:5000-6000')",
                ["peaks"],
            ),
            (
                "SELECT DISTANCE(a.interval, b.interval) AS d "
                "FROM peaks a, genes b",
                ["peaks", "genes"],
            ),
        ],
        ids=[
            "intersects_literal",
            "contains",
            "within",
            "nearest",
            "join",
            "any",
            "distance",
        ],
    )
    def test_transpile_datafusion_matches_generic_output(self, query, tables):
        """Test that datafusion is currently a pure alias for the generic target.

        Given:
            A GIQL query transpiled with dialect=None and with
            dialect="datafusion". The operator spread is a regression net,
            not codepath coverage: datafusion has no engine-specific
            expansion yet, so every query funnels through the same generic
            path.
        When:
            Comparing the two outputs.
        Then:
            It should produce byte-identical SQL for every query, pinning
            datafusion as a generic alias until later steps diverge it.
        """
        # Act
        generic_sql = transpile(query, tables=tables)
        datafusion_sql = transpile(query, tables=tables, dialect="datafusion")

        # Assert
        assert datafusion_sql == generic_sql

    @pytest.mark.parametrize(
        "query, tables",
        [
            (
                "SELECT DISTANCE(a.interval, b.interval) AS d FROM peaks a, genes b",
                ["peaks", "genes"],
            ),
            (
                "SELECT DISTANCE(a.interval, b.interval, stranded := true) AS d "
                "FROM peaks a, genes b",
                ["peaks", "genes"],
            ),
            (
                "SELECT DISTANCE(a.interval, b.interval, signed := true) AS d "
                "FROM peaks a, genes b",
                ["peaks", "genes"],
            ),
            (
                "SELECT DISTANCE(a.interval, b.interval, stranded := true, "
                "signed := true) AS d FROM peaks a, genes b",
                ["peaks", "genes"],
            ),
        ],
        ids=["unsigned", "stranded", "signed", "stranded_signed"],
    )
    def test_transpile_distance_is_byte_identical_across_targets(self, query, tables):
        """Test that the migrated DISTANCE operator emits identical SQL everywhere.

        Given:
            A DISTANCE query in each of the four shapes (unsigned/signed x
            non-stranded/stranded), which migrated onto the AST-expansion pass
            with a single generic expander.
        When:
            Transpiling it with dialect=None, "duckdb", and "datafusion".
        Then:
            The three outputs should be byte-identical — the generic expander
            covers every target and DISTANCE has no engine-specific divergence.
        """
        # Act
        generic_sql = transpile(query, tables=tables, dialect=None)
        duckdb_sql = transpile(query, tables=tables, dialect="duckdb")
        datafusion_sql = transpile(query, tables=tables, dialect="datafusion")

        # Assert
        assert generic_sql == duckdb_sql == datafusion_sql

    def test_transpile_datafusion_accepts_intersects_bin_size(self):
        """Test that datafusion honours the binned-join bin size identically.

        Given:
            A column-to-column INTERSECTS join with an explicit
            intersects_bin_size, transpiled with dialect=None and
            dialect="datafusion".
        When:
            Comparing the two outputs.
        Then:
            It should produce byte-identical SQL — datafusion takes the same
            binned-join path and honours the bin size.
        """
        # Arrange
        query = (
            "SELECT a.chrom, b.chrom FROM peaks a "
            "JOIN genes b ON a.interval INTERSECTS b.interval"
        )

        # Act
        generic_sql = transpile(
            query, tables=["peaks", "genes"], intersects_bin_size=50000
        )
        datafusion_sql = transpile(
            query,
            tables=["peaks", "genes"],
            dialect="datafusion",
            intersects_bin_size=50000,
        )

        # Assert
        assert datafusion_sql == generic_sql

    def test_transpile_duckdb_returns_sql_for_column_join(self):
        """Test the duckdb happy path for a column-to-column INTERSECTS join.

        Given:
            A column-to-column INTERSECTS join and dialect="duckdb".
        When:
            Transpiling.
        Then:
            It should return SQL referencing the joined tables without raising.
        """
        # Act
        sql = transpile(
            "SELECT a.chrom, a.start, b.start FROM peaks a "
            "JOIN genes b ON a.interval INTERSECTS b.interval",
            tables=["peaks", "genes"],
            dialect="duckdb",
        )

        # Assert
        assert "SELECT" in sql
        assert "peaks" in sql
        assert "genes" in sql

    def test_transpile_duckdb_rejects_intersects_bin_size(self):
        """Test that duckdb and intersects_bin_size are mutually exclusive.

        Given:
            dialect="duckdb" combined with an explicit intersects_bin_size.
        When:
            Transpiling.
        Then:
            It should raise ValueError explaining the IEJoin plan ignores
            the bin size.
        """
        # Act & assert
        with pytest.raises(
            ValueError,
            match=r"intersects_bin_size has no effect.*Pass one or the other, not both\.",
        ):
            transpile(
                "SELECT a.chrom, b.chrom FROM peaks a "
                "JOIN genes b ON a.interval INTERSECTS b.interval",
                tables=["peaks", "genes"],
                dialect="duckdb",
                intersects_bin_size=50000,  # type: ignore[call-overload]
            )

    def test_transpile_duckdb_rejects_zero_intersects_bin_size(self):
        """Test that a falsy-but-set bin size still triggers the duckdb guard.

        Given:
            dialect="duckdb" combined with intersects_bin_size=0.
        When:
            Transpiling.
        Then:
            It should raise ValueError — the guard is `is not None`, so 0
            (falsy) still conflicts with the IEJoin plan.
        """
        # Act & assert
        with pytest.raises(ValueError, match="intersects_bin_size has no effect"):
            transpile(
                "SELECT a.chrom, b.chrom FROM peaks a "
                "JOIN genes b ON a.interval INTERSECTS b.interval",
                tables=["peaks", "genes"],
                dialect="duckdb",
                intersects_bin_size=0,  # type: ignore[call-overload]
            )


class TestModuleExports:
    """Tests for module-level exports."""

    def test_transpile_exported(self):
        """
        GIVEN the giql module
        WHEN accessing transpile
        THEN should be available at module level
        """
        assert hasattr(giql, "transpile")
        assert callable(giql.transpile)

    def test_table_exported(self):
        """
        GIVEN the giql module
        WHEN accessing Table
        THEN should be available at module level
        """
        assert hasattr(giql, "Table")
        assert giql.Table is Table
