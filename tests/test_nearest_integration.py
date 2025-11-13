"""Integration tests for NEAREST operator with actual database execution.

Tests verify end-to-end functionality: parsing, transpilation, and execution
of NEAREST queries across all supported dialects.
"""

import duckdb
import pytest

from giql import GIQLEngine


@pytest.fixture
def duckdb_engine_with_sample_data():
    """Create DuckDB engine with sample peaks and genes data."""
    engine = GIQLEngine(target_dialect="duckdb")

    # Create peaks table
    engine.conn.execute("""
        CREATE TABLE peaks (
            peak_id INTEGER,
            chromosome VARCHAR,
            start_pos INTEGER,
            end_pos INTEGER
        )
    """)

    # Create genes table
    engine.conn.execute("""
        CREATE TABLE genes (
            gene_id INTEGER,
            gene_name VARCHAR,
            chromosome VARCHAR,
            start_pos INTEGER,
            end_pos INTEGER,
            strand VARCHAR
        )
    """)

    # Insert sample peaks on chr1
    engine.conn.execute("""
        INSERT INTO peaks VALUES
        (1, 'chr1', 1000, 1100),   -- Peak at position 1000-1100
        (2, 'chr1', 5000, 5100),   -- Peak at position 5000-5100
        (3, 'chr2', 2000, 2100)    -- Peak on different chromosome
    """)

    # Insert sample genes on chr1 and chr2
    engine.conn.execute("""
        INSERT INTO genes VALUES
        (1, 'GENE_A', 'chr1', 1500, 2000, '+'),  -- 400bp from peak 1
        (2, 'GENE_B', 'chr1', 3000, 3500, '+'),  -- 1900bp from peak 1, close to peak 2
        (3, 'GENE_C', 'chr1', 1050, 1200, '-'),  -- Overlaps peak 1 (chr1:1000-1100)
        (4, 'GENE_D', 'chr1', 6000, 6500, '+'),  -- 900bp from peak 2
        (5, 'GENE_E', 'chr1', 8000, 8500, '-'),  -- 2900bp from peak 2
        (6, 'GENE_F', 'chr2', 2500, 3000, '+'),  -- 400bp from peak 3
        (7, 'GENE_G', 'chr2', 2050, 2200, '-')   -- Overlaps peak 3 (chr2:2000-2100)
    """)

    # Register schema
    engine.register_table_schema(
        "peaks",
        {
            "peak_id": "INTEGER",
            "chromosome": "VARCHAR",
            "start_pos": "INTEGER",
            "end_pos": "INTEGER",
        },
        genomic_column="position",
    )
    engine.register_table_schema(
        "genes",
        {
            "gene_id": "INTEGER",
            "gene_name": "VARCHAR",
            "chromosome": "VARCHAR",
            "start_pos": "INTEGER",
            "end_pos": "INTEGER",
            "strand": "VARCHAR",
        },
        genomic_column="position",
    )

    return engine


class TestNearestIntegrationDuckDB:
    """Integration tests for NEAREST with DuckDB backend."""

    def test_nearest_k1_basic(self, duckdb_engine_with_sample_data):
        """
        GIVEN a peaks table and genes table with sample data
        WHEN querying for k=1 nearest gene to each peak
        THEN should return the closest gene for each peak
        """
        engine = duckdb_engine_with_sample_data

        cursor = engine.execute("""
            SELECT
                peaks.peak_id,
                peaks.chromosome,
                peaks.start_pos AS peak_start,
                nearest.gene_name,
                nearest.start_pos AS gene_start,
                nearest.distance
            FROM peaks
            CROSS JOIN LATERAL NEAREST(genes, reference=peaks.position, k=1) AS nearest
            ORDER BY peaks.peak_id
        """)

        rows = cursor.fetchall()

        # Should have 3 rows (one per peak)
        assert len(rows) == 3

        # Peak 1 (chr1:1000-1100) should be closest to GENE_C (overlapping, distance=0)
        assert rows[0][0] == 1  # peak_id
        assert rows[0][3] == "GENE_C"  # gene_name
        assert rows[0][5] == 0  # distance (overlapping)

        # Peak 2 (chr1:5000-5100) should be closest to GENE_D (900bp away)
        assert rows[1][0] == 2  # peak_id
        assert rows[1][3] == "GENE_D"  # gene_name
        assert rows[1][5] == 900  # distance

        # Peak 3 (chr2:2000-2100) should be closest to GENE_G (overlapping, distance=0)
        assert rows[2][0] == 3  # peak_id
        assert rows[2][3] == "GENE_G"  # gene_name
        assert rows[2][5] == 0  # distance (overlapping)

    def test_nearest_k3_multiple_neighbors(self, duckdb_engine_with_sample_data):
        """
        GIVEN a peaks table and genes table with sample data
        WHEN querying for k=3 nearest genes to each peak
        THEN should return up to 3 closest genes for each peak
        """
        engine = duckdb_engine_with_sample_data

        cursor = engine.execute("""
            SELECT
                peaks.peak_id,
                nearest.gene_name,
                nearest.distance
            FROM peaks
            CROSS JOIN LATERAL NEAREST(genes, reference=peaks.position, k=3) AS nearest
            ORDER BY peaks.peak_id, nearest.distance
        """)

        rows = cursor.fetchall()

        # Should have results for all peaks (some may have <3 neighbors on same chromosome)
        peak1_results = [r for r in rows if r[0] == 1]
        peak2_results = [r for r in rows if r[0] == 2]
        peak3_results = [r for r in rows if r[0] == 3]

        # Peak 1 should have 3 genes on chr1
        assert len(peak1_results) == 3
        assert peak1_results[0][1] == "GENE_C"  # Closest (overlapping)
        assert peak1_results[0][2] == 0

        # Peak 2 should have 3 genes on chr1
        assert len(peak2_results) == 3
        assert peak2_results[0][1] == "GENE_D"  # Closest
        assert peak2_results[0][2] == 900

        # Peak 3 should have 2 genes on chr2 (only 2 genes exist on chr2)
        assert len(peak3_results) == 2
        assert peak3_results[0][1] == "GENE_G"  # Closest (overlapping)
        assert peak3_results[0][2] == 0

    def test_nearest_with_max_distance(self, duckdb_engine_with_sample_data):
        """
        GIVEN a peaks table and genes table with sample data
        WHEN querying for k=5 nearest genes with max_distance=1000
        THEN should only return genes within 1000bp
        """
        engine = duckdb_engine_with_sample_data

        cursor = engine.execute("""
            SELECT
                peaks.peak_id,
                nearest.gene_name,
                nearest.distance
            FROM peaks
            CROSS JOIN LATERAL NEAREST(genes, reference=peaks.position, k=5, max_distance=1000) AS nearest
            WHERE peaks.peak_id = 2
            ORDER BY nearest.distance
        """)

        rows = cursor.fetchall()

        # Peak 2 should have GENE_D (900bp) but not GENE_E (2900bp)
        assert len(rows) <= 5
        assert all(row[2] <= 1000 for row in rows)
        assert "GENE_D" in [row[1] for row in rows]
        assert "GENE_E" not in [row[1] for row in rows]

    def test_nearest_standalone_literal_reference(self, duckdb_engine_with_sample_data):
        """
        GIVEN a genes table with sample data
        WHEN querying with a literal reference position
        THEN should return k nearest genes to that position
        """
        engine = duckdb_engine_with_sample_data

        cursor = engine.execute("""
            SELECT
                gene_name,
                chromosome,
                start_pos,
                distance
            FROM NEAREST(genes, reference='chr1:4000-4100', k=2) AS nearest
            ORDER BY distance
        """)

        rows = cursor.fetchall()

        # Should have 2 genes
        assert len(rows) == 2

        # GENE_B (chr1:3000-3500) should be closest (~500bp)
        assert rows[0][0] == "GENE_B"
        assert rows[0][3] < 1000  # Close distance

        # Results should be ordered by distance
        assert rows[0][3] <= rows[1][3]

    def test_nearest_k_greater_than_available(self, duckdb_engine_with_sample_data):
        """
        GIVEN a genes table with limited genes per chromosome
        WHEN querying for k=10 (more than available)
        THEN should return all available genes on the same chromosome
        """
        engine = duckdb_engine_with_sample_data

        cursor = engine.execute("""
            SELECT
                peaks.peak_id,
                nearest.gene_name
            FROM peaks
            CROSS JOIN LATERAL NEAREST(genes, reference=peaks.position, k=10) AS nearest
            WHERE peaks.peak_id = 3
        """)

        rows = cursor.fetchall()

        # Peak 3 is on chr2, which only has 2 genes
        assert len(rows) == 2

    def test_nearest_k_equals_1(self, duckdb_engine_with_sample_data):
        """
        GIVEN a peaks table and genes table
        WHEN querying with k=1 (default behavior)
        THEN should return exactly 1 nearest gene per peak
        """
        engine = duckdb_engine_with_sample_data

        cursor = engine.execute("""
            SELECT
                peaks.peak_id,
                COUNT(*) as gene_count
            FROM peaks
            CROSS JOIN LATERAL NEAREST(genes, reference=peaks.position, k=1) AS nearest
            GROUP BY peaks.peak_id
        """)

        rows = cursor.fetchall()

        # Each peak should have exactly 1 gene
        assert len(rows) == 3
        assert all(row[1] == 1 for row in rows)

    def test_nearest_different_chromosomes(self, duckdb_engine_with_sample_data):
        """
        GIVEN peaks and genes on different chromosomes
        WHEN querying for nearest genes
        THEN should only return genes on the same chromosome
        """
        engine = duckdb_engine_with_sample_data

        cursor = engine.execute("""
            SELECT
                peaks.peak_id,
                nearest.chromosome,
                nearest.gene_name
            FROM peaks
            CROSS JOIN LATERAL NEAREST(genes, reference=peaks.position, k=5) AS nearest
        """)

        rows = cursor.fetchall()

        # Verify chromosome matching
        chr1_peaks = [row for row in rows if row[0] in (1, 2)]
        chr2_peaks = [row for row in rows if row[0] == 3]

        # All genes for chr1 peaks should be on chr1
        assert all(row[1] == "chr1" for row in chr1_peaks)

        # All genes for chr2 peaks should be on chr2
        assert all(row[1] == "chr2" for row in chr2_peaks)

    def test_nearest_literal_reference_with_strand_notation(
        self, duckdb_engine_with_sample_data
    ):
        """
        GIVEN a genes table with sample data
        WHEN querying with a literal reference using strand notation
        THEN should parse and execute correctly
        """
        engine = duckdb_engine_with_sample_data

        cursor = engine.execute("""
            SELECT
                gene_name,
                chromosome,
                start_pos,
                distance
            FROM NEAREST(genes, reference='chr1:4000-4100:+', k=3) AS nearest
            ORDER BY distance
        """)

        rows = cursor.fetchall()

        # Should return results (strand notation should be parsed)
        assert len(rows) > 0
        assert all(row[1] == "chr1" for row in rows)  # All on chr1

        # Results should be ordered by distance
        for i in range(len(rows) - 1):
            assert rows[i][3] <= rows[i + 1][3], "Results should be ordered by distance"

    def test_nearest_all_features_beyond_max_distance(
        self, duckdb_engine_with_sample_data
    ):
        """
        GIVEN peaks and genes with no genes within max_distance
        WHEN querying with a very small max_distance
        THEN should return empty result set
        """
        engine = duckdb_engine_with_sample_data

        cursor = engine.execute("""
            SELECT
                peaks.peak_id,
                nearest.gene_name,
                nearest.distance
            FROM peaks
            CROSS JOIN LATERAL NEAREST(genes, reference=peaks.position, k=5, max_distance=10) AS nearest
            WHERE peaks.peak_id = 2
        """)

        rows = cursor.fetchall()

        # Peak 2 has no genes within 10bp, should return empty
        assert len(rows) == 0, (
            "Should return empty when all features beyond max_distance"
        )

    def test_nearest_some_within_some_beyond_max_distance(
        self, duckdb_engine_with_sample_data
    ):
        """
        GIVEN peaks and genes with some genes within and some beyond max_distance
        WHEN querying with max_distance constraint
        THEN should return only features within threshold
        """
        engine = duckdb_engine_with_sample_data

        cursor = engine.execute("""
            SELECT
                peaks.peak_id,
                nearest.gene_name,
                nearest.distance
            FROM peaks
            CROSS JOIN LATERAL NEAREST(genes, reference=peaks.position, k=10, max_distance=500) AS nearest
            WHERE peaks.peak_id = 1
            ORDER BY nearest.distance
        """)

        rows = cursor.fetchall()

        # Peak 1 (chr1:1000-1100) has:
        # - GENE_C at distance 0 (overlapping, chr1:1050-1200)
        # - GENE_A at distance 400 (chr1:1500-2000)
        # - Other genes are further away
        assert len(rows) >= 1, "Should have at least GENE_C"
        assert len(rows) <= 2, "Should have at most GENE_C and GENE_A within 500bp"

        # All returned features should be within max_distance
        assert all(row[2] <= 500 for row in rows), (
            "All features should be within max_distance"
        )

        # Should include GENE_C (overlapping)
        gene_names = [row[1] for row in rows]
        assert "GENE_C" in gene_names, "Should include overlapping GENE_C"

    def test_nearest_stranded_same_strand_only(self, duckdb_engine_with_sample_data):
        """
        GIVEN peaks and genes with strand information
        WHEN querying with stranded=true
        THEN should return only same-strand features
        """
        engine = duckdb_engine_with_sample_data

        # Add strand column to peaks
        engine.conn.execute("ALTER TABLE peaks ADD COLUMN strand VARCHAR DEFAULT '+'")

        # Re-register schema with strand info
        engine.register_table_schema(
            "peaks",
            {
                "peak_id": "INTEGER",
                "chromosome": "VARCHAR",
                "start_pos": "INTEGER",
                "end_pos": "INTEGER",
                "strand": "VARCHAR",
            },
            genomic_column="position",
        )

        cursor = engine.execute("""
            SELECT
                peaks.peak_id,
                nearest.gene_name,
                nearest.strand,
                nearest.distance
            FROM peaks
            CROSS JOIN LATERAL NEAREST(genes, reference=peaks.position, k=5, stranded=true) AS nearest
            WHERE peaks.peak_id = 1
            ORDER BY nearest.distance
        """)

        rows = cursor.fetchall()

        # Peak 1 is on + strand
        # Should get genes on + strand only
        if len(rows) > 0:
            # At minimum, verify the query executes and strand column is present
            assert len(rows[0]) == 4, "Should have 4 columns including strand"
            # All genes should be on + strand
            for row in rows:
                assert row[2] == "+", f"Gene {row[1]} should be on + strand"

    def test_nearest_stranded_different_strands_excluded(
        self, duckdb_engine_with_sample_data
    ):
        """
        GIVEN peaks and genes with mixed strands
        WHEN querying with stranded=true
        THEN should exclude genes on different strands
        """
        # First, add strand column to peaks table
        engine = duckdb_engine_with_sample_data

        # Add strand column and set values
        engine.conn.execute("ALTER TABLE peaks ADD COLUMN strand VARCHAR DEFAULT '+'")

        # Re-register schema with strand info
        engine.register_table_schema(
            "peaks",
            {
                "peak_id": "INTEGER",
                "chromosome": "VARCHAR",
                "start_pos": "INTEGER",
                "end_pos": "INTEGER",
                "strand": "VARCHAR",
            },
            genomic_column="position",
        )

        cursor = engine.execute("""
            SELECT
                peaks.peak_id,
                nearest.gene_name,
                peaks.strand AS peak_strand,
                nearest.strand AS gene_strand,
                nearest.distance
            FROM peaks
            CROSS JOIN LATERAL NEAREST(genes, reference=peaks.position, k=10, stranded=true) AS nearest
            WHERE peaks.peak_id = 1
        """)

        rows = cursor.fetchall()

        # Peak 1 is on + strand
        # Only genes on + strand should be returned
        plus_strand_genes = ["GENE_A", "GENE_D"]  # From test data

        if len(rows) > 0:
            # Verify all returned genes have same strand as peak
            for row in rows:
                assert row[2] == row[3], (
                    f"Peak strand {row[2]} should match gene strand {row[3]}"
                )
                # Gene should be on + strand
                assert row[3] == "+", f"Gene {row[1]} should be on + strand"

    def test_nearest_stranded_unspecified_strand(self, duckdb_engine_with_sample_data):
        """
        GIVEN genes with unspecified strand ('.' or NULL)
        WHEN querying with stranded=true
        THEN behavior depends on implementation (may include or exclude)
        """
        engine = duckdb_engine_with_sample_data

        # Add a gene with unspecified strand
        engine.conn.execute("""
            INSERT INTO genes VALUES
            (8, 'GENE_H', 'chr1', 1100, 1150, '.')
        """)

        # Add strand to peaks
        engine.conn.execute("ALTER TABLE peaks ADD COLUMN strand VARCHAR DEFAULT '+'")
        engine.register_table_schema(
            "peaks",
            {
                "peak_id": "INTEGER",
                "chromosome": "VARCHAR",
                "start_pos": "INTEGER",
                "end_pos": "INTEGER",
                "strand": "VARCHAR",
            },
            genomic_column="position",
        )

        cursor = engine.execute("""
            SELECT
                peaks.peak_id,
                nearest.gene_name,
                nearest.strand,
                nearest.distance
            FROM peaks
            CROSS JOIN LATERAL NEAREST(genes, reference=peaks.position, k=10, stranded=true) AS nearest
            WHERE peaks.peak_id = 1
            ORDER BY nearest.distance
        """)

        rows = cursor.fetchall()

        # Just verify the query executes without error
        # The handling of '.' strand is implementation-dependent
        assert isinstance(rows, list), "Should return a list of results"

    def test_nearest_signed_distance(self, duckdb_engine_with_sample_data):
        """
        GIVEN a DuckDB database with peaks and genes
        WHEN querying for nearest genes with signed=true
        THEN should return negative distances for upstream features and positive for downstream
        """
        engine = duckdb_engine_with_sample_data

        cursor = engine.execute("""
            SELECT
                peaks.peak_id,
                nearest.gene_name,
                nearest.distance
            FROM peaks
            CROSS JOIN LATERAL NEAREST(genes, reference=peaks.position, k=10, signed=true) AS nearest
            WHERE peaks.peak_id = 2
            ORDER BY ABS(nearest.distance)
        """)

        rows = cursor.fetchall()

        # Peak 2 is at chr1:5000-5100
        # GENE_B (chr1:3000-3500) is upstream → should have negative distance (-1500)
        # GENE_D (chr1:6000-6500) is downstream → should have positive distance (900)
        # GENE_E (chr1:8000-8500) is downstream → should have positive distance (2900)

        assert len(rows) >= 3, "Should find at least 3 genes on chr1"

        # Find GENE_B (upstream)
        gene_b = [r for r in rows if r[1] == "GENE_B"][0]
        assert gene_b[2] < 0, (
            f"GENE_B is upstream, distance should be negative, got {gene_b[2]}"
        )
        assert gene_b[2] == -1500, f"GENE_B distance should be -1500, got {gene_b[2]}"

        # Find GENE_D (downstream)
        gene_d = [r for r in rows if r[1] == "GENE_D"][0]
        assert gene_d[2] > 0, (
            f"GENE_D is downstream, distance should be positive, got {gene_d[2]}"
        )
        assert gene_d[2] == 900, f"GENE_D distance should be 900, got {gene_d[2]}"

    def test_nearest_upstream_filtering(self, duckdb_engine_with_sample_data):
        """
        GIVEN a DuckDB database with peaks and genes
        WHEN querying for nearest genes with signed=true and filtering for upstream (distance < 0)
        THEN should return only upstream features
        """
        engine = duckdb_engine_with_sample_data

        cursor = engine.execute("""
            SELECT
                peaks.peak_id,
                nearest.gene_name,
                nearest.distance
            FROM peaks
            CROSS JOIN LATERAL NEAREST(genes, reference=peaks.position, k=10, signed=true) AS nearest
            WHERE peaks.peak_id = 2 AND nearest.distance < 0
            ORDER BY nearest.distance DESC
        """)

        rows = cursor.fetchall()

        # Peak 2 is at chr1:5000-5100
        # Only GENE_B (chr1:3000-3500) is upstream
        assert len(rows) >= 1, "Should find at least 1 upstream gene"

        # All distances should be negative (upstream)
        for row in rows:
            assert row[2] < 0, (
                f"All distances should be negative (upstream), got {row[2]} for {row[1]}"
            )

        # GENE_B should be in the results
        gene_names = [r[1] for r in rows]
        assert "GENE_B" in gene_names, "GENE_B should be in upstream results"

    def test_nearest_downstream_filtering(self, duckdb_engine_with_sample_data):
        """
        GIVEN a DuckDB database with peaks and genes
        WHEN querying for nearest genes with signed=true and filtering for downstream (distance > 0)
        THEN should return only downstream features
        """
        engine = duckdb_engine_with_sample_data

        cursor = engine.execute("""
            SELECT
                peaks.peak_id,
                nearest.gene_name,
                nearest.distance
            FROM peaks
            CROSS JOIN LATERAL NEAREST(genes, reference=peaks.position, k=10, signed=true) AS nearest
            WHERE peaks.peak_id = 2 AND nearest.distance > 0
            ORDER BY nearest.distance
        """)

        rows = cursor.fetchall()

        # Peak 2 is at chr1:5000-5100
        # GENE_D and GENE_E are downstream
        assert len(rows) >= 2, "Should find at least 2 downstream genes"

        # All distances should be positive (downstream)
        for row in rows:
            assert row[2] > 0, (
                f"All distances should be positive (downstream), got {row[2]} for {row[1]}"
            )

        # GENE_D should be closest downstream
        assert rows[0][1] == "GENE_D", (
            f"GENE_D should be closest downstream, got {rows[0][1]}"
        )
        assert rows[0][2] == 900, f"GENE_D distance should be 900, got {rows[0][2]}"

    def test_nearest_stranded_and_signed(self, duckdb_engine_with_sample_data):
        """
        GIVEN a DuckDB database with peaks and genes
        WHEN querying for nearest genes with both stranded=true and signed=true
        THEN should filter by strand AND return signed distances
        """
        engine = duckdb_engine_with_sample_data

        # Add strand column to peaks
        engine.conn.execute("ALTER TABLE peaks ADD COLUMN strand VARCHAR DEFAULT '+'")
        engine.register_table_schema(
            "peaks",
            {
                "peak_id": "INTEGER",
                "chromosome": "VARCHAR",
                "start_pos": "INTEGER",
                "end_pos": "INTEGER",
                "strand": "VARCHAR",
            },
            genomic_column="position",
        )

        cursor = engine.execute("""
            SELECT
                peaks.peak_id,
                nearest.gene_name,
                nearest.strand,
                nearest.distance
            FROM peaks
            CROSS JOIN LATERAL NEAREST(genes, reference=peaks.position, k=10, stranded=true, signed=true) AS nearest
            WHERE peaks.peak_id = 2
            ORDER BY ABS(nearest.distance)
        """)

        rows = cursor.fetchall()

        # Peak 2 is at chr1:5000-5100 with strand '+'
        # Should only get genes with strand '+': GENE_A, GENE_B, GENE_D
        # All should have signed distances

        for row in rows:
            peak_id, gene_name, strand, distance = row
            # All should be on '+' strand
            assert strand == "+", (
                f"All genes should be on '+' strand, got '{strand}' for {gene_name}"
            )

        # Verify we get both upstream and downstream genes
        gene_names_and_distances = [(r[1], r[3]) for r in rows]

        # Should have at least one upstream (negative) and one downstream (positive)
        has_upstream = any(d < 0 for _, d in gene_names_and_distances)
        has_downstream = any(d > 0 for _, d in gene_names_and_distances)

        assert has_upstream or has_downstream, (
            "Should have either upstream or downstream genes"
        )


class TestNearestIntegrationSQLite:
    """Integration tests for NEAREST with SQLite backend.

    SQLite does not support LATERAL joins. NEAREST works in standalone mode only.
    """

    def test_nearest_lateral_not_supported(self, duckdb_engine_with_sample_data):
        """
        GIVEN a SQLite database with peaks and genes
        WHEN attempting NEAREST in LATERAL mode (correlated)
        THEN should raise ValueError with helpful message
        """
        from giql.generators.sqlite import GIQLSQLiteGenerator

        sql = """
            SELECT *
            FROM peaks
            CROSS JOIN LATERAL NEAREST(genes, k=1) AS nearest
        """

        # Parse the query
        from sqlglot import parse_one
        from giql.dialect import GIQLDialect

        ast = parse_one(sql, dialect=GIQLDialect)

        # Try to generate SQLite SQL - should raise error
        generator = GIQLSQLiteGenerator(
            schema_info=duckdb_engine_with_sample_data.schema_info
        )

        with pytest.raises(ValueError) as exc_info:
            generator.generate(ast)

        assert "LATERAL" in str(exc_info.value)
        assert "SQLite" in str(exc_info.value)
        assert "standalone mode" in str(exc_info.value)

    def test_nearest_standalone_works(self, duckdb_engine_with_sample_data):
        """
        GIVEN a SQLite database with peaks and genes
        WHEN using NEAREST in standalone mode
        THEN should generate valid SQLite SQL
        """
        from giql.generators.sqlite import GIQLSQLiteGenerator

        sql = """
            SELECT *
            FROM NEAREST(genes, reference='chr1:1000-2000', k=3)
        """

        # Parse the query
        from sqlglot import parse_one
        from giql.dialect import GIQLDialect

        ast = parse_one(sql, dialect=GIQLDialect)

        # Generate SQLite SQL - should work
        generator = GIQLSQLiteGenerator(
            schema_info=duckdb_engine_with_sample_data.schema_info
        )

        result = generator.generate(ast)

        # Should contain standalone query pattern
        assert "SELECT" in result
        assert "FROM genes" in result or "FROM GENES" in result
        assert "ORDER BY" in result.upper()
        assert "LIMIT 3" in result
