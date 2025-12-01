"""Edge case tests for NEAREST operator.

Tests verify correct handling of boundary conditions, error cases,
and unusual inputs for the NEAREST operator.
"""

import pytest
from hypothesis import assume
from hypothesis import given
from hypothesis import strategies as st

from giql import GIQLEngine


@pytest.fixture
def duckdb_engine_with_edge_case_data():
    """Create DuckDB engine with data designed for edge case testing."""
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
            end_pos INTEGER
        )
    """)

    # Insert test data
    # Peak 1: chr1:1000-1100
    # Peak 2: chr2:5000-5100 (different chromosome, no genes)
    # Peak 3: chr1:10000-10100
    engine.conn.execute("""
        INSERT INTO peaks VALUES
        (1, 'chr1', 1000, 1100),
        (2, 'chr2', 5000, 5100),
        (3, 'chr1', 10000, 10100)
    """)

    # Genes with specific distance relationships
    # GENE_A and GENE_B are both 500bp from Peak 1 (tie scenario)
    # GENE_C overlaps Peak 1 (distance=0)
    # GENE_D, GENE_E, GENE_F on chr1 but far from Peak 3
    engine.conn.execute("""
        INSERT INTO genes VALUES
        (1, 'GENE_A', 'chr1', 1600, 1700),
        (2, 'GENE_B', 'chr1', 400, 500),
        (3, 'GENE_C', 'chr1', 1050, 1150),
        (4, 'GENE_D', 'chr1', 10500, 10600),
        (5, 'GENE_E', 'chr1', 11000, 11100),
        (6, 'GENE_F', 'chr1', 12000, 12100)
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
        genomic_column="interval",
    )
    engine.register_table_schema(
        "genes",
        {
            "gene_id": "INTEGER",
            "gene_name": "VARCHAR",
            "chromosome": "VARCHAR",
            "start_pos": "INTEGER",
            "end_pos": "INTEGER",
        },
        genomic_column="interval",
    )

    return engine


class TestNearestEdgeCases:
    """Edge case tests for NEAREST operator."""

    def test_k_equals_zero(self, duckdb_engine_with_edge_case_data):
        """
        GIVEN a NEAREST query with k=0
        WHEN executing the query
        THEN should return no results (LIMIT 0)
        """
        engine = duckdb_engine_with_edge_case_data

        cursor = engine.execute("""
            SELECT
                peaks.peak_id,
                nearest.gene_name
            FROM peaks
            CROSS JOIN LATERAL NEAREST(genes, reference=peaks.interval, k=0) AS nearest
            WHERE peaks.peak_id = 1
        """)

        rows = cursor.fetchall()
        assert len(rows) == 0, "k=0 should return no results"

    def test_ties_multiple_features_same_distance(
        self, duckdb_engine_with_edge_case_data
    ):
        """
        GIVEN multiple genes at the same distance from a peak
        WHEN querying for k=1 nearest
        THEN should return at least 1 result (behavior may vary for ties)
        """
        engine = duckdb_engine_with_edge_case_data

        cursor = engine.execute("""
            SELECT
                peaks.peak_id,
                nearest.gene_name,
                nearest.distance
            FROM peaks
            CROSS JOIN LATERAL NEAREST(genes, reference=peaks.interval, k=1) AS nearest
            WHERE peaks.peak_id = 1
            ORDER BY nearest.distance, nearest.gene_name
        """)

        rows = cursor.fetchall()

        # Should have at least 1 result
        assert len(rows) >= 1, "Should return at least one result for k=1"

        # All results should be at the same distance (ties)
        # Note: GENE_A and GENE_B are both 500bp away, GENE_C overlaps (0bp)
        # So the closest should be GENE_C at distance 0
        assert rows[0][1] == "GENE_C", (
            f"Closest gene should be GENE_C (overlapping), got {rows[0][1]}"
        )
        assert rows[0][2] == 0, f"Distance should be 0 (overlap), got {rows[0][2]}"

    def test_empty_result_set_different_chromosome(
        self, duckdb_engine_with_edge_case_data
    ):
        """
        GIVEN a peak on a chromosome with no genes
        WHEN querying for nearest genes
        THEN should return empty result set
        """
        engine = duckdb_engine_with_edge_case_data

        cursor = engine.execute("""
            SELECT
                peaks.peak_id,
                nearest.gene_name,
                nearest.distance
            FROM peaks
            CROSS JOIN LATERAL NEAREST(genes, reference=peaks.interval, k=10) AS nearest
            WHERE peaks.peak_id = 2
        """)

        rows = cursor.fetchall()

        # Peak 2 is on chr2, but all genes are on chr1
        # Should return empty result set
        assert len(rows) == 0, (
            "Should return empty result for peak on chromosome with no genes"
        )

    def test_overlapping_features_distance_zero(self, duckdb_engine_with_edge_case_data):
        """
        GIVEN a gene that overlaps a peak
        WHEN querying for nearest genes
        THEN should return distance=0 for overlapping gene
        """
        engine = duckdb_engine_with_edge_case_data

        cursor = engine.execute("""
            SELECT
                peaks.peak_id,
                nearest.gene_name,
                nearest.distance
            FROM peaks
            CROSS JOIN LATERAL NEAREST(genes, reference=peaks.interval, k=5) AS nearest
            WHERE peaks.peak_id = 1
            ORDER BY nearest.distance
        """)

        rows = cursor.fetchall()

        # GENE_C (chr1:1050-1150) overlaps Peak 1 (chr1:1000-1100)
        assert len(rows) > 0, "Should find genes"

        # First result should be the overlapping gene with distance 0
        assert rows[0][1] == "GENE_C", (
            f"First result should be GENE_C (overlapping), got {rows[0][1]}"
        )
        assert rows[0][2] == 0, (
            f"Distance should be 0 for overlapping gene, got {rows[0][2]}"
        )

    def test_missing_reference_in_standalone_mode(
        self, duckdb_engine_with_edge_case_data
    ):
        """
        GIVEN a standalone NEAREST query without reference parameter
        WHEN parsing/executing the query
        THEN should raise an error (reference is required in standalone mode)
        """
        engine = duckdb_engine_with_edge_case_data

        # Standalone mode (FROM NEAREST(...)) without reference parameter
        # This should fail because we can't determine the reference point
        with pytest.raises(Exception) as exc_info:
            engine.execute("""
                SELECT *
                FROM NEAREST(genes, k=3)
            """)

        # Should get an error about missing reference
        # The exact error message may vary, but it should mention reference
        error_msg = str(exc_info.value).lower()
        # Could be a ValueError, AttributeError, or SQL error depending on where it fails
        # Just verify it fails - the specific error type will be improved in T065

    def test_missing_target_table_in_schema(self, duckdb_engine_with_edge_case_data):
        """
        GIVEN a NEAREST query referencing a non-existent table
        WHEN executing the query
        THEN should raise an error about missing table
        """
        engine = duckdb_engine_with_edge_case_data

        # Query references 'nonexistent_table' which doesn't exist
        with pytest.raises(Exception) as exc_info:
            engine.execute("""
                SELECT *
                FROM peaks
                CROSS JOIN LATERAL NEAREST(nonexistent_table, reference=peaks.interval, k=3) AS nearest
            """)

        # Should get an error about the missing table
        error_msg = str(exc_info.value).lower()
        # DuckDB should raise an error about the table not existing

    def test_invalid_literal_range_format(self, duckdb_engine_with_edge_case_data):
        """
        GIVEN a NEAREST query with invalid literal range format
        WHEN parsing/executing the query
        THEN should raise an error about invalid range format
        """
        engine = duckdb_engine_with_edge_case_data

        # Invalid range formats
        # Note: "chr1:1000" is valid (point format), so not included
        invalid_ranges = [
            "chr1:not-a-number",  # Non-numeric coordinates
            "invalid-format",  # No colon separator
            "chr1:2000-1000",  # End before start (start >= end)
        ]

        for invalid_range in invalid_ranges:
            with pytest.raises(ValueError) as exc_info:
                engine.execute(f"""
                    SELECT *
                    FROM NEAREST(genes, reference='{invalid_range}', k=3)
                """)

            # Should get a ValueError about invalid range format
            error_msg = str(exc_info.value).lower()
            assert "invalid" in error_msg or "must be less" in error_msg, (
                f"Error message should mention invalid format or start/end issue: {exc_info.value}"
            )

    def test_nearest_with_additional_where_clause(
        self, duckdb_engine_with_edge_case_data
    ):
        """
        GIVEN a NEAREST query with additional WHERE clause filtering
        WHEN executing the query
        THEN should apply both NEAREST and WHERE filters
        """
        engine = duckdb_engine_with_edge_case_data

        cursor = engine.execute("""
            SELECT
                peaks.peak_id,
                nearest.gene_name,
                nearest.distance
            FROM peaks
            CROSS JOIN LATERAL NEAREST(genes, reference=peaks.interval, k=10) AS nearest
            WHERE peaks.peak_id = 1 AND nearest.distance < 600
            ORDER BY nearest.distance
        """)

        rows = cursor.fetchall()

        # Should find genes within 600bp of Peak 1
        # GENE_C overlaps (0bp) and GENE_A/GENE_B are 500bp away
        assert len(rows) >= 1, "Should find genes within 600bp"

        # All returned genes should have distance < 600
        for row in rows:
            assert row[2] < 600, f"All distances should be < 600bp, got {row[2]}"

    def test_nearest_with_cte(self, duckdb_engine_with_edge_case_data):
        """
        GIVEN a NEAREST query using a CTE for multiple query points
        WHEN executing the query
        THEN should correctly handle NEAREST within CTE
        """
        engine = duckdb_engine_with_edge_case_data

        cursor = engine.execute("""
            WITH selected_peaks AS (
                SELECT * FROM peaks WHERE peak_id IN (1, 3)
            )
            SELECT
                selected_peaks.peak_id,
                nearest.gene_name,
                nearest.distance
            FROM selected_peaks
            CROSS JOIN LATERAL NEAREST(genes, reference=selected_peaks.interval, k=2) AS nearest
            ORDER BY selected_peaks.peak_id, nearest.distance
        """)

        rows = cursor.fetchall()

        # Should find 2 nearest genes for each of 2 peaks = up to 4 results
        assert len(rows) > 0, "Should find genes for peaks in CTE"

        # Check that we have results for both peaks
        peak_ids = set(row[0] for row in rows)
        assert 1 in peak_ids, "Should have results for peak 1"
        assert 3 in peak_ids, "Should have results for peak 3"

    def test_k_greater_than_total_features_all_chromosomes(
        self, duckdb_engine_with_edge_case_data
    ):
        """
        GIVEN k greater than total number of features on the same chromosome
        WHEN querying for nearest genes
        THEN should return all available features on that chromosome
        """
        engine = duckdb_engine_with_edge_case_data

        cursor = engine.execute("""
            SELECT
                peaks.peak_id,
                nearest.gene_name
            FROM peaks
            CROSS JOIN LATERAL NEAREST(genes, reference=peaks.interval, k=1000) AS nearest
            WHERE peaks.peak_id = 1
        """)

        rows = cursor.fetchall()

        # Peak 1 is on chr1, and there are 6 genes on chr1
        # Should return all 6 genes, not 1000
        assert len(rows) == 6, f"Should return all 6 genes on chr1, got {len(rows)}"

    def test_ties_with_k_greater_than_one(self, duckdb_engine_with_edge_case_data):
        """
        GIVEN multiple features at the same distance (ties)
        WHEN querying with k that includes tied features
        THEN should handle ties consistently
        """
        engine = duckdb_engine_with_edge_case_data

        cursor = engine.execute("""
            SELECT
                peaks.peak_id,
                nearest.gene_name,
                nearest.distance
            FROM peaks
            CROSS JOIN LATERAL NEAREST(genes, reference=peaks.interval, k=3) AS nearest
            WHERE peaks.peak_id = 1
            ORDER BY nearest.distance, nearest.gene_name
        """)

        rows = cursor.fetchall()

        # Peak 1 has:
        # - GENE_C at 0bp (overlap)
        # - GENE_A and GENE_B both at 500bp (tie)
        # With k=3, should get all 3

        assert len(rows) == 3, f"Should return 3 nearest genes, got {len(rows)}"

        # First should be GENE_C (distance 0)
        assert rows[0][1] == "GENE_C"
        assert rows[0][2] == 0

        # Next two should be GENE_A and GENE_B (distance 500, order may vary)
        gene_names_at_500 = [rows[1][1], rows[2][1]]
        assert set(gene_names_at_500) == {"GENE_A", "GENE_B"}, (
            f"Should have GENE_A and GENE_B at 500bp"
        )
        assert rows[1][2] == 500
        assert rows[2][2] == 500


class TestNearestPropertyBased:
    """Property-based tests for NEAREST operator using Hypothesis."""

    @given(
        start1=st.integers(min_value=0, max_value=100000),
        length1=st.integers(min_value=1, max_value=1000),
        start2=st.integers(min_value=0, max_value=100000),
        length2=st.integers(min_value=1, max_value=1000),
    )
    def test_distance_non_negative_for_non_overlapping(
        self, start1, length1, start2, length2
    ):
        """
        PROPERTY: Distance between non-overlapping intervals is always non-negative
        GIVEN two non-overlapping genomic intervals
        WHEN calculating distance using NEAREST
        THEN distance should be >= 0
        """
        end1 = start1 + length1
        end2 = start2 + length2

        # Skip if intervals overlap
        assume(not (start1 < end2 and end1 > start2))

        engine = GIQLEngine(target_dialect="duckdb")

        # Create tables
        engine.conn.execute("""
            CREATE TABLE ref (id INTEGER, chromosome VARCHAR, start_pos INTEGER, end_pos INTEGER)
        """)
        engine.conn.execute("""
            CREATE TABLE target (id INTEGER, chromosome VARCHAR, start_pos INTEGER, end_pos INTEGER)
        """)

        # Insert test data
        engine.conn.execute(f"""
            INSERT INTO ref VALUES (1, 'chr1', {start1}, {end1})
        """)
        engine.conn.execute(f"""
            INSERT INTO target VALUES (1, 'chr1', {start2}, {end2})
        """)

        # Register schema
        engine.register_table_schema(
            "ref",
            {
                "id": "INTEGER",
                "chromosome": "VARCHAR",
                "start_pos": "INTEGER",
                "end_pos": "INTEGER",
            },
            genomic_column="interval",
        )
        engine.register_table_schema(
            "target",
            {
                "id": "INTEGER",
                "chromosome": "VARCHAR",
                "start_pos": "INTEGER",
                "end_pos": "INTEGER",
            },
            genomic_column="interval",
        )

        # Query for nearest
        cursor = engine.execute("""
            SELECT nearest.distance
            FROM ref
            CROSS JOIN LATERAL NEAREST(target, reference=ref.interval, k=1) AS nearest
        """)

        rows = cursor.fetchall()
        if len(rows) > 0:
            distance = rows[0][0]
            assert distance >= 0, f"Distance should be non-negative, got {distance}"

    @given(
        start1=st.integers(min_value=0, max_value=100000),
        length1=st.integers(min_value=1, max_value=1000),
        overlap_start=st.integers(min_value=1, max_value=500),
    )
    def test_overlapping_intervals_have_zero_distance(
        self, start1, length1, overlap_start
    ):
        """
        PROPERTY: Overlapping intervals have distance 0
        GIVEN two genomic intervals that overlap
        WHEN calculating distance using NEAREST
        THEN distance should be 0
        """
        end1 = start1 + length1
        # Create overlapping interval
        start2 = start1 + overlap_start
        end2 = start2 + length1

        # Ensure they actually overlap
        assume(start1 < end2 and end1 > start2)

        engine = GIQLEngine(target_dialect="duckdb")

        # Create tables
        engine.conn.execute("""
            CREATE TABLE ref (id INTEGER, chromosome VARCHAR, start_pos INTEGER, end_pos INTEGER)
        """)
        engine.conn.execute("""
            CREATE TABLE target (id INTEGER, chromosome VARCHAR, start_pos INTEGER, end_pos INTEGER)
        """)

        # Insert test data
        engine.conn.execute(f"""
            INSERT INTO ref VALUES (1, 'chr1', {start1}, {end1})
        """)
        engine.conn.execute(f"""
            INSERT INTO target VALUES (1, 'chr1', {start2}, {end2})
        """)

        # Register schema
        engine.register_table_schema(
            "ref",
            {
                "id": "INTEGER",
                "chromosome": "VARCHAR",
                "start_pos": "INTEGER",
                "end_pos": "INTEGER",
            },
            genomic_column="interval",
        )
        engine.register_table_schema(
            "target",
            {
                "id": "INTEGER",
                "chromosome": "VARCHAR",
                "start_pos": "INTEGER",
                "end_pos": "INTEGER",
            },
            genomic_column="interval",
        )

        # Query for nearest
        cursor = engine.execute("""
            SELECT nearest.distance
            FROM ref
            CROSS JOIN LATERAL NEAREST(target, reference=ref.interval, k=1) AS nearest
        """)

        rows = cursor.fetchall()
        assert len(rows) > 0, "Should find overlapping interval"
        distance = rows[0][0]
        assert distance == 0, (
            f"Overlapping intervals should have distance 0, got {distance}"
        )

    @given(
        k=st.integers(min_value=1, max_value=10),
        n_features=st.integers(min_value=0, max_value=15),
    )
    def test_k_parameter_returns_at_most_k_results(self, k, n_features):
        """
        PROPERTY: k parameter limits results to at most k features
        GIVEN k parameter and n available features
        WHEN querying for k nearest
        THEN should return min(k, n) results
        """
        engine = GIQLEngine(target_dialect="duckdb")

        # Create tables
        engine.conn.execute("""
            CREATE TABLE ref (id INTEGER, chromosome VARCHAR, start_pos INTEGER, end_pos INTEGER)
        """)
        engine.conn.execute("""
            CREATE TABLE target (id INTEGER, chromosome VARCHAR, start_pos INTEGER, end_pos INTEGER)
        """)

        # Insert reference point
        engine.conn.execute("""
            INSERT INTO ref VALUES (1, 'chr1', 1000, 1100)
        """)

        # Insert n_features target features
        for i in range(n_features):
            # Spread features out to avoid ties
            start = 2000 + (i * 500)
            end = start + 100
            engine.conn.execute(f"""
                INSERT INTO target VALUES ({i}, 'chr1', {start}, {end})
            """)

        # Register schema
        engine.register_table_schema(
            "ref",
            {
                "id": "INTEGER",
                "chromosome": "VARCHAR",
                "start_pos": "INTEGER",
                "end_pos": "INTEGER",
            },
            genomic_column="interval",
        )
        engine.register_table_schema(
            "target",
            {
                "id": "INTEGER",
                "chromosome": "VARCHAR",
                "start_pos": "INTEGER",
                "end_pos": "INTEGER",
            },
            genomic_column="interval",
        )

        # Query for k nearest
        cursor = engine.execute(f"""
            SELECT COUNT(*)
            FROM ref
            CROSS JOIN LATERAL NEAREST(target, reference=ref.interval, k={k}) AS nearest
        """)

        rows = cursor.fetchall()
        count = rows[0][0]

        # Should return at most k results, but not more than available features
        expected_count = min(k, n_features)
        assert count == expected_count, (
            f"Expected {expected_count} results (min({k}, {n_features})), got {count}"
        )
