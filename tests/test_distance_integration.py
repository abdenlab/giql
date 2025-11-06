"""Integration tests for DISTANCE operator end-to-end execution.

Tests verify that DISTANCE queries execute correctly against actual databases.
"""

import duckdb
import pytest

from giql import GIQLEngine


class TestDistanceIntegration:
    """Integration tests for DISTANCE query execution."""

    def test_basic_distance_query_execution(self):
        """
        GIVEN two tables with genomic features
        WHEN executing a query with DISTANCE()
        THEN should return correct distance values
        """
        # Create engine
        engine = GIQLEngine(target_dialect="duckdb")

        # Create test data using engine's connection
        # Table A: peaks
        engine.conn.execute("""
            CREATE TABLE peaks (
                name VARCHAR,
                chromosome VARCHAR,
                start_pos INTEGER,
                end_pos INTEGER
            )
        """)
        engine.conn.execute("""
            INSERT INTO peaks VALUES
                ('peak1', 'chr1', 100, 200),
                ('peak2', 'chr1', 500, 600),
                ('peak3', 'chr2', 1000, 1100)
        """)

        # Table B: genes
        engine.conn.execute("""
            CREATE TABLE genes (
                name VARCHAR,
                chromosome VARCHAR,
                start_pos INTEGER,
                end_pos INTEGER
            )
        """)
        engine.conn.execute("""
            INSERT INTO genes VALUES
                ('geneA', 'chr1', 150, 300),
                ('geneB', 'chr1', 700, 800),
                ('geneC', 'chr2', 2000, 2100)
        """)

        # Register schemas
        engine.register_table_schema(
            "peaks",
            {
                "name": "VARCHAR",
                "chromosome": "VARCHAR",
                "start_pos": "INTEGER",
                "end_pos": "INTEGER",
            },
            genomic_column="position",
        )
        engine.register_table_schema(
            "genes",
            {
                "name": "VARCHAR",
                "chromosome": "VARCHAR",
                "start_pos": "INTEGER",
                "end_pos": "INTEGER",
            },
            genomic_column="position",
        )

        # Execute DISTANCE query
        query = """
        SELECT
            p.name as peak,
            g.name as gene,
            DISTANCE(p.position, g.position) as distance
        FROM peaks p
        CROSS JOIN genes g
        ORDER BY p.name, g.name
        """

        cursor = engine.execute(query)
        results = cursor.fetchall()

        # Expected results:
        # peak1 (chr1:100-200) vs geneA (chr1:150-300) = 0 (overlap)
        # peak1 (chr1:100-200) vs geneB (chr1:700-800) = 500 (gap: 700-200)
        # peak1 (chr1:100-200) vs geneC (chr2:2000-2100) = NULL (different chrom)
        # peak2 (chr1:500-600) vs geneA (chr1:150-300) = 200 (gap: 500-300)
        # peak2 (chr1:500-600) vs geneB (chr1:700-800) = 100 (gap: 700-600)
        # peak2 (chr1:500-600) vs geneC (chr2:2000-2100) = NULL (different chrom)
        # peak3 (chr2:1000-1100) vs geneA (chr1:150-300) = NULL (different chrom)
        # peak3 (chr2:1000-1100) vs geneB (chr1:700-800) = NULL (different chrom)
        # peak3 (chr2:1000-1100) vs geneC (chr2:2000-2100) = 900 (gap: 2000-1100)

        expected = [
            ('peak1', 'geneA', 0),
            ('peak1', 'geneB', 500),
            ('peak1', 'geneC', None),
            ('peak2', 'geneA', 200),
            ('peak2', 'geneB', 100),
            ('peak2', 'geneC', None),
            ('peak3', 'geneA', None),
            ('peak3', 'geneB', None),
            ('peak3', 'geneC', 900),
        ]

        assert results == expected, f"Expected:\n{expected}\n\nGot:\n{results}"

        engine.close()
