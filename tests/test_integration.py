"""
Integration tests for GIQL engine.
"""


class TestIntegration:
    """Test GIQL queries work correctly across different databases."""

    def test_simple_intersects(self, engine_with_variants):
        """
        GIVEN variants data loaded
        WHEN querying with simple INTERSECTS
        THEN should return overlapping variants
        """
        result = engine_with_variants.query(
            "SELECT * FROM variants WHERE position INTERSECTS 'chr1:1000-2000'"
        )

        assert len(result) == 1
        assert result.iloc[0]["id"] == 1

    def test_intersects_any(self, engine_with_variants):
        """
        GIVEN variants data
        WHEN querying with INTERSECTS ANY
        THEN should return variants overlapping any range
        """
        result = engine_with_variants.query(
            """
            SELECT * FROM variants
            WHERE position INTERSECTS ANY('chr1:1000-2000', 'chr1:10000-11000')
            """
        )

        assert len(result) == 2
        assert set(result["id"]) == {1, 2}

    def test_intersects_any_cross_chromosome(self, engine_with_variants):
        """
        GIVEN variants data
        WHEN querying with INTERSECTS ANY across chromosomes
        THEN should return variants from both chromosomes
        """
        result = engine_with_variants.query(
            """
            SELECT * FROM variants
            WHERE position INTERSECTS ANY('chr1:1000-2000', 'chr2:400-700')
            """
        )

        assert len(result) == 2
        assert set(result["id"]) == {1, 4}

    def test_contains_point(self, engine_with_variants):
        """
        GIVEN variants data
        WHEN querying with CONTAINS on a point
        THEN should return variants containing that point
        """
        result = engine_with_variants.query(
            "SELECT * FROM variants WHERE position CONTAINS 'chr1:1550'"
        )

        assert len(result) == 1
        assert result.iloc[0]["id"] == 1

    def test_within(self, engine_with_variants):
        """
        GIVEN variants data
        WHEN querying with WITHIN
        THEN should return variants within the range
        """
        result = engine_with_variants.query(
            "SELECT * FROM variants WHERE position WITHIN 'chr1:0-20000'"
        )

        assert len(result) == 3
        assert set(result["id"]) == {1, 2, 3}

    def test_intersects_all(self, engine_with_variants):
        """
        GIVEN variants data
        WHEN querying with INTERSECTS ALL
        THEN should return variants overlapping all ranges
        """
        result = engine_with_variants.query(
            """
            SELECT * FROM variants
            WHERE position INTERSECTS ALL('chr1:1000-2000', 'chr1:1400-1700')
            """
        )

        assert len(result) == 1
        assert result.iloc[0]["id"] == 1

    def test_filter_after_giql(self, engine_with_variants):
        """
        GIVEN variants data overlapping a region
        WHEN filtering results manually
        THEN should get correct subset
        """
        # Get all overlapping variants first
        result = engine_with_variants.query(
            """
            SELECT * FROM variants
            WHERE position INTERSECTS ANY('chr1:1000-2000', 'chr1:10000-11000')
            ORDER BY id
            """
        )

        # Filter in application code
        filtered = result[result["quality"] > 30]
        assert len(filtered) == 1
        assert filtered.iloc[0]["id"] == 2

    def test_with_aggregation(self, engine_with_variants):
        """
        GIVEN variants data
        WHEN using aggregation with GIQL
        THEN should correctly group and aggregate
        """
        result = engine_with_variants.query(
            """
            SELECT chromosome, COUNT(*) as cnt
            FROM variants
            WHERE position INTERSECTS ANY('chr1:0-20000', 'chr2:0-10000')
            GROUP BY chromosome
            ORDER BY chromosome
            """
        )

        assert len(result) == 2
        chr1_count = result[result["chromosome"] == "chr1"]["cnt"].iloc[0]
        chr2_count = result[result["chromosome"] == "chr2"]["cnt"].iloc[0]
        assert chr1_count == 3
        assert chr2_count == 2

    def test_empty_result(self, engine_with_variants):
        """
        GIVEN variants data
        WHEN querying region with no variants
        THEN should return empty result
        """
        result = engine_with_variants.query(
            "SELECT * FROM variants WHERE position INTERSECTS 'chr3:1000-2000'"
        )

        assert len(result) == 0

    def test_invalid_range_format(self, engine_with_variants):
        """
        GIVEN invalid range format
        WHEN executing query
        THEN should raise ValueError
        """
        import pytest

        with pytest.raises(ValueError, match="Could not parse genomic range"):
            engine_with_variants.query(
                "SELECT * FROM variants WHERE position INTERSECTS 'invalid'"
            )

    def test_overlapping_ranges_in_all(self, engine_with_variants):
        """
        GIVEN overlapping ranges in ALL
        WHEN executing query
        THEN should only return variants overlapping all ranges
        """
        result = engine_with_variants.query(
            """
            SELECT * FROM variants
            WHERE position INTERSECTS ALL('chr1:1400-1700', 'chr1:1450-1650')
            """
        )

        assert len(result) == 1
        assert result.iloc[0]["id"] == 1


class TestComplexQueries:
    """Test complex SQL features with GIQL."""

    def test_count_with_giql(self, duckdb_engine_with_data):
        """
        GIVEN variants data
        WHEN using COUNT with GIQL in WHERE
        THEN should work correctly
        """
        result = duckdb_engine_with_data.query(
            """
            SELECT COUNT(*) as cnt FROM variants
            WHERE position INTERSECTS 'chr1:0-100000'
            """
        )

        assert result.iloc[0]["cnt"] == 3

    def test_join_with_genes(self, duckdb_engine_with_data):
        """
        GIVEN variants and genes data
        WHEN joining on chromosome and checking overlaps manually
        THEN should find overlapping features
        """
        # Note: GIQL in JOIN conditions requires more complex parsing
        # For now, use simpler approach with WHERE clause
        result = duckdb_engine_with_data.query(
            """
            SELECT v.id, v.ref, v.alt, g.name
            FROM variants v, genes g
            WHERE v.chromosome = g.chromosome
              AND v.chromosome = 'chr1'
              AND v.position INTERSECTS 'chr1:1000-2000'
              AND g.name = 'GENE1'
            ORDER BY v.id
            """
        )

        # Variant 1 (1500-1600) overlaps GENE1 (1000-2000)
        assert len(result) >= 1
        assert 1 in result["id"].values

    def test_order_and_limit(self, duckdb_engine_with_data):
        """
        GIVEN variants data
        WHEN using ORDER BY and LIMIT with GIQL
        THEN should work correctly
        """
        result = duckdb_engine_with_data.query(
            """
            SELECT * FROM variants
            WHERE position INTERSECTS 'chr1:0-20000'
            ORDER BY quality DESC
            LIMIT 2
            """
        )

        assert len(result) == 2
        # Should get highest quality variants
        assert result.iloc[0]["quality"] == 40.0  # id=2
        assert result.iloc[1]["quality"] == 30.0  # id=1

    def test_multiple_giql_conditions(self, duckdb_engine_with_data):
        """
        GIVEN variants data
        WHEN using multiple GIQL conditions
        THEN should apply all conditions
        """
        result = duckdb_engine_with_data.query(
            """
            SELECT * FROM variants
            WHERE position INTERSECTS 'chr1:0-20000'
              AND position WITHIN 'chr1:1000-11000'
            ORDER BY id
            """
        )

        assert len(result) == 2
        assert set(result["id"]) == {1, 2}

    def test_giql_with_standard_sql_conditions(self, duckdb_engine_with_data):
        """
        GIVEN variants data
        WHEN mixing GIQL and standard SQL conditions
        THEN should apply all conditions correctly
        """
        result = duckdb_engine_with_data.query(
            """
            SELECT * FROM variants
            WHERE position INTERSECTS 'chr1:0-20000'
              AND quality > 30
            ORDER BY id
            """
        )

        assert len(result) == 1
        assert result.iloc[0]["id"] == 2
        assert result.iloc[0]["quality"] == 40.0

    def test_cte_with_giql(self, duckdb_engine_with_data):
        """
        GIVEN variants data
        WHEN using CTE with GIQL operators
        THEN should work correctly
        """
        result = duckdb_engine_with_data.query(
            """
            WITH overlapping AS (
                SELECT * FROM variants
                WHERE position INTERSECTS 'chr1:0-20000'
            )
            SELECT * FROM overlapping
            WHERE quality > 30
            ORDER BY id
            """
        )

        assert len(result) == 1
        assert result.iloc[0]["id"] == 2
        assert result.iloc[0]["quality"] == 40.0


class TestMultiDialect:
    """Test that same queries work across different databases."""

    def test_same_query_different_dialects(self, sample_variants_csv):
        """
        GIVEN the same query
        WHEN executed on different dialects
        THEN should produce same results
        """
        query = """
            SELECT id FROM variants
            WHERE position INTERSECTS 'chr1:1000-2000'
            ORDER BY id
        """

        results = {}
        for dialect in ["duckdb", "sqlite"]:
            from giql import GIQLEngine

            with GIQLEngine(target_dialect=dialect) as engine:
                engine.load_csv("variants", sample_variants_csv)
                engine.register_table_schema(
                    "variants",
                    {
                        "id": "INTEGER",
                        "chromosome": "VARCHAR",
                        "start_pos": "BIGINT",
                        "end_pos": "BIGINT",
                    },
                    genomic_column="position",
                )
                result = engine.query(query)
                results[dialect] = set(result["id"])

        # All dialects should produce same results
        assert len(set(map(frozenset, results.values()))) == 1
        assert results["duckdb"] == {1}
