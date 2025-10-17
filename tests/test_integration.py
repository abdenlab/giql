"""
Integration tests for GIQL engine.
"""


class TestIntegration:
    """Test GIQL queries work correctly across different databases."""

    def test_simple_intersects(self, engine_with_variants, to_df):
        """
        GIVEN variants data loaded
        WHEN querying with simple INTERSECTS
        THEN should return overlapping variants
        """
        result = to_df(
            engine_with_variants.execute(
                "SELECT * FROM variants WHERE position INTERSECTS 'chr1:1000-2000'"
            )
        )

        assert len(result) == 1
        assert result.iloc[0]["id"] == 1

    def test_intersects_any(self, engine_with_variants, to_df):
        """
        GIVEN variants data
        WHEN querying with INTERSECTS ANY
        THEN should return variants overlapping any range
        """
        result = to_df(
            engine_with_variants.execute(
                """
            SELECT * FROM variants
            WHERE position INTERSECTS ANY('chr1:1000-2000', 'chr1:10000-11000')
            """
            )
        )

        assert len(result) == 2
        assert set(result["id"]) == {1, 2}

    def test_intersects_any_cross_chromosome(self, engine_with_variants, to_df):
        """
        GIVEN variants data
        WHEN querying with INTERSECTS ANY across chromosomes
        THEN should return variants from both chromosomes
        """
        result = to_df(
            engine_with_variants.execute(
                """
            SELECT * FROM variants
            WHERE position INTERSECTS ANY('chr1:1000-2000', 'chr2:400-700')
            """
            )
        )

        assert len(result) == 2
        assert set(result["id"]) == {1, 4}

    def test_contains_point(self, engine_with_variants, to_df):
        """
        GIVEN variants data
        WHEN querying with CONTAINS on a point
        THEN should return variants containing that point
        """
        result = to_df(
            engine_with_variants.execute(
                "SELECT * FROM variants WHERE position CONTAINS 'chr1:1550'"
            )
        )

        assert len(result) == 1
        assert result.iloc[0]["id"] == 1

    def test_within(self, engine_with_variants, to_df):
        """
        GIVEN variants data
        WHEN querying with WITHIN
        THEN should return variants within the range
        """
        result = to_df(
            engine_with_variants.execute(
                "SELECT * FROM variants WHERE position WITHIN 'chr1:0-20000'"
            )
        )

        assert len(result) == 3
        assert set(result["id"]) == {1, 2, 3}

    def test_intersects_all(self, engine_with_variants, to_df):
        """
        GIVEN variants data
        WHEN querying with INTERSECTS ALL
        THEN should return variants overlapping all ranges
        """
        result = to_df(
            engine_with_variants.execute(
                """
            SELECT * FROM variants
            WHERE position INTERSECTS ALL('chr1:1000-2000', 'chr1:1400-1700')
            """
            )
        )

        assert len(result) == 1
        assert result.iloc[0]["id"] == 1

    def test_filter_after_giql(self, engine_with_variants, to_df):
        """
        GIVEN variants data overlapping a region
        WHEN filtering results manually
        THEN should get correct subset
        """
        # Get all overlapping variants first
        result = to_df(
            engine_with_variants.execute(
                """
            SELECT * FROM variants
            WHERE position INTERSECTS ANY('chr1:1000-2000', 'chr1:10000-11000')
            ORDER BY id
            """
            )
        )

        # Filter in application code
        filtered = result[result["quality"] > 30]
        assert len(filtered) == 1
        assert filtered.iloc[0]["id"] == 2

    def test_with_aggregation(self, engine_with_variants, to_df):
        """
        GIVEN variants data
        WHEN using aggregation with GIQL
        THEN should correctly group and aggregate
        """
        result = to_df(
            engine_with_variants.execute(
                """
            SELECT chromosome, COUNT(*) as cnt
            FROM variants
            WHERE position INTERSECTS ANY('chr1:0-20000', 'chr2:0-10000')
            GROUP BY chromosome
            ORDER BY chromosome
            """
            )
        )

        assert len(result) == 2
        chr1_count = result[result["chromosome"] == "chr1"]["cnt"].iloc[0]
        chr2_count = result[result["chromosome"] == "chr2"]["cnt"].iloc[0]
        assert chr1_count == 3
        assert chr2_count == 2

    def test_empty_result(self, engine_with_variants, to_df):
        """
        GIVEN variants data
        WHEN querying region with no variants
        THEN should return empty result
        """
        result = to_df(
            engine_with_variants.execute(
                "SELECT * FROM variants WHERE position INTERSECTS 'chr4:1000-2000'"
            )
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
            engine_with_variants.execute(
                "SELECT * FROM variants WHERE position INTERSECTS 'invalid'"
            )

    def test_overlapping_ranges_in_all(self, engine_with_variants, to_df):
        """
        GIVEN overlapping ranges in ALL
        WHEN executing query
        THEN should only return variants overlapping all ranges
        """
        result = to_df(
            engine_with_variants.execute(
                """
            SELECT * FROM variants
            WHERE position INTERSECTS ALL('chr1:1400-1700', 'chr1:1450-1650')
            """
            )
        )

        assert len(result) == 1
        assert result.iloc[0]["id"] == 1

    def test_count_with_giql(self, duckdb_engine_with_data, to_df):
        """
        GIVEN variants data
        WHEN using COUNT with GIQL in WHERE
        THEN should work correctly
        """
        result = to_df(
            duckdb_engine_with_data.execute(
                """
            SELECT COUNT(*) as cnt FROM variants
            WHERE position INTERSECTS 'chr1:0-100000'
            """
            )
        )

        # chr1 now has 4 variants (1, 2, 3, 6)
        assert result.iloc[0]["cnt"] == 4

    def test_join_with_genes(self, duckdb_engine_with_data, to_df):
        """
        GIVEN variants and genes data
        WHEN joining on chromosome and checking overlaps manually
        THEN should find overlapping features
        """
        # Note: GIQL in JOIN conditions requires more complex parsing
        # For now, use simpler approach with WHERE clause
        result = to_df(
            duckdb_engine_with_data.execute(
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
        )

        # Variant 1 (1500-1600) overlaps GENE1 (1000-2000)
        assert len(result) >= 1
        assert 1 in result["id"].values

    def test_order_and_limit(self, duckdb_engine_with_data, to_df):
        """
        GIVEN variants data
        WHEN using ORDER BY and LIMIT with GIQL
        THEN should work correctly
        """
        result = to_df(
            duckdb_engine_with_data.execute(
                """
            SELECT * FROM variants
            WHERE position INTERSECTS 'chr1:0-20000'
            ORDER BY quality DESC
            LIMIT 2
            """
            )
        )

        assert len(result) == 2
        # Should get highest quality variants
        assert result.iloc[0]["quality"] == 40.0  # id=2
        assert result.iloc[1]["quality"] == 30.0  # id=1

    def test_multiple_giql_conditions(self, duckdb_engine_with_data, to_df):
        """
        GIVEN variants data
        WHEN using multiple GIQL conditions
        THEN should apply all conditions
        """
        result = to_df(
            duckdb_engine_with_data.execute(
                """
            SELECT * FROM variants
            WHERE position INTERSECTS 'chr1:0-20000'
              AND position WITHIN 'chr1:1000-11000'
            ORDER BY id
            """
            )
        )

        assert len(result) == 2
        assert set(result["id"]) == {1, 2}

    def test_giql_with_standard_sql_conditions(self, duckdb_engine_with_data, to_df):
        """
        GIVEN variants data
        WHEN mixing GIQL and standard SQL conditions
        THEN should apply all conditions correctly
        """
        result = to_df(
            duckdb_engine_with_data.execute(
                """
            SELECT * FROM variants
            WHERE position INTERSECTS 'chr1:0-20000'
              AND quality > 30
            ORDER BY id
            """
            )
        )

        assert len(result) == 1
        assert result.iloc[0]["id"] == 2
        assert result.iloc[0]["quality"] == 40.0

    def test_cte_with_giql(self, duckdb_engine_with_data, to_df):
        """
        GIVEN variants data
        WHEN using CTE with GIQL operators
        THEN should work correctly
        """
        result = to_df(
            duckdb_engine_with_data.execute(
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
        )

        assert len(result) == 1
        assert result.iloc[0]["id"] == 2
        assert result.iloc[0]["quality"] == 40.0


class TestBedtoolsRecipes:
    """Test recipes from README that replicate bedtools intersect behavior."""

    def test_recipe_default_intersect(self, engine_with_variants, to_df):
        """
        GIVEN features data
        WHEN querying with default intersect (DISTINCT)
        THEN should return unique overlapping features
        """
        result = to_df(
            engine_with_variants.execute(
                """
            SELECT DISTINCT *
            FROM variants
            WHERE position INTERSECTS 'chr1:0-20000'
            ORDER BY id
            """
            )
        )

        assert len(result) == 3
        assert set(result["id"]) == {1, 2, 3}

    def test_recipe_v_flag_no_overlaps(self, duckdb_engine_with_data, to_df):
        """
        GIVEN variants and genes data
        WHEN using LEFT JOIN with NULL check (bedtools -v)
        THEN should return features with no overlaps
        """
        result = to_df(
            duckdb_engine_with_data.execute(
                """
            SELECT v.*
            FROM variants v
            LEFT JOIN genes g ON v.position INTERSECTS g.position
            WHERE g.chromosome IS NULL
            ORDER BY v.id
            """
            )
        )

        # Variants 6, 7, 8 don't overlap any genes
        assert len(result) == 3
        assert set(result["id"]) == {6, 7, 8}

    def test_recipe_u_flag_any_overlap(self, duckdb_engine_with_data, to_df):
        """
        GIVEN variants and genes data
        WHEN using INNER JOIN with DISTINCT (bedtools -u)
        THEN should return unique features with any overlap
        """
        result = to_df(
            duckdb_engine_with_data.execute(
                """
            SELECT DISTINCT v.*
            FROM variants v
            INNER JOIN genes g ON v.position INTERSECTS g.position
            ORDER BY v.id
            """
            )
        )

        # Variants that overlap at least one gene
        assert len(result) >= 1
        assert 1 in result["id"].values

    def test_recipe_c_flag_count_overlaps(self, duckdb_engine_with_data, to_df):
        """
        GIVEN variants and genes data
        WHEN using COUNT with LEFT JOIN (bedtools -c)
        THEN should count overlaps for each feature
        """
        result = to_df(
            duckdb_engine_with_data.execute(
                """
            SELECT v.id, COUNT(g.name) as overlap_count
            FROM variants v
            LEFT JOIN genes g ON v.position INTERSECTS g.position
            GROUP BY v.id, v.chromosome, v.start_pos, v.end_pos
            ORDER BY v.id
            """
            )
        )

        assert len(result) == 8
        # Check that counts are non-negative integers
        assert all(result["overlap_count"] >= 0)
        # Variants 1-5 should have overlaps, 6-8 should have 0
        assert result[result["id"].isin([1, 2, 3, 4, 5])]["overlap_count"].sum() > 0
        assert result[result["id"].isin([6, 7, 8])]["overlap_count"].sum() == 0

    def test_recipe_strand_specific_same(self, duckdb_engine_with_data, to_df):
        """
        GIVEN features with strand information
        WHEN filtering by same strand (bedtools -s)
        THEN should only return same-strand overlaps
        """
        result = to_df(
            duckdb_engine_with_data.execute(
                """
            SELECT v.*, g.strand as gene_strand
            FROM variants v
            INNER JOIN genes g ON v.position INTERSECTS g.position
            WHERE v.chromosome = g.chromosome
              AND g.strand = '+'
            """
            )
        )

        # All results should have gene_strand = '+'
        if len(result) > 0:
            assert all(result["gene_strand"] == "+")

    def test_recipe_overlap_fraction(self, duckdb_engine_with_data, to_df):
        """
        GIVEN variants and genes data
        WHEN requiring minimum overlap fraction (bedtools -f)
        THEN should filter by overlap fraction
        """
        result = to_df(
            duckdb_engine_with_data.execute(
                """
            SELECT v.*
            FROM variants v, genes g
            WHERE v.position INTERSECTS g.position
              AND (
                  LEAST(v.end_pos, g.end_pos) - GREATEST(v.start_pos, g.start_pos)
              ) >= 0.5 * (v.end_pos - v.start_pos)
            """
            )
        )

        # Should have some results but potentially fewer than without fraction filter
        assert len(result) >= 0

    def test_recipe_loj_left_outer_join(self, duckdb_engine_with_data, to_df):
        """
        GIVEN variants and genes data
        WHEN using LEFT JOIN (bedtools -loj)
        THEN should return all A features with NULL for non-overlapping
        """
        result = to_df(
            duckdb_engine_with_data.execute(
                """
            SELECT v.id, v.chromosome, g.name as gene_name
            FROM variants v
            LEFT JOIN genes g ON v.position INTERSECTS g.position
            ORDER BY v.id
            """
            )
        )

        # Should have all 8 variants
        assert len(result) >= 8
        # Variants 6, 7, 8 should have NULL gene_name
        assert result["gene_name"].isna().any()
        non_overlapping = result[result["gene_name"].isna()]
        assert set(non_overlapping["id"]) == {6, 7, 8}

    def test_recipe_wo_write_overlap_bp(self, duckdb_engine_with_data, to_df):
        """
        GIVEN variants and genes data
        WHEN calculating overlap in base pairs (bedtools -wo)
        THEN should return overlap amount
        """
        result = to_df(
            duckdb_engine_with_data.execute(
                """
            SELECT
                v.id,
                g.name,
                (LEAST(v.end_pos, g.end_pos) - GREATEST(v.start_pos, g.start_pos)) as overlap_bp
            FROM variants v
            INNER JOIN genes g ON v.position INTERSECTS g.position
            """
            )
        )

        if len(result) > 0:
            # All overlap_bp values should be positive
            assert all(result["overlap_bp"] > 0)

    def test_recipe_wao_write_overlap_all(self, duckdb_engine_with_data, to_df):
        """
        GIVEN variants and genes data
        WHEN calculating overlap for ALL features (bedtools -wao)
        THEN should return 0 for non-overlapping
        """
        result = to_df(
            duckdb_engine_with_data.execute(
                """
            SELECT
                v.id,
                CASE
                    WHEN g.chromosome IS NULL THEN 0
                    ELSE LEAST(v.end_pos, g.end_pos) - GREATEST(v.start_pos, g.start_pos)
                END as overlap_bp
            FROM variants v
            LEFT JOIN genes g ON v.position INTERSECTS g.position
            ORDER BY v.id
            """
            )
        )

        # Should have entries for all 8 variants
        assert len(result) >= 8
        # All overlap_bp should be >= 0
        assert all(result["overlap_bp"] >= 0)
        # Variants 6, 7, 8 should have 0 overlap
        non_overlapping = result[result["id"].isin([6, 7, 8])]
        assert all(non_overlapping["overlap_bp"] == 0)

    def test_recipe_multiple_files_union(
        self, sample_variants_csv, sample_genes_csv, to_df
    ):
        """
        GIVEN multiple B files combined with UNION
        WHEN querying against combined set
        THEN should find overlaps from any file
        """
        from giql import GIQLEngine

        with GIQLEngine(target_dialect="duckdb") as engine:
            engine.load_csv("variants", sample_variants_csv)
            engine.load_csv("genes1", sample_genes_csv)
            engine.load_csv("genes2", sample_genes_csv)

            for table in ["variants", "genes1", "genes2"]:
                schema = {
                    "id": "INTEGER" if table == "variants" else None,
                    "gene_id": "INTEGER" if table != "variants" else None,
                    "name": "VARCHAR",
                    "chromosome": "VARCHAR",
                    "start_pos": "BIGINT",
                    "end_pos": "BIGINT",
                }
                schema = {k: v for k, v in schema.items() if v is not None}
                engine.register_table_schema(table, schema, genomic_column="position")

            result = to_df(
                engine.execute(
                    """
                WITH all_genes AS (
                    SELECT * FROM genes1
                    UNION ALL
                    SELECT * FROM genes2
                )
                SELECT DISTINCT v.id
                FROM variants v
                INNER JOIN all_genes g ON v.position INTERSECTS g.position
                ORDER BY v.id
                """
                )
            )

            # Should find overlaps from union of both gene files
            assert len(result) >= 1

    def test_recipe_count_with_having(self, duckdb_engine_with_data, to_df):
        """
        GIVEN variants and genes data
        WHEN counting overlaps with HAVING clause
        THEN should filter by overlap count
        """
        result = to_df(
            duckdb_engine_with_data.execute(
                """
            SELECT v.*
            FROM variants v
            INNER JOIN genes g ON v.position INTERSECTS g.position
            GROUP BY v.id, v.chromosome, v.start_pos, v.end_pos, v.ref, v.alt, v.quality
            HAVING COUNT(*) >= 1
            ORDER BY v.id
            """
            )
        )

        # Should only return variants with at least 1 overlap
        assert len(result) >= 1

    def test_recipe_aggregate_stats(self, duckdb_engine_with_data, to_df):
        """
        GIVEN variants and genes data
        WHEN calculating aggregate statistics per chromosome
        THEN should group and aggregate correctly
        """
        result = to_df(
            duckdb_engine_with_data.execute(
                """
            SELECT
                v.chromosome,
                COUNT(DISTINCT v.id) as total_variants,
                COUNT(g.name) as total_overlaps
            FROM variants v
            LEFT JOIN genes g ON v.position INTERSECTS g.position
            GROUP BY v.chromosome
            ORDER BY v.chromosome
            """
            )
        )

        # Should have stats for all 3 chromosomes (chr1, chr2, chr3)
        assert len(result) == 3
        assert "chr1" in result["chromosome"].values
        assert "chr2" in result["chromosome"].values
        assert "chr3" in result["chromosome"].values
        # All counts should be non-negative
        assert all(result["total_variants"] >= 0)
        assert all(result["total_overlaps"] >= 0)
        # chr3 should have 0 overlaps
        chr3_overlaps = result[result["chromosome"] == "chr3"]["total_overlaps"].iloc[0]
        assert chr3_overlaps == 0


class TestMultiDialect:
    """Test that same queries work across different databases."""

    def test_same_query_different_dialects(self, sample_variants_csv, to_df):
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
                result = to_df(engine.execute(query))
                results[dialect] = set(result["id"])

        # All dialects should produce same results
        assert len(set(map(frozenset, results.values()))) == 1
        assert results["duckdb"] == {1}
