import tempfile

from hypothesis import given
from hypothesis import settings
from hypothesis import strategies as st

from giql import GIQLEngine


class TestGIQLEngine:
    def test_engine_initialization_duckdb(self):
        """
        GIVEN GIQLEngine with duckdb dialect
        WHEN initializing engine
        THEN should create connection successfully
        """
        engine = GIQLEngine(target_dialect="duckdb")
        assert engine.target_dialect == "duckdb"
        assert engine.conn is not None
        engine.close()

    def test_engine_initialization_sqlite(self):
        """
        GIVEN GIQLEngine with sqlite dialect
        WHEN initializing engine
        THEN should create connection successfully
        """
        with tempfile.NamedTemporaryFile(suffix=".db", delete=False) as tmp:
            engine = GIQLEngine(target_dialect="sqlite", db_path=tmp.name)
            assert engine.target_dialect == "sqlite"
            assert engine.conn is not None
            engine.close()

    def test_engine_context_manager(self):
        """
        GIVEN GIQLEngine used as context manager
        WHEN exiting context
        THEN should close connection automatically
        """
        with GIQLEngine() as engine:
            assert engine.conn is not None

    def test_load_csv_and_query_duckdb(self, tmp_path):
        """
        GIVEN CSV data loaded into DuckDB
        WHEN executing GIQL query
        THEN should return correct results
        """
        # Create sample CSV
        csv_content = """id,chromosome,start_pos,end_pos,ref,alt
1,chr1,1500,1600,A,T
2,chr1,10500,10600,G,C
3,chr2,500,600,C,G
"""
        csv_path = tmp_path / "variants.csv"
        csv_path.write_text(csv_content)

        with GIQLEngine(target_dialect="duckdb") as engine:
            engine.load_csv("variants", str(csv_path))

            # Query using INTERSECTS
            result = engine.query(
                "SELECT * FROM variants WHERE position INTERSECTS 'chr1:1000-2000'"
            )

            assert len(result) == 1
            assert result.iloc[0]["id"] == 1

    def test_load_csv_and_query_sqlite(self, tmp_path):
        """
        GIVEN CSV data loaded into SQLite
        WHEN executing GIQL query
        THEN should return correct results
        """
        # Create sample CSV
        csv_content = """id,chromosome,start_pos,end_pos,ref,alt
1,chr1,1500,1600,A,T
2,chr1,10500,10600,G,C
3,chr2,500,600,C,G
"""
        csv_path = tmp_path / "variants.csv"
        csv_path.write_text(csv_content)

        with GIQLEngine(target_dialect="sqlite") as engine:
            engine.load_csv("variants", str(csv_path))

            # Query using INTERSECTS
            result = engine.query(
                "SELECT * FROM variants WHERE position INTERSECTS 'chr1:1000-2000'"
            )

            assert len(result) == 1
            assert result.iloc[0]["id"] == 1

    def test_intersects_any_query(self, tmp_path):
        """
        GIVEN variants data
        WHEN querying with INTERSECTS ANY
        THEN should return variants overlapping any range
        """
        csv_content = """id,chromosome,start_pos,end_pos
1,chr1,1500,1600
2,chr1,10500,10600
3,chr2,500,600
"""
        csv_path = tmp_path / "variants.csv"
        csv_path.write_text(csv_content)

        with GIQLEngine(target_dialect="duckdb") as engine:
            engine.load_csv("variants", str(csv_path))

            result = engine.query(
                "SELECT * FROM variants "
                "WHERE position INTERSECTS ANY('chr1:1000-2000', 'chr2:400-700')"
            )

            assert len(result) == 2
            assert set(result["id"]) == {1, 3}

    def test_contains_query(self, tmp_path):
        """
        GIVEN variants data
        WHEN querying with CONTAINS
        THEN should return variants containing the point
        """
        csv_content = """id,chromosome,start_pos,end_pos
1,chr1,1500,1600
2,chr1,10500,10600
"""
        csv_path = tmp_path / "variants.csv"
        csv_path.write_text(csv_content)

        with GIQLEngine(target_dialect="duckdb") as engine:
            engine.load_csv("variants", str(csv_path))

            result = engine.query(
                "SELECT * FROM variants WHERE position CONTAINS 'chr1:1550'"
            )

            assert len(result) == 1
            assert result.iloc[0]["id"] == 1

    def test_within_query(self, tmp_path):
        """
        GIVEN variants data
        WHEN querying with WITHIN
        THEN should return variants within the range
        """
        csv_content = """id,chromosome,start_pos,end_pos
1,chr1,1500,1600
2,chr1,10500,10600
3,chr1,15000,15100
"""
        csv_path = tmp_path / "variants.csv"
        csv_path.write_text(csv_content)

        with GIQLEngine(target_dialect="duckdb") as engine:
            engine.load_csv("variants", str(csv_path))

            result = engine.query(
                "SELECT * FROM variants WHERE position WITHIN 'chr1:1000-11000'"
            )

            assert len(result) == 2
            assert set(result["id"]) == {1, 2}

    def test_verbose_mode(self, tmp_path):
        """
        GIVEN engine with verbose mode
        WHEN executing query
        THEN should print transpiled SQL
        """
        csv_content = """id,chromosome,start_pos,end_pos
1,chr1,1500,1600
"""
        csv_path = tmp_path / "variants.csv"
        csv_path.write_text(csv_content)

        with GIQLEngine(target_dialect="duckdb", verbose=True) as engine:
            engine.load_csv("variants", str(csv_path))
            result = engine.query(
                "SELECT * FROM variants WHERE position INTERSECTS 'chr1:1000-2000'"
            )
            assert len(result) == 1

    @given(
        chrom_col=st.sampled_from(["chromosome", "chr", "chrom", "contig", "seqname"]),
        start_col=st.sampled_from(["start_pos", "start", "begin", "pos", "chromStart"]),
        end_col=st.sampled_from(["end_pos", "end", "stop", "chromEnd"]),
        strand_col=st.sampled_from(["strand", "str", "orientation", "direction"]),
    )
    def test_custom_genomic_columns(self, chrom_col, start_col, end_col, strand_col):
        """
        GIVEN CSV data with custom genomic column names
        WHEN registering schema with custom column mappings
        THEN queries should work correctly with any valid column names
        """
        # Create temporary directory and CSV with custom column names
        with tempfile.TemporaryDirectory() as tmp_dir:
            csv_content = f"""id,{chrom_col},{start_col},{end_col},{strand_col},name
1,chr1,1500,1600,+,variant1
2,chr1,10500,10600,-,variant2
3,chr2,500,600,+,variant3
4,chr1,1400,1700,+,variant4
"""
            csv_path = f"{tmp_dir}/custom_variants.csv"
            with open(csv_path, "w") as f:
                f.write(csv_content)

            with GIQLEngine(target_dialect="duckdb", verbose=False) as engine:
                engine.load_csv("variants", csv_path)

                # Register schema with custom column names
                engine.register_table_schema(
                    "variants",
                    {
                        "id": "INTEGER",
                        chrom_col: "VARCHAR",
                        start_col: "BIGINT",
                        end_col: "BIGINT",
                        strand_col: "VARCHAR",
                        "name": "VARCHAR",
                    },
                    genomic_column="position",
                    chrom_col=chrom_col,
                    start_col=start_col,
                    end_col=end_col,
                    strand_col=strand_col,
                )

                # Test INTERSECTS query
                result = engine.query(
                    "SELECT * FROM variants WHERE position INTERSECTS 'chr1:1000-2000'"
                )
                assert len(result) == 2
                assert set(result["id"]) == {1, 4}

                # Test CLUSTER query (uses genomic columns internally)
                result = engine.query(
                    "SELECT *, CLUSTER(position) AS cluster_id FROM variants ORDER BY id"
                )
                assert len(result) == 4
                # Variants 1 and 4 should cluster together (overlapping on chr1)
                assert result.iloc[0]["cluster_id"] == result.iloc[3]["cluster_id"]
                # Variant 2 should be in different cluster (no overlap)
                assert result.iloc[1]["cluster_id"] != result.iloc[0]["cluster_id"]

                # Test stranded CLUSTER query
                result = engine.query(
                    """SELECT *, CLUSTER(position, stranded=TRUE) AS cluster_id
                       FROM variants ORDER BY id"""
                )
                assert len(result) == 4
                # With stranded=TRUE, variants 1 and 4 should cluster together (both + and overlapping)
                assert result.iloc[0]["cluster_id"] == result.iloc[3]["cluster_id"]
                # Note: cluster_ids are independent per (chromosome, strand) partition
                # So variants on different strands CAN have the same cluster_id number
                assert "cluster_id" in result.columns

                # Test MERGE query
                result = engine.query("SELECT MERGE(position) FROM variants")
                # Should merge overlapping intervals
                assert len(result) >= 1

    @given(
        # Table 1 (variants) column names
        v_chrom_col=st.sampled_from(["chromosome", "chr", "chrom"]),
        v_start_col=st.sampled_from(["start_pos", "start", "begin"]),
        v_end_col=st.sampled_from(["end_pos", "end", "stop"]),
        # Table 2 (features) column names (use different names to ensure they're distinct)
        f_chrom_col=st.sampled_from(["seqname", "contig", "chr_name"]),
        f_start_col=st.sampled_from(["pos", "chromStart", "feature_start"]),
        f_end_col=st.sampled_from(["chromEnd", "feature_end", "terminus"]),
    )
    @settings(deadline=None)
    def test_join_with_different_schemas(
        self, v_chrom_col, v_start_col, v_end_col, f_chrom_col, f_start_col, f_end_col
    ):
        """
        GIVEN two tables with different custom genomic column schemas
        WHEN joining them using INTERSECTS
        THEN queries should correctly use each table's custom column names
        """
        with tempfile.TemporaryDirectory() as tmp_dir:
            # Create variants table CSV
            variants_csv = f"""id,{v_chrom_col},{v_start_col},{v_end_col},name
1,chr1,1500,1600,var1
2,chr1,10500,10600,var2
3,chr2,500,600,var3
"""
            variants_path = f"{tmp_dir}/variants.csv"
            with open(variants_path, "w") as f:
                f.write(variants_csv)

            # Create features table CSV with DIFFERENT column names
            features_csv = f"""id,{f_chrom_col},{f_start_col},{f_end_col},type
1,chr1,1000,2000,exon
2,chr1,10000,11000,intron
3,chr2,400,700,promoter
"""
            features_path = f"{tmp_dir}/features.csv"
            with open(features_path, "w") as f:
                f.write(features_csv)

            with GIQLEngine(target_dialect="duckdb", verbose=False) as engine:
                # Load both tables
                engine.load_csv("variants", variants_path)
                engine.load_csv("features", features_path)

                # Register schemas with different column names
                engine.register_table_schema(
                    "variants",
                    {
                        "id": "INTEGER",
                        v_chrom_col: "VARCHAR",
                        v_start_col: "BIGINT",
                        v_end_col: "BIGINT",
                        "name": "VARCHAR",
                    },
                    genomic_column="position",
                    chrom_col=v_chrom_col,
                    start_col=v_start_col,
                    end_col=v_end_col,
                )

                engine.register_table_schema(
                    "features",
                    {
                        "id": "INTEGER",
                        f_chrom_col: "VARCHAR",
                        f_start_col: "BIGINT",
                        f_end_col: "BIGINT",
                        "type": "VARCHAR",
                    },
                    genomic_column="region",
                    chrom_col=f_chrom_col,
                    start_col=f_start_col,
                    end_col=f_end_col,
                )

                # Test JOIN with INTERSECTS on both tables
                result = engine.query(
                    """
                    SELECT v.name, f.type
                    FROM variants v
                    JOIN features f ON v.position INTERSECTS f.region
                    ORDER BY v.id
                    """
                )

                # Variant 1 (chr1:1500-1600) intersects Feature 1 (chr1:1000-2000)
                # Variant 2 (chr1:10500-10600) intersects Feature 2 (chr1:10000-11000)
                # Variant 3 (chr2:500-600) intersects Feature 3 (chr2:400-700)
                assert len(result) == 3
                assert list(result["name"]) == ["var1", "var2", "var3"]
                assert list(result["type"]) == ["exon", "intron", "promoter"]

                # Test LEFT JOIN to verify schema resolution works
                result = engine.query(
                    """
                    SELECT v.id, v.name, f.type
                    FROM variants v
                    LEFT JOIN features f ON v.position INTERSECTS f.region
                    WHERE v.id = 1
                    """
                )
                assert len(result) == 1
                assert result.iloc[0]["name"] == "var1"
                assert result.iloc[0]["type"] == "exon"

                # Test WHERE clause with INTERSECTS on specific table
                result = engine.query(
                    """
                    SELECT v.id, v.name
                    FROM variants v, features f
                    WHERE v.position INTERSECTS f.region
                    AND v.position INTERSECTS 'chr1:1000-2000'
                    """
                )
                # Only variant 1 intersects both feature and the specified range
                assert len(result) == 1
                assert result.iloc[0]["name"] == "var1"
