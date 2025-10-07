import tempfile

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
