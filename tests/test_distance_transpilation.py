"""Transpilation tests for DISTANCE operator SQL generation.

Tests verify that DISTANCE() is correctly transpiled to SQL CASE expressions
across different SQL dialects (DuckDB, SQLite, PostgreSQL).
"""

from sqlglot import parse_one

from giql.dialect import GIQLDialect
from giql.generators import BaseGIQLGenerator
from giql.generators import GIQLDuckDBGenerator


class TestDistanceTranspilation:
    """Tests for DISTANCE SQL generation across dialects."""

    def test_distance_transpilation_duckdb(self):
        """
        GIVEN a GIQL query with DISTANCE()
        WHEN transpiling to DuckDB SQL
        THEN should generate complete CASE expression with distance logic
        """
        sql = """
        SELECT DISTANCE(a.interval, b.interval) as dist
        FROM features_a a CROSS JOIN features_b b
        """

        ast = parse_one(sql, dialect=GIQLDialect)
        generator = GIQLDuckDBGenerator()
        output = generator.generate(ast)

        expected = """SELECT CASE WHEN a."chromosome" != b."chromosome" THEN NULL WHEN a."start_pos" < b."end_pos" AND a."end_pos" > b."start_pos" THEN 0 WHEN a."end_pos" <= b."start_pos" THEN (b."start_pos" - a."end_pos") ELSE (a."start_pos" - b."end_pos") END AS dist FROM features_a AS a CROSS JOIN features_b AS b"""

        assert output == expected, f"Expected:\n{expected}\n\nGot:\n{output}"

    def test_distance_transpilation_sqlite(self):
        """
        GIVEN a GIQL query with DISTANCE()
        WHEN transpiling to SQLite SQL
        THEN should generate complete compatible SQL CASE expression
        """
        sql = """
        SELECT DISTANCE(a.interval, b.interval) as dist
        FROM features_a a, features_b b
        """

        ast = parse_one(sql, dialect=GIQLDialect)
        generator = BaseGIQLGenerator()
        output = generator.generate(ast)

        expected = """SELECT CASE WHEN a."chromosome" != b."chromosome" THEN NULL WHEN a."start_pos" < b."end_pos" AND a."end_pos" > b."start_pos" THEN 0 WHEN a."end_pos" <= b."start_pos" THEN (b."start_pos" - a."end_pos") ELSE (a."start_pos" - b."end_pos") END AS dist FROM features_a AS a, features_b AS b"""

        assert output == expected, f"Expected:\n{expected}\n\nGot:\n{output}"

    def test_distance_transpilation_postgres(self):
        """
        GIVEN a GIQL query with DISTANCE()
        WHEN transpiling to PostgreSQL SQL
        THEN should generate complete compatible SQL CASE expression
        """
        sql = """
        SELECT DISTANCE(a.interval, b.interval) as dist
        FROM features_a a CROSS JOIN features_b b
        """

        ast = parse_one(sql, dialect=GIQLDialect)
        generator = BaseGIQLGenerator()
        output = generator.generate(ast)

        expected = """SELECT CASE WHEN a."chromosome" != b."chromosome" THEN NULL WHEN a."start_pos" < b."end_pos" AND a."end_pos" > b."start_pos" THEN 0 WHEN a."end_pos" <= b."start_pos" THEN (b."start_pos" - a."end_pos") ELSE (a."start_pos" - b."end_pos") END AS dist FROM features_a AS a CROSS JOIN features_b AS b"""

        assert output == expected, f"Expected:\n{expected}\n\nGot:\n{output}"

    def test_distance_transpilation_signed_duckdb(self):
        """
        GIVEN a GIQL query with DISTANCE(..., signed=true)
        WHEN transpiling to DuckDB SQL
        THEN should generate CASE expression with signed distances
            (negative for upstream, positive for downstream)
        """
        sql = """
        SELECT DISTANCE(a.interval, b.interval, signed=true) as dist
        FROM features_a a CROSS JOIN features_b b
        """

        ast = parse_one(sql, dialect=GIQLDialect)
        generator = GIQLDuckDBGenerator()
        output = generator.generate(ast)

        # Signed distance: upstream (B before A) returns negative,
        # downstream (B after A) returns positive
        expected = (
            'SELECT CASE WHEN a."chromosome" != b."chromosome" THEN NULL '
            'WHEN a."start_pos" < b."end_pos" AND a."end_pos" > b."start_pos" THEN 0 '
            'WHEN a."end_pos" <= b."start_pos" THEN (b."start_pos" - a."end_pos") '
            'ELSE -(a."start_pos" - b."end_pos") END AS dist '
            "FROM features_a AS a CROSS JOIN features_b AS b"
        )

        assert output == expected, f"Expected:\n{expected}\n\nGot:\n{output}"
