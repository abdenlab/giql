"""Transpilation tests for DISTANCE operator SQL generation.

Tests verify that DISTANCE() is correctly transpiled to SQL CASE expressions.
"""

from sqlglot import parse_one

import giql  # noqa: F401  (ensures the built-in expanders are registered)
from giql import transpile
from giql.canonicalizer import canonicalize_coordinates
from giql.dialect import GIQLDialect
from giql.expander import ExpandOperators
from giql.resolver import resolve_operator_refs
from giql.table import Tables
from giql.targets import GenericTarget


def _generate(sql: str, tables: Tables | None = None) -> str:
    """Parse, run normalization passes 1-3, then generate SQL.

    DISTANCE operand resolution and coordinate canonicalization moved out of the
    emitter and into the ResolveOperatorRefs / CanonicalizeCoordinates passes
    (epic #114, issues #119 / #123). DISTANCE generation itself then moved onto
    the registry's AST-expansion pass (epic #137, issue #140), so the operator
    node must be expanded before generation too. Emitter-level tests run all
    three passes, exactly as :func:`giql.transpile.transpile` does, rather than
    calling ``generate`` on a bare parsed AST.
    """
    tables = tables or Tables()
    ast = parse_one(sql, dialect=GIQLDialect)
    ast = resolve_operator_refs(ast, tables)
    ast = canonicalize_coordinates(ast)
    ast = ExpandOperators(GenericTarget(), tables).transform(ast)
    return ast.sql()


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

        output = _generate(sql)

        expected = """SELECT CASE WHEN a."chrom" <> b."chrom" THEN NULL WHEN a."start" < b."end" AND a."end" > b."start" THEN 0 WHEN a."end" <= b."start" THEN (b."start" - a."end" + 1) ELSE (a."start" - b."end" + 1) END AS dist FROM features_a AS a CROSS JOIN features_b AS b"""

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

        output = _generate(sql)

        expected = """SELECT CASE WHEN a."chrom" <> b."chrom" THEN NULL WHEN a."start" < b."end" AND a."end" > b."start" THEN 0 WHEN a."end" <= b."start" THEN (b."start" - a."end" + 1) ELSE (a."start" - b."end" + 1) END AS dist FROM features_a AS a, features_b AS b"""

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

        output = _generate(sql)

        expected = """SELECT CASE WHEN a."chrom" <> b."chrom" THEN NULL WHEN a."start" < b."end" AND a."end" > b."start" THEN 0 WHEN a."end" <= b."start" THEN (b."start" - a."end" + 1) ELSE (a."start" - b."end" + 1) END AS dist FROM features_a AS a CROSS JOIN features_b AS b"""

        assert output == expected, f"Expected:\n{expected}\n\nGot:\n{output}"

    def test_distance_resolver_path_matches_direct_generation(self):
        """
        GIVEN a DISTANCE query over default-convention tables
        WHEN transpiling through the full pipeline (with the tables registered)
            versus running the two normalization passes with no tables registered
        THEN both paths should emit byte-identical SQL, proving the qualified
            ResolvedColumn metadata resolves identically whether or not the
            operand's relation is registered
        """
        query = (
            "SELECT DISTANCE(a.interval, b.interval) AS dist "
            "FROM features_a a, features_b b"
        )

        via_transpile = transpile(query, tables=["features_a", "features_b"])
        via_generate = _generate(query)

        assert via_transpile == via_generate

    def test_distance_transpilation_signed_duckdb(self):
        """
        GIVEN a GIQL query with DISTANCE(..., signed := true)
        WHEN transpiling to DuckDB SQL
        THEN should generate CASE expression with signed distances
            (negative for upstream, positive for downstream)
        """
        sql = """
        SELECT DISTANCE(a.interval, b.interval, signed := true) as dist
        FROM features_a a CROSS JOIN features_b b
        """

        output = _generate(sql)

        # Signed distance: upstream (B before A) returns negative,
        # downstream (B after A) returns positive
        expected = (
            'SELECT CASE WHEN a."chrom" <> b."chrom" THEN NULL '
            'WHEN a."start" < b."end" AND a."end" > b."start" THEN 0 '
            'WHEN a."end" <= b."start" THEN (b."start" - a."end" + 1) '
            'ELSE -(a."start" - b."end" + 1) END AS dist '
            "FROM features_a AS a CROSS JOIN features_b AS b"
        )

        assert output == expected, f"Expected:\n{expected}\n\nGot:\n{output}"
