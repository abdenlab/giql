"""Parser tests for DISTANCE operator syntax.

Tests verify that the GIQL parser correctly recognizes and parses
DISTANCE function calls with various argument patterns.
"""

from sqlglot import parse_one

from giql.dialect import GIQLDialect
from giql.expressions import GIQLDistance


class TestDistanceParsing:
    """Tests for parsing DISTANCE function syntax."""

    def test_parse_distance_with_column_to_column(self):
        """
        GIVEN a GIQL query with DISTANCE(a.interval, b.interval)
        WHEN parsing the query
        THEN should create GIQLDistance AST node with correct arguments
        """
        sql = "SELECT DISTANCE(a.interval, b.interval) FROM features_a a, features_b b"

        ast = parse_one(sql, dialect=GIQLDialect)

        # Find the DISTANCE expression in the SELECT clause
        select_expr = ast.expressions[0]
        assert isinstance(select_expr, GIQLDistance), (
            f"Expected GIQLDistance node, got {type(select_expr)}"
        )

        # Verify arguments are present
        assert select_expr.this is not None, "Missing first argument (this)"
        assert select_expr.args.get("expression") is not None, (
            "Missing second argument (expression)"
        )

    def test_parse_distance_with_literal_range(self):
        """
        GIVEN a GIQL query with DISTANCE(a.interval, 'chr1:100-200')
        WHEN parsing the query
        THEN should create GIQLDistance node with column and literal range
        """
        sql = "SELECT DISTANCE(a.interval, 'chr1:100-200') FROM features a"

        ast = parse_one(sql, dialect=GIQLDialect)

        # Find the DISTANCE expression
        select_expr = ast.expressions[0]
        assert isinstance(select_expr, GIQLDistance), (
            f"Expected GIQLDistance node, got {type(select_expr)}"
        )

        # Verify both arguments present
        assert select_expr.this is not None, "Missing first argument"
        assert select_expr.args.get("expression") is not None, "Missing second argument"

        # Second argument should be a literal string
        second_arg = select_expr.args["expression"]
        assert "chr1" in str(second_arg).lower(), "Expected chromosome in literal range"

    def test_from_arg_list_with_eq_ignores_named_param(self):
        """Test that = syntax is not treated as named parameter assignment.

        Given:
            A GIQL query with DISTANCE(a.interval, b.interval, stranded=true) using = syntax
        When:
            Parsing the query
        Then:
            It should not treat stranded as a named parameter
        """
        # Act
        ast = parse_one(
            "SELECT DISTANCE(a.interval, b.interval, stranded=true) FROM features_a a, features_b b",
            dialect=GIQLDialect,
        )

        # Assert
        select_expr = ast.expressions[0]
        assert isinstance(select_expr, GIQLDistance)
        assert select_expr.args.get("stranded") is None, (
            "= should not be treated as named parameter assignment"
        )

    def test_from_arg_list_with_kwarg_syntax(self):
        """Test that => (SQL-standard) syntax works for named parameters.

        Given:
            A GIQL query with DISTANCE(a.interval, b.interval, stranded => true) using => syntax
        When:
            Parsing the query
        Then:
            It should parse stranded as a named parameter
        """
        # Act
        ast = parse_one(
            "SELECT DISTANCE(a.interval, b.interval, stranded => true) FROM features_a a, features_b b",
            dialect=GIQLDialect,
        )

        # Assert
        select_expr = ast.expressions[0]
        assert isinstance(select_expr, GIQLDistance)
        assert select_expr.args.get("stranded") is not None, (
            "Missing stranded parameter with => syntax"
        )
