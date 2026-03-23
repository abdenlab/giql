"""Parser tests for MERGE operator syntax.

Tests verify that the GIQL parser correctly recognizes and parses
MERGE function calls, and that only := and => are accepted for
named parameter binding.
"""

from sqlglot import parse_one

from giql.dialect import GIQLDialect
from giql.expressions import GIQLMerge


class TestMergeParsing:
    """Tests for parsing MERGE function syntax."""

    def test_from_arg_list_with_eq_ignores_named_param(self):
        """Test that = syntax is not treated as named parameter assignment.

        Given:
            A GIQL query with MERGE(interval, stranded=true) using = syntax
        When:
            Parsing the query
        Then:
            It should not treat stranded as a named parameter
        """
        # Act
        ast = parse_one(
            "SELECT MERGE(interval, stranded=true) FROM peaks",
            dialect=GIQLDialect,
        )

        # Assert
        merge_expr = ast.expressions[0]
        assert isinstance(merge_expr, GIQLMerge)
        assert merge_expr.args.get("stranded") is None, (
            "= should not be treated as named parameter assignment"
        )

    def test_from_arg_list_with_kwarg_syntax(self):
        """Test that => (SQL-standard) syntax works for named parameters.

        Given:
            A GIQL query with MERGE(interval, stranded => true) using => syntax
        When:
            Parsing the query
        Then:
            It should parse stranded as a named parameter
        """
        # Act
        ast = parse_one(
            "SELECT MERGE(interval, stranded => true) FROM peaks",
            dialect=GIQLDialect,
        )

        # Assert
        merge_expr = ast.expressions[0]
        assert isinstance(merge_expr, GIQLMerge)
        assert merge_expr.args.get("stranded") is not None, (
            "Missing stranded parameter with => syntax"
        )
