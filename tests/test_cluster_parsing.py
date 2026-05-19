"""Parser tests for CLUSTER operator syntax.

Tests verify that the GIQL parser correctly recognizes and parses
CLUSTER function calls with various argument patterns, and that
only := and => are accepted for named parameter binding.
"""

import pytest
from sqlglot import parse_one
from sqlglot.errors import ParseError

from giql.dialect import GIQLDialect
from giql.expressions import GIQLCluster


class TestClusterParsing:
    """Tests for parsing CLUSTER function syntax."""

    def test_from_arg_list_with_property_eq_syntax(self):
        """
        GIVEN a GIQL query with CLUSTER(interval, stranded := true)
        WHEN parsing the query
        THEN should create GIQLCluster node with stranded in args
        """
        # Act
        ast = parse_one(
            "SELECT *, CLUSTER(interval, stranded := true) AS cluster_id FROM peaks",
            dialect=GIQLDialect,
        )

        # Assert
        alias_expr = ast.expressions[1]
        cluster_expr = alias_expr.this
        assert isinstance(cluster_expr, GIQLCluster), (
            f"Expected GIQLCluster node, got {type(cluster_expr)}"
        )
        assert cluster_expr.args.get("stranded") is not None, (
            "Missing stranded parameter"
        )

    def test_from_arg_list_with_eq_as_positional(self):
        """
        GIVEN a GIQL query with CLUSTER(interval, stranded=true) using = syntax
        WHEN parsing the query
        THEN should not treat stranded as a named parameter
        """
        # Act
        ast = parse_one(
            "SELECT *, CLUSTER(interval, stranded=true) AS cluster_id FROM peaks",
            dialect=GIQLDialect,
        )

        # Assert
        alias_expr = ast.expressions[1]
        cluster_expr = alias_expr.this
        assert isinstance(cluster_expr, GIQLCluster)
        assert cluster_expr.args.get("stranded") is None, (
            "= should not be treated as named parameter assignment"
        )

    def test_from_arg_list_with_kwarg_syntax(self):
        """
        GIVEN a GIQL query with CLUSTER(interval, stranded => true) using => syntax
        WHEN parsing the query
        THEN should parse stranded as a named parameter
        """
        # Act
        ast = parse_one(
            "SELECT *, CLUSTER(interval, stranded => true) AS cluster_id FROM peaks",
            dialect=GIQLDialect,
        )

        # Assert
        alias_expr = ast.expressions[1]
        cluster_expr = alias_expr.this
        assert isinstance(cluster_expr, GIQLCluster)
        assert cluster_expr.args.get("stranded") is not None, (
            "Missing stranded parameter with => syntax"
        )

    def test_from_arg_list_should_reject_missing_target(self):
        """Test that a CLUSTER call with no interval argument is rejected.

        Given:
            A GIQL query with CLUSTER(stranded := true) supplying only a
            named argument and no positional genomic interval column.
        When:
            Parsing the query.
        Then:
            It should raise a ParseError naming the required column argument.
        """
        # Arrange, act, & assert
        with pytest.raises(ParseError, match="requires a genomic interval"):
            parse_one(
                "SELECT CLUSTER(stranded := true) AS cluster_id FROM peaks",
                dialect=GIQLDialect,
            )
