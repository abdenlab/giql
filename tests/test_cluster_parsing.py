"""Parser tests for CLUSTER operator syntax.

Tests verify that the GIQL parser correctly recognizes and parses
CLUSTER function calls with various argument patterns, and that
only := and => are accepted for named parameter binding.
"""

import pytest
from sqlglot import exp
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

    def test_from_arg_list_with_predicate(self):
        """Test that a predicate named argument is captured on the node.

        Given:
            A GIQL query with CLUSTER(interval, predicate := depth = prev.depth).
        When:
            Parsing the query.
        Then:
            It should attach the predicate as an equality expression in args.
        """
        # Act
        ast = parse_one(
            "SELECT *, CLUSTER(interval, predicate := depth = prev.depth) AS cid "
            "FROM peaks",
            dialect=GIQLDialect,
        )

        # Assert
        cluster_expr = ast.expressions[1].this
        assert isinstance(cluster_expr, GIQLCluster)
        assert isinstance(cluster_expr.args.get("predicate"), exp.EQ)

    def test_from_arg_list_with_predicate_prev_qualifier(self):
        """Test that a prev. qualifier parses as a column on the predecessor.

        Given:
            A predicate CLUSTER(interval, predicate := depth = prev.depth)
            referencing the predecessor row with the prev. qualifier.
        When:
            Parsing the query.
        Then:
            It should parse prev.depth as a column whose table is prev.
        """
        # Act
        ast = parse_one(
            "SELECT *, CLUSTER(interval, predicate := depth = prev.depth) AS cid "
            "FROM peaks",
            dialect=GIQLDialect,
        )

        # Assert
        predicate = ast.expressions[1].this.args["predicate"]
        prev_ref = predicate.expression
        assert isinstance(prev_ref, exp.Column)
        assert prev_ref.table == "prev"
        assert prev_ref.name == "depth"

    def test_from_arg_list_with_predicate_kwarg_syntax(self):
        """Test that the => kwarg form also binds the predicate argument.

        Given:
            A CLUSTER call using CLUSTER(interval, predicate => score = prev.score)
            with the => kwarg form rather than :=.
        When:
            Parsing the query.
        Then:
            It should attach the predicate as an equality expression in args.
        """
        # Act
        ast = parse_one(
            "SELECT *, CLUSTER(interval, predicate => score = prev.score) AS cid "
            "FROM peaks",
            dialect=GIQLDialect,
        )

        # Assert
        cluster_expr = ast.expressions[1].this
        assert isinstance(cluster_expr, GIQLCluster)
        assert isinstance(cluster_expr.args.get("predicate"), exp.EQ)

    def test_from_arg_list_with_predicate_and_positional_distance(self):
        """Test that a predicate composes with positional distance and stranded.

        Given:
            A query mixing a positional distance, stranded :=, and predicate :=
            on a single CLUSTER call.
        When:
            Parsing the query.
        Then:
            It should retain distance, stranded, and predicate together in args.
        """
        # Act
        ast = parse_one(
            "SELECT *, CLUSTER(interval, 1000, stranded := true, "
            "predicate := name = prev.name) AS cid FROM peaks",
            dialect=GIQLDialect,
        )

        # Assert
        cluster_expr = ast.expressions[1].this
        assert isinstance(cluster_expr, GIQLCluster)
        assert cluster_expr.args.get("distance") is not None
        assert cluster_expr.args.get("stranded") is not None
        assert cluster_expr.args.get("predicate") is not None
