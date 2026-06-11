"""Parser tests for MERGE operator syntax.

Tests verify that the GIQL parser correctly recognizes and parses
MERGE function calls, and that only := and => are accepted for
named parameter binding.
"""

import pytest
from sqlglot import exp
from sqlglot import parse_one
from sqlglot.errors import ParseError

from giql.dialect import GIQLDialect
from giql.expressions import GIQLMerge


class TestMergeParsing:
    """Tests for parsing MERGE function syntax."""

    def test_from_arg_list_with_property_eq_syntax(self):
        """
        GIVEN a GIQL query with MERGE(interval, stranded := true)
        WHEN parsing the query
        THEN should parse stranded as a named parameter
        """
        # Act
        ast = parse_one(
            "SELECT MERGE(interval, stranded := true) FROM peaks",
            dialect=GIQLDialect,
        )

        # Assert
        merge_expr = ast.expressions[0]
        assert isinstance(merge_expr, GIQLMerge)
        assert merge_expr.args.get("stranded") is not None, "Missing stranded parameter"

    def test_from_arg_list_with_eq_as_positional(self):
        """
        GIVEN a GIQL query with MERGE(interval, stranded=true) using = syntax
        WHEN parsing the query
        THEN should not treat stranded as a named parameter
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
        """
        GIVEN a GIQL query with MERGE(interval, stranded => true) using => syntax
        WHEN parsing the query
        THEN should parse stranded as a named parameter
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

    def test_from_arg_list_should_reject_missing_target(self):
        """Test that a MERGE call with no interval argument is rejected.

        Given:
            A GIQL query with MERGE(stranded := true) supplying only a named
            argument and no positional genomic interval column.
        When:
            Parsing the query.
        Then:
            It should raise a ParseError naming the required column argument.
        """
        # Arrange, act, & assert
        with pytest.raises(ParseError, match="requires a genomic interval"):
            parse_one("SELECT MERGE(stranded := true) FROM peaks", dialect=GIQLDialect)

    def test_from_arg_list_with_predicate(self):
        """Test that a predicate named argument is captured on the node.

        Given:
            A GIQL query with MERGE(interval, predicate := depth = prev.depth).
        When:
            Parsing the query.
        Then:
            It should attach the predicate as an equality expression in args.
        """
        # Act
        ast = parse_one(
            "SELECT MERGE(interval, predicate := depth = prev.depth) FROM peaks",
            dialect=GIQLDialect,
        )

        # Assert
        merge_expr = ast.expressions[0]
        assert isinstance(merge_expr, GIQLMerge)
        assert isinstance(merge_expr.args.get("predicate"), exp.EQ)

    def test_from_arg_list_with_predicate_kwarg_syntax(self):
        """Test that the => kwarg form also binds the predicate argument.

        Given:
            A MERGE call using MERGE(interval, predicate => score = prev.score)
            with the => kwarg form rather than :=.
        When:
            Parsing the query.
        Then:
            It should attach the predicate as an equality expression in args.
        """
        # Act
        ast = parse_one(
            "SELECT MERGE(interval, predicate => score = prev.score) FROM peaks",
            dialect=GIQLDialect,
        )

        # Assert
        merge_expr = ast.expressions[0]
        assert isinstance(merge_expr, GIQLMerge)
        assert isinstance(merge_expr.args.get("predicate"), exp.EQ)

    def test_from_arg_list_with_compound_predicate(self):
        """Test that a multi-term predicate parses as a conjunction.

        Given:
            A query MERGE(interval, predicate := strand = prev.strand AND
            name = prev.name) joining two pairwise comparisons with AND.
        When:
            Parsing the query.
        Then:
            It should attach the predicate as an AND of two prev. comparisons.
        """
        # Act
        ast = parse_one(
            "SELECT MERGE(interval, predicate := strand = prev.strand "
            "AND name = prev.name) FROM peaks",
            dialect=GIQLDialect,
        )

        # Assert
        predicate = ast.expressions[0].args["predicate"]
        assert isinstance(predicate, exp.And)
        prev_tables = {col.table for col in predicate.find_all(exp.Column)}
        assert "prev" in prev_tables
