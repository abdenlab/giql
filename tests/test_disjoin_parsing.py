"""Parser tests for the DISJOIN operator syntax.

Tests verify that the GIQL parser recognizes DISJOIN table-function calls
and maps positional and named arguments onto the GIQLDisjoin AST node.
"""

import pytest
from sqlglot import exp
from sqlglot import parse_one
from sqlglot.errors import ParseError

from giql.dialect import GIQLDialect
from giql.expressions import GIQLDisjoin


def _disjoin_node(sql: str) -> GIQLDisjoin:
    """Parse a GIQL query and return its first GIQLDisjoin node."""
    ast = parse_one(sql, dialect=GIQLDialect)
    node = ast.find(GIQLDisjoin)
    assert node is not None, f"Expected a GIQLDisjoin node in: {sql}"
    return node


class TestDisjoinParsing:
    """Tests for parsing DISJOIN function syntax."""

    def test_from_arg_list_should_expose_documented_arg_types(self):
        """Test that GIQLDisjoin declares the documented argument schema.

        Given:
            The GIQLDisjoin expression class.
        When:
            Inspecting its arg_types mapping.
        Then:
            It should require `this` and make `reference` optional.
        """
        # Arrange, act, & assert
        assert GIQLDisjoin.arg_types == {"this": True, "reference": False}

    def test_from_arg_list_should_map_positional_arg_to_target(self):
        """Test that a positional argument binds to the target set.

        Given:
            A GIQL query with DISJOIN(features).
        When:
            Parsing the query.
        Then:
            It should create a GIQLDisjoin node with the target set as `this`.
        """
        # Arrange & act
        node = _disjoin_node("SELECT * FROM DISJOIN(features)")

        # Assert
        assert isinstance(node, GIQLDisjoin)
        assert node.this is not None
        assert node.this.name == "features"

    def test_from_arg_list_should_leave_reference_unset_when_omitted(self):
        """Test that an omitted reference argument stays unset.

        Given:
            A GIQL query with DISJOIN(features) and no reference argument.
        When:
            Parsing the query.
        Then:
            It should leave the GIQLDisjoin node's reference argument unset.
        """
        # Arrange & act
        node = _disjoin_node("SELECT * FROM DISJOIN(features)")

        # Assert
        assert node.args.get("reference") is None

    def test_from_arg_list_should_reject_extra_positional_argument(self):
        """Test that a bare second positional argument is rejected.

        Given:
            A GIQL query with DISJOIN(features, refs) passing the reference
            set as a bare second positional argument.
        When:
            Parsing the query.
        Then:
            It should raise a ParseError naming the positional-argument limit.
        """
        # Arrange, act, & assert
        with pytest.raises(ParseError, match="at most one positional"):
            parse_one("SELECT * FROM DISJOIN(features, refs)", dialect=GIQLDialect)

    def test_from_arg_list_should_reject_unknown_named_argument(self):
        """Test that an unrecognized named argument is rejected.

        Given:
            A GIQL query with DISJOIN(features, mask := refs) naming an
            argument that is not part of the DISJOIN schema.
        When:
            Parsing the query.
        Then:
            It should raise a ParseError naming the unexpected argument.
        """
        # Arrange, act, & assert
        with pytest.raises(ParseError, match="unexpected named argument"):
            parse_one(
                "SELECT * FROM DISJOIN(features, mask := refs)", dialect=GIQLDialect
            )

    def test_from_arg_list_should_reject_missing_target(self):
        """Test that a DISJOIN call with no target argument is rejected.

        Given:
            A GIQL query with DISJOIN(reference := refs) supplying only a
            named reference and no positional target.
        When:
            Parsing the query.
        Then:
            It should raise a ParseError naming the required target argument.
        """
        # Arrange, act, & assert
        with pytest.raises(ParseError, match="requires a target table"):
            parse_one(
                "SELECT * FROM DISJOIN(reference := refs)", dialect=GIQLDialect
            )

    def test_from_arg_list_should_map_reference_when_named_with_walrus(self):
        """Test that a walrus-named reference argument is carried.

        Given:
            A GIQL query with DISJOIN(features, reference := refs).
        When:
            Parsing the query.
        Then:
            It should carry the reference argument on the GIQLDisjoin node.
        """
        # Arrange & act
        node = _disjoin_node("SELECT * FROM DISJOIN(features, reference := refs)")

        # Assert
        reference = node.args.get("reference")
        assert reference is not None
        assert reference.name == "refs"

    def test_from_arg_list_should_map_reference_when_named_with_arrow(self):
        """Test that an arrow-named reference argument is carried.

        Given:
            A GIQL query with DISJOIN(features, reference => refs).
        When:
            Parsing the query.
        Then:
            It should carry the reference argument on the GIQLDisjoin node.
        """
        # Arrange & act
        node = _disjoin_node("SELECT * FROM DISJOIN(features, reference => refs)")

        # Assert
        assert node.args.get("reference") is not None

    def test_from_arg_list_should_carry_subquery_reference(self):
        """Test that a subquery reference is carried as a nested SELECT.

        Given:
            A GIQL query whose DISJOIN reference is a subquery.
        When:
            Parsing the query.
        Then:
            It should carry a reference argument containing a nested SELECT.
        """
        # Arrange & act
        node = _disjoin_node(
            "SELECT * FROM DISJOIN(features, reference := (SELECT * FROM refs))"
        )

        # Assert
        reference = node.args.get("reference")
        assert reference is not None
        assert reference.find(exp.Select) is not None
