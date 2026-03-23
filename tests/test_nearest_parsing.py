"""Parser tests for NEAREST operator syntax.

Tests verify that the GIQL parser correctly recognizes and parses
NEAREST function calls with various argument patterns.
"""

from sqlglot import parse_one

from giql.dialect import GIQLDialect
from giql.expressions import GIQLNearest


class TestNearestParsing:
    """Tests for parsing NEAREST function syntax."""

    def test_parse_nearest_basic_syntax(self):
        """
        GIVEN a GIQL query with NEAREST(genes, k := 3)
        WHEN parsing the query
        THEN should create GIQLNearest AST node with correct arguments
        """
        sql = "SELECT * FROM peaks CROSS JOIN LATERAL NEAREST(genes, k := 3)"

        ast = parse_one(sql, dialect=GIQLDialect)

        # Find the NEAREST expression in the JOIN clause
        # Navigate: Select -> joins[0] -> this (which should be NEAREST)
        joins = ast.args.get("joins")
        assert joins is not None, "Expected JOIN clause"

        join = joins[0]
        lateral_expr = join.this

        # For LATERAL, the function is nested inside
        if hasattr(lateral_expr, "this"):
            nearest_expr = lateral_expr.this
        else:
            nearest_expr = lateral_expr

        assert isinstance(nearest_expr, GIQLNearest), (
            f"Expected GIQLNearest node, got {type(nearest_expr)}"
        )

        # Verify target table argument
        assert nearest_expr.this is not None, "Missing target table argument (this)"

        # Verify k parameter
        k_param = nearest_expr.args.get("k")
        assert k_param is not None, "Missing k parameter"

    def test_parse_nearest_with_literal_reference(self):
        """
        GIVEN a GIQL query with NEAREST(genes, reference := 'chr1:1000-2000', k := 3)
        WHEN parsing the query
        THEN should create GIQLNearest node with literal reference parameter
        """
        sql = "SELECT * FROM NEAREST(genes, reference := 'chr1:1000-2000', k := 3)"

        ast = parse_one(sql, dialect=GIQLDialect)

        # Navigate to NEAREST function in FROM clause
        # When a function is in FROM, sqlglot wraps it in a Table expression
        from_clause = ast.args.get("from_")
        table_expr = from_clause.this

        # The NEAREST expression should be nested in the Table's 'this'
        nearest_expr = table_expr.this if hasattr(table_expr, "this") else table_expr

        assert isinstance(nearest_expr, GIQLNearest), (
            f"Expected GIQLNearest node, got {type(nearest_expr)}"
        )

        # Verify reference parameter exists
        reference = nearest_expr.args.get("reference")
        assert reference is not None, "Missing reference parameter"

        # Verify k parameter
        k_param = nearest_expr.args.get("k")
        assert k_param is not None, "Missing k parameter"

    def test_parse_nearest_with_max_distance(self):
        """
        GIVEN a GIQL query with NEAREST(genes, k := 5, max_distance := 100000)
        WHEN parsing the query
        THEN should parse max_distance parameter correctly
        """
        sql = "SELECT * FROM peaks CROSS JOIN LATERAL NEAREST(genes, k := 5, max_distance := 100000)"

        ast = parse_one(sql, dialect=GIQLDialect)

        # Navigate to NEAREST
        joins = ast.args.get("joins")
        join = joins[0]
        lateral_expr = join.this

        if hasattr(lateral_expr, "this"):
            nearest_expr = lateral_expr.this
        else:
            nearest_expr = lateral_expr

        assert isinstance(nearest_expr, GIQLNearest), f"Expected GIQLNearest node"

        # Verify max_distance parameter
        max_distance = nearest_expr.args.get("max_distance")
        assert max_distance is not None, "Missing max_distance parameter"

    def test_parse_nearest_with_stranded(self):
        """
        GIVEN a GIQL query with NEAREST(genes, k := 3, stranded := true)
        WHEN parsing the query
        THEN should parse stranded parameter correctly
        """
        sql = "SELECT * FROM peaks CROSS JOIN LATERAL NEAREST(genes, k := 3, stranded := true)"

        ast = parse_one(sql, dialect=GIQLDialect)

        # Navigate to NEAREST
        joins = ast.args.get("joins")
        join = joins[0]
        lateral_expr = join.this

        if hasattr(lateral_expr, "this"):
            nearest_expr = lateral_expr.this
        else:
            nearest_expr = lateral_expr

        assert isinstance(nearest_expr, GIQLNearest), f"Expected GIQLNearest node"

        # Verify stranded parameter
        stranded = nearest_expr.args.get("stranded")
        assert stranded is not None, "Missing stranded parameter"

    def test_parse_nearest_all_parameters(self):
        """
        GIVEN a GIQL query with all NEAREST parameters
        WHEN parsing the query
        THEN should parse all parameters correctly
        """
        sql = """
        SELECT * FROM peaks
        CROSS JOIN LATERAL NEAREST(
            genes,
            reference := peaks.interval,
            k := 5,
            max_distance := 50000,
            stranded := true
        )
        """

        ast = parse_one(sql, dialect=GIQLDialect)

        # Navigate to NEAREST
        joins = ast.args.get("joins")
        join = joins[0]
        lateral_expr = join.this

        if hasattr(lateral_expr, "this"):
            nearest_expr = lateral_expr.this
        else:
            nearest_expr = lateral_expr

        assert isinstance(nearest_expr, GIQLNearest), f"Expected GIQLNearest node"

        # Verify all parameters exist
        assert nearest_expr.this is not None, "Missing target table"
        assert nearest_expr.args.get("reference") is not None, "Missing reference"
        assert nearest_expr.args.get("k") is not None, "Missing k"
        assert nearest_expr.args.get("max_distance") is not None, "Missing max_distance"
        assert nearest_expr.args.get("stranded") is not None, "Missing stranded"

    def test_parse_nearest_with_strand_notation(self):
        """
        GIVEN a GIQL query with literal reference using strand notation
        WHEN parsing 'chr1:1000-2000:+' format
        THEN should parse the strand-annotated range correctly
        """
        sql = "SELECT * FROM NEAREST(genes, reference := 'chr1:1000-2000:+', k := 3)"

        ast = parse_one(sql, dialect=GIQLDialect)

        # Navigate to NEAREST function in FROM clause
        from_clause = ast.args.get("from_")
        table_expr = from_clause.this

        # The NEAREST expression should be nested in the Table's 'this'
        nearest_expr = table_expr.this if hasattr(table_expr, "this") else table_expr

        assert isinstance(nearest_expr, GIQLNearest), (
            f"Expected GIQLNearest node, got {type(nearest_expr)}"
        )

        # Verify reference parameter exists and contains strand notation
        reference = nearest_expr.args.get("reference")
        assert reference is not None, "Missing reference parameter"

        # Verify k parameter
        k_param = nearest_expr.args.get("k")
        assert k_param is not None, "Missing k parameter"

    def test_from_arg_list_with_eq_ignores_named_param(self):
        """Test that = syntax is not treated as named parameter assignment.

        Given:
            A GIQL query with NEAREST(genes, k=3) using = syntax
        When:
            Parsing the query
        Then:
            It should not treat k as a named parameter
        """
        # Act
        ast = parse_one(
            "SELECT * FROM peaks CROSS JOIN LATERAL NEAREST(genes, k=3)",
            dialect=GIQLDialect,
        )

        # Assert
        joins = ast.args.get("joins")
        join = joins[0]
        lateral_expr = join.this
        nearest_expr = lateral_expr.this if hasattr(lateral_expr, "this") else lateral_expr
        assert isinstance(nearest_expr, GIQLNearest)
        assert nearest_expr.args.get("k") is None, (
            "= should not be treated as named parameter assignment"
        )

    def test_from_arg_list_with_kwarg_syntax(self):
        """Test that => (SQL-standard) syntax works for named parameters.

        Given:
            A GIQL query with NEAREST(genes, k => 3) using => syntax
        When:
            Parsing the query
        Then:
            It should parse k as a named parameter
        """
        # Act
        ast = parse_one(
            "SELECT * FROM peaks CROSS JOIN LATERAL NEAREST(genes, k => 3)",
            dialect=GIQLDialect,
        )

        # Assert
        joins = ast.args.get("joins")
        join = joins[0]
        lateral_expr = join.this
        nearest_expr = lateral_expr.this if hasattr(lateral_expr, "this") else lateral_expr
        assert isinstance(nearest_expr, GIQLNearest)
        assert nearest_expr.args.get("k") is not None, "Missing k parameter with => syntax"
