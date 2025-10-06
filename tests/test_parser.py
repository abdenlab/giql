import pytest
from sqlglot import parse_one

from giql.dialect import GIQLDialect
from giql.expressions import Contains
from giql.expressions import Intersects
from giql.expressions import SpatialSetPredicate
from giql.expressions import Within


class TestParser:
    def test_parse_simple_intersects(self):
        """
        GIVEN a SQL query with INTERSECTS operator
        WHEN parsing the query
        THEN should create an Intersects AST node
        """
        sql = "SELECT * FROM variants WHERE position INTERSECTS 'chr1:1000-2000'"
        ast = parse_one(sql, dialect=GIQLDialect)

        # Find the INTERSECTS node
        intersects_node = None
        for node in ast.walk():
            if isinstance(node, Intersects):
                intersects_node = node
                break

        assert intersects_node is not None

    def test_parse_contains(self):
        """
        GIVEN a SQL query with CONTAINS operator
        WHEN parsing the query
        THEN should create a Contains AST node
        """
        sql = "SELECT * FROM variants WHERE position CONTAINS 'chr1:1500'"
        ast = parse_one(sql, dialect=GIQLDialect)

        contains_node = None
        for node in ast.walk():
            if isinstance(node, Contains):
                contains_node = node
                break

        assert contains_node is not None

    def test_parse_within(self):
        """
        GIVEN a SQL query with WITHIN operator
        WHEN parsing the query
        THEN should create a Within AST node
        """
        sql = "SELECT * FROM variants WHERE position WITHIN 'chr1:1000-5000'"
        ast = parse_one(sql, dialect=GIQLDialect)

        within_node = None
        for node in ast.walk():
            if isinstance(node, Within):
                within_node = node
                break

        assert within_node is not None

    def test_parse_intersects_any(self):
        """
        GIVEN a SQL query with INTERSECTS ANY operator
        WHEN parsing the query
        THEN should create a SpatialSetPredicate with ANY quantifier
        """
        sql = "SELECT * FROM v WHERE position INTERSECTS ANY('chr1:1000-2000', 'chr1:5000-6000')"
        ast = parse_one(sql, dialect=GIQLDialect)

        spatial_set = None
        for node in ast.walk():
            if isinstance(node, SpatialSetPredicate):
                spatial_set = node
                break

        assert spatial_set is not None
        assert spatial_set.args["operator"] == "INTERSECTS"
        assert spatial_set.args["quantifier"] == "ANY"

    def test_parse_intersects_all(self):
        """
        GIVEN a SQL query with INTERSECTS ALL operator
        WHEN parsing the query
        THEN should create a SpatialSetPredicate with ALL quantifier
        """
        sql = "SELECT * FROM v WHERE position INTERSECTS ALL('chr1:1000-2000', 'chr1:1500-1800')"
        ast = parse_one(sql, dialect=GIQLDialect)

        spatial_set = None
        for node in ast.walk():
            if isinstance(node, SpatialSetPredicate):
                spatial_set = node
                break

        assert spatial_set is not None
        assert spatial_set.args["operator"] == "INTERSECTS"
        assert spatial_set.args["quantifier"] == "ALL"

    def test_parse_contains_any(self):
        """
        GIVEN a SQL query with CONTAINS ANY operator
        WHEN parsing the query
        THEN should create a SpatialSetPredicate with CONTAINS operator
        """
        sql = "SELECT * FROM v WHERE position CONTAINS ANY('chr1:1500', 'chr1:1600')"
        ast = parse_one(sql, dialect=GIQLDialect)

        spatial_set = None
        for node in ast.walk():
            if isinstance(node, SpatialSetPredicate):
                spatial_set = node
                break

        assert spatial_set is not None
        assert spatial_set.args["operator"] == "CONTAINS"
        assert spatial_set.args["quantifier"] == "ANY"
