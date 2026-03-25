"""Tests for giql.dialect module."""

from sqlglot import exp
from sqlglot import parse_one
from sqlglot.tokens import TokenType

from giql.dialect import CONTAINS
from giql.dialect import INTERSECTS
from giql.dialect import WITHIN
from giql.dialect import GIQLDialect
from giql.expressions import Contains
from giql.expressions import GIQLCluster
from giql.expressions import GIQLCoverage
from giql.expressions import GIQLDistance
from giql.expressions import GIQLMerge
from giql.expressions import GIQLNearest
from giql.expressions import Intersects
from giql.expressions import SpatialPredicate
from giql.expressions import SpatialSetPredicate
from giql.expressions import Within


class TestDialectConstants:
    """Tests for module-level constants and token registration."""

    def test_dc_001_constant_values(self):
        """GIVEN the module is imported
        WHEN INTERSECTS, CONTAINS, WITHIN constants are accessed
        THEN they equal "INTERSECTS", "CONTAINS", "WITHIN" respectively.
        """
        assert INTERSECTS == "INTERSECTS"
        assert CONTAINS == "CONTAINS"
        assert WITHIN == "WITHIN"

    def test_dc_002_token_type_attributes(self):
        """GIVEN the module is imported
        WHEN TokenType attributes are checked
        THEN TokenType has INTERSECTS, CONTAINS, WITHIN attributes.
        """
        assert hasattr(TokenType, "INTERSECTS")
        assert hasattr(TokenType, "CONTAINS")
        assert hasattr(TokenType, "WITHIN")


class TestGIQLDialect:
    """Tests for GIQLDialect parsing of spatial predicates and GIQL functions."""

    def test_gd_001_intersects_predicate(self):
        """GIVEN a query string with `column INTERSECTS 'chr1:1000-2000'`
        WHEN the query is parsed with GIQLDialect
        THEN the AST contains an Intersects node with correct left and right expressions.
        """
        ast = parse_one(
            "SELECT * FROM t WHERE column INTERSECTS 'chr1:1000-2000'",
            dialect=GIQLDialect,
        )
        nodes = list(ast.find_all(Intersects))
        assert len(nodes) == 1
        node = nodes[0]
        assert node.this.name == "column"
        assert node.expression.this == "chr1:1000-2000"

    def test_gd_002_contains_predicate(self):
        """GIVEN a query string with `column CONTAINS 'chr1:1500'`
        WHEN the query is parsed with GIQLDialect
        THEN the AST contains a Contains node.
        """
        ast = parse_one(
            "SELECT * FROM t WHERE column CONTAINS 'chr1:1500'",
            dialect=GIQLDialect,
        )
        nodes = list(ast.find_all(Contains))
        assert len(nodes) == 1

    def test_gd_003_within_predicate(self):
        """GIVEN a query string with `column WITHIN 'chr1:1000-5000'`
        WHEN the query is parsed with GIQLDialect
        THEN the AST contains a Within node.
        """
        ast = parse_one(
            "SELECT * FROM t WHERE column WITHIN 'chr1:1000-5000'",
            dialect=GIQLDialect,
        )
        nodes = list(ast.find_all(Within))
        assert len(nodes) == 1

    def test_gd_004_intersects_any(self):
        """GIVEN a query with `column INTERSECTS ANY('chr1:1000-2000', 'chr1:5000-6000')`
        WHEN the query is parsed
        THEN the AST contains a SpatialSetPredicate with quantifier=ANY.
        """
        ast = parse_one(
            "SELECT * FROM t WHERE column INTERSECTS ANY('chr1:1000-2000', 'chr1:5000-6000')",
            dialect=GIQLDialect,
        )
        nodes = list(ast.find_all(SpatialSetPredicate))
        assert len(nodes) == 1
        node = nodes[0]
        assert node.args["quantifier"] == "ANY"

    def test_gd_005_intersects_all(self):
        """GIVEN a query with `column INTERSECTS ALL('chr1:1000-2000', 'chr1:5000-6000')`
        WHEN the query is parsed
        THEN the AST contains a SpatialSetPredicate with quantifier=ALL.
        """
        ast = parse_one(
            "SELECT * FROM t WHERE column INTERSECTS ALL('chr1:1000-2000', 'chr1:5000-6000')",
            dialect=GIQLDialect,
        )
        nodes = list(ast.find_all(SpatialSetPredicate))
        assert len(nodes) == 1
        node = nodes[0]
        assert node.args["quantifier"] == "ALL"

    def test_gd_006_plain_sql_fallback(self):
        """GIVEN a query with no spatial operators (plain SQL)
        WHEN the query is parsed with GIQLDialect
        THEN the AST is a standard SELECT without spatial nodes.
        """
        ast = parse_one(
            "SELECT id, name FROM t WHERE id = 1",
            dialect=GIQLDialect,
        )
        spatial_nodes = list(ast.find_all(SpatialPredicate, SpatialSetPredicate))
        assert len(spatial_nodes) == 0
        assert ast.find(exp.Select) is not None

    def test_gd_007_cluster_basic(self):
        """GIVEN a query with `CLUSTER(interval)`
        WHEN the query is parsed
        THEN the AST contains a GIQLCluster node.
        """
        ast = parse_one(
            "SELECT CLUSTER(interval) FROM t",
            dialect=GIQLDialect,
        )
        nodes = list(ast.find_all(GIQLCluster))
        assert len(nodes) == 1

    def test_gd_008_cluster_with_distance(self):
        """GIVEN a query with `CLUSTER(interval, 1000)`
        WHEN the query is parsed
        THEN the GIQLCluster node has distance arg set.
        """
        ast = parse_one(
            "SELECT CLUSTER(interval, 1000) FROM t",
            dialect=GIQLDialect,
        )
        nodes = list(ast.find_all(GIQLCluster))
        assert len(nodes) == 1
        node = nodes[0]
        assert node.args.get("distance") is not None

    def test_gd_009_merge_basic(self):
        """GIVEN a query with `MERGE(interval)`
        WHEN the query is parsed
        THEN the AST contains a GIQLMerge node.
        """
        ast = parse_one(
            "SELECT MERGE(interval) FROM t",
            dialect=GIQLDialect,
        )
        nodes = list(ast.find_all(GIQLMerge))
        assert len(nodes) == 1

    def test_gd_010_coverage_with_resolution(self):
        """GIVEN a query with `COVERAGE(interval, 1000)`
        WHEN the query is parsed
        THEN the AST contains a GIQLCoverage node with resolution set.
        """
        ast = parse_one(
            "SELECT COVERAGE(interval, 1000) FROM t",
            dialect=GIQLDialect,
        )
        nodes = list(ast.find_all(GIQLCoverage))
        assert len(nodes) == 1
        node = nodes[0]
        assert node.args.get("resolution") is not None

    def test_gd_011_coverage_with_stat(self):
        """GIVEN a query with `COVERAGE(interval, 500, stat := 'mean')`
        WHEN the query is parsed
        THEN the GIQLCoverage node has stat arg set.
        """
        ast = parse_one(
            "SELECT COVERAGE(interval, 500, stat := 'mean') FROM t",
            dialect=GIQLDialect,
        )
        nodes = list(ast.find_all(GIQLCoverage))
        assert len(nodes) == 1
        node = nodes[0]
        assert node.args.get("stat") is not None
        assert node.args["stat"].this == "mean"

    def test_gd_012_coverage_with_kwarg_resolution(self):
        """GIVEN a query with `COVERAGE(interval, resolution => 1000)`
        WHEN the query is parsed
        THEN the GIQLCoverage node has resolution set via Kwarg.
        """
        ast = parse_one(
            "SELECT COVERAGE(interval, resolution => 1000) FROM t",
            dialect=GIQLDialect,
        )
        nodes = list(ast.find_all(GIQLCoverage))
        assert len(nodes) == 1
        node = nodes[0]
        assert node.args.get("resolution") is not None

    def test_gd_013_coverage_with_stat_and_target(self):
        """GIVEN a query with `COVERAGE(interval, 1000, stat := 'mean', target := 'score')`
        WHEN the query is parsed
        THEN the GIQLCoverage node has stat and target args set.
        """
        ast = parse_one(
            "SELECT COVERAGE(interval, 1000, stat := 'mean', target := 'score') FROM t",
            dialect=GIQLDialect,
        )
        nodes = list(ast.find_all(GIQLCoverage))
        assert len(nodes) == 1
        node = nodes[0]
        assert node.args.get("stat") is not None
        assert node.args["stat"].this == "mean"
        assert node.args.get("target") is not None
        assert node.args["target"].this == "score"

    def test_gd_014_distance_function(self):
        """GIVEN a query with `DISTANCE(a.interval, b.interval)`
        WHEN the query is parsed
        THEN the AST contains a GIQLDistance node.
        """
        ast = parse_one(
            "SELECT DISTANCE(a.interval, b.interval) FROM t",
            dialect=GIQLDialect,
        )
        nodes = list(ast.find_all(GIQLDistance))
        assert len(nodes) == 1

    def test_gd_015_nearest_with_k(self):
        """GIVEN a query with `NEAREST(genes, k=3)`
        WHEN the query is parsed
        THEN the AST contains a GIQLNearest node with k arg set.
        """
        ast = parse_one(
            "SELECT NEAREST(genes, k=3) FROM t",
            dialect=GIQLDialect,
        )
        nodes = list(ast.find_all(GIQLNearest))
        assert len(nodes) == 1
        node = nodes[0]
        assert node.args.get("k") is not None
