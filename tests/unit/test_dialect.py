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

    def test_constants_should_equal_their_uppercase_names(self):
        """Test module-level spatial-operator constants expose their uppercase names.

        Given:
            The giql.dialect module is imported
        When:
            INTERSECTS, CONTAINS, and WITHIN constants are accessed
        Then:
            It should equal "INTERSECTS", "CONTAINS", and "WITHIN" respectively
        """
        # DC-001
        # Arrange / Act / Assert
        assert INTERSECTS == "INTERSECTS"
        assert CONTAINS == "CONTAINS"
        assert WITHIN == "WITHIN"

    def test_TokenType_should_expose_spatial_operator_attributes(self):
        """Test that TokenType is extended with spatial-operator attributes.

        Given:
            The giql.dialect module is imported
        When:
            TokenType attributes are checked for spatial operators
        Then:
            It should expose INTERSECTS, CONTAINS, and WITHIN attributes
        """
        # DC-002
        # Arrange / Act / Assert
        assert hasattr(TokenType, "INTERSECTS")
        assert hasattr(TokenType, "CONTAINS")
        assert hasattr(TokenType, "WITHIN")


class TestGIQLDialect:
    """Tests for GIQLDialect parsing of spatial predicates and GIQL functions."""

    def test_parse_one_should_produce_Intersects_node_for_intersects_predicate(self):
        """Test parsing `column INTERSECTS 'chr1:1000-2000'` yields an Intersects node.

        Given:
            A SELECT query containing `column INTERSECTS 'chr1:1000-2000'`
        When:
            The query is parsed with GIQLDialect
        Then:
            It should produce an Intersects node whose left side is the column
            and whose right side is the literal range string
        """
        # GD-001
        # Arrange
        query = "SELECT * FROM t WHERE column INTERSECTS 'chr1:1000-2000'"

        # Act
        ast = parse_one(
            query,
            dialect=GIQLDialect,
        )

        # Assert
        nodes = list(ast.find_all(Intersects))
        assert len(nodes) == 1
        node = nodes[0]
        assert node.this.name == "column"
        assert node.expression.this == "chr1:1000-2000"

    def test_parse_one_should_produce_Contains_node_for_contains_predicate(self):
        """Test parsing `column CONTAINS 'chr1:1500'` yields a Contains node.

        Given:
            A SELECT query containing `column CONTAINS 'chr1:1500'`
        When:
            The query is parsed with GIQLDialect
        Then:
            It should produce exactly one Contains node in the AST
        """
        # GD-002
        # Arrange
        query = "SELECT * FROM t WHERE column CONTAINS 'chr1:1500'"

        # Act
        ast = parse_one(
            query,
            dialect=GIQLDialect,
        )

        # Assert
        nodes = list(ast.find_all(Contains))
        assert len(nodes) == 1

    def test_parse_one_should_produce_Within_node_for_within_predicate(self):
        """Test parsing `column WITHIN 'chr1:1000-5000'` yields a Within node.

        Given:
            A SELECT query containing `column WITHIN 'chr1:1000-5000'`
        When:
            The query is parsed with GIQLDialect
        Then:
            It should produce exactly one Within node in the AST
        """
        # GD-003
        # Arrange
        query = "SELECT * FROM t WHERE column WITHIN 'chr1:1000-5000'"

        # Act
        ast = parse_one(
            query,
            dialect=GIQLDialect,
        )

        # Assert
        nodes = list(ast.find_all(Within))
        assert len(nodes) == 1

    def test_parse_one_should_set_quantifier_to_ANY_for_intersects_any(self):
        """Test `INTERSECTS ANY(...)` produces a SpatialSetPredicate with quantifier ANY.

        Given:
            A SELECT query containing `column INTERSECTS ANY('chr1:1000-2000', 'chr1:5000-6000')`
        When:
            The query is parsed with GIQLDialect
        Then:
            It should produce a SpatialSetPredicate whose quantifier argument is "ANY"
        """
        # GD-004
        # Arrange
        query = (
            "SELECT * FROM t WHERE column INTERSECTS ANY('chr1:1000-2000', 'chr1:5000-6000')"
        )

        # Act
        ast = parse_one(
            query,
            dialect=GIQLDialect,
        )

        # Assert
        nodes = list(ast.find_all(SpatialSetPredicate))
        assert len(nodes) == 1
        node = nodes[0]
        assert node.args["quantifier"] == "ANY"

    def test_parse_one_should_set_quantifier_to_ALL_for_intersects_all(self):
        """Test `INTERSECTS ALL(...)` produces a SpatialSetPredicate with quantifier ALL.

        Given:
            A SELECT query containing `column INTERSECTS ALL('chr1:1000-2000', 'chr1:5000-6000')`
        When:
            The query is parsed with GIQLDialect
        Then:
            It should produce a SpatialSetPredicate whose quantifier argument is "ALL"
        """
        # GD-005
        # Arrange
        query = (
            "SELECT * FROM t WHERE column INTERSECTS ALL('chr1:1000-2000', 'chr1:5000-6000')"
        )

        # Act
        ast = parse_one(
            query,
            dialect=GIQLDialect,
        )

        # Assert
        nodes = list(ast.find_all(SpatialSetPredicate))
        assert len(nodes) == 1
        node = nodes[0]
        assert node.args["quantifier"] == "ALL"

    def test_parse_one_should_produce_plain_select_when_no_spatial_operators_are_used(self):
        """Test plain SQL parses without any spatial nodes under GIQLDialect.

        Given:
            A SELECT query with no spatial operators
        When:
            The query is parsed with GIQLDialect
        Then:
            It should produce a standard Select AST containing no
            SpatialPredicate or SpatialSetPredicate nodes
        """
        # GD-006
        # Arrange
        query = "SELECT id, name FROM t WHERE id = 1"

        # Act
        ast = parse_one(
            query,
            dialect=GIQLDialect,
        )

        # Assert
        spatial_nodes = list(ast.find_all(SpatialPredicate, SpatialSetPredicate))
        assert len(spatial_nodes) == 0
        assert ast.find(exp.Select) is not None

    def test_parse_one_should_produce_GIQLCluster_node_for_cluster_call(self):
        """Test `CLUSTER(interval)` parses into a GIQLCluster AST node.

        Given:
            A SELECT query containing `CLUSTER(interval)`
        When:
            The query is parsed with GIQLDialect
        Then:
            It should produce exactly one GIQLCluster node in the AST
        """
        # GD-007
        # Arrange
        query = "SELECT CLUSTER(interval) FROM t"

        # Act
        ast = parse_one(
            query,
            dialect=GIQLDialect,
        )

        # Assert
        nodes = list(ast.find_all(GIQLCluster))
        assert len(nodes) == 1

    def test_parse_one_should_set_distance_arg_on_GIQLCluster_when_distance_is_given(self):
        """Test `CLUSTER(interval, 1000)` sets the distance argument on GIQLCluster.

        Given:
            A SELECT query containing `CLUSTER(interval, 1000)`
        When:
            The query is parsed with GIQLDialect
        Then:
            It should produce a GIQLCluster node whose distance argument is set
        """
        # GD-008
        # Arrange
        query = "SELECT CLUSTER(interval, 1000) FROM t"

        # Act
        ast = parse_one(
            query,
            dialect=GIQLDialect,
        )

        # Assert
        nodes = list(ast.find_all(GIQLCluster))
        assert len(nodes) == 1
        node = nodes[0]
        assert node.args.get("distance") is not None

    def test_parse_one_should_produce_GIQLMerge_node_for_merge_call(self):
        """Test `MERGE(interval)` parses into a GIQLMerge AST node.

        Given:
            A SELECT query containing `MERGE(interval)`
        When:
            The query is parsed with GIQLDialect
        Then:
            It should produce exactly one GIQLMerge node in the AST
        """
        # GD-009
        # Arrange
        query = "SELECT MERGE(interval) FROM t"

        # Act
        ast = parse_one(
            query,
            dialect=GIQLDialect,
        )

        # Assert
        nodes = list(ast.find_all(GIQLMerge))
        assert len(nodes) == 1

    def test_parse_one_should_set_resolution_arg_on_GIQLCoverage_when_resolution_is_positional(self):
        """Test `COVERAGE(interval, 1000)` sets the resolution argument on GIQLCoverage.

        Given:
            A SELECT query containing `COVERAGE(interval, 1000)`
        When:
            The query is parsed with GIQLDialect
        Then:
            It should produce a GIQLCoverage node whose resolution argument is set
        """
        # GD-010
        # Arrange
        query = "SELECT COVERAGE(interval, 1000) FROM t"

        # Act
        ast = parse_one(
            query,
            dialect=GIQLDialect,
        )

        # Assert
        nodes = list(ast.find_all(GIQLCoverage))
        assert len(nodes) == 1
        node = nodes[0]
        assert node.args.get("resolution") is not None

    def test_parse_one_should_set_resolution_arg_on_GIQLCoverage_when_resolution_is_passed_as_kwarg(self):
        """Test `COVERAGE(interval, resolution => 1000)` sets resolution via Kwarg syntax.

        Given:
            A SELECT query containing `COVERAGE(interval, resolution => 1000)`
        When:
            The query is parsed with GIQLDialect
        Then:
            It should produce a GIQLCoverage node whose resolution argument is set
        """
        # GD-012
        # Arrange
        query = "SELECT COVERAGE(interval, resolution => 1000) FROM t"

        # Act
        ast = parse_one(
            query,
            dialect=GIQLDialect,
        )

        # Assert
        nodes = list(ast.find_all(GIQLCoverage))
        assert len(nodes) == 1
        node = nodes[0]
        assert node.args.get("resolution") is not None

    def test_parse_one_should_produce_GIQLDistance_node_for_distance_call(self):
        """Test `DISTANCE(a.interval, b.interval)` parses into a GIQLDistance AST node.

        Given:
            A SELECT query containing `DISTANCE(a.interval, b.interval)`
        When:
            The query is parsed with GIQLDialect
        Then:
            It should produce exactly one GIQLDistance node in the AST
        """
        # GD-014
        # Arrange
        query = "SELECT DISTANCE(a.interval, b.interval) FROM t"

        # Act
        ast = parse_one(
            query,
            dialect=GIQLDialect,
        )

        # Assert
        nodes = list(ast.find_all(GIQLDistance))
        assert len(nodes) == 1

    def test_parse_one_should_set_k_arg_on_GIQLNearest_when_k_named_param_is_given(self):
        """Test `NEAREST(genes, k := 3)` sets the k argument on GIQLNearest.

        Given:
            A SELECT query containing `NEAREST(genes, k := 3)`
        When:
            The query is parsed with GIQLDialect
        Then:
            It should produce a GIQLNearest node whose k argument is set
        """
        # GD-015
        # Arrange
        query = "SELECT NEAREST(genes, k := 3) FROM t"

        # Act
        ast = parse_one(
            query,
            dialect=GIQLDialect,
        )

        # Assert
        nodes = list(ast.find_all(GIQLNearest))
        assert len(nodes) == 1
        node = nodes[0]
        assert node.args.get("k") is not None
