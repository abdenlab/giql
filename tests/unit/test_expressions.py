"""Tests for custom AST expression nodes.

Test specification: specs/test_expressions.md
"""

from hypothesis import HealthCheck
from hypothesis import given
from hypothesis import settings
from hypothesis import strategies as st
from sqlglot import exp
from sqlglot import parse_one

from giql.dialect import GIQLDialect
from giql.expressions import Contains
from giql.expressions import GenomicRange
from giql.expressions import GIQLCluster
from giql.expressions import GIQLCoverage
from giql.expressions import GIQLDistance
from giql.expressions import GIQLMerge
from giql.expressions import GIQLNearest
from giql.expressions import Intersects
from giql.expressions import SpatialPredicate
from giql.expressions import SpatialSetPredicate
from giql.expressions import Within

VALID_STATS = ["count", "mean", "sum", "min", "max"]


class TestGenomicRange:
    """Tests for GenomicRange expression node."""

    def test___init___should_succeed_when_required_args_supplied(self):
        """Test GenomicRange instantiates with just required args.

        Given:
            All required args (chromosome, start, end)
        When:
            GenomicRange is instantiated
        Then:
            It should have correct chromosome, start, and end args
        """
        # Arrange
        chrom = exp.Literal.string("chr1")
        start = exp.Literal.number(1000)
        end = exp.Literal.number(2000)

        # Act
        gr = GenomicRange(chromosome=chrom, start=start, end=end)

        # Assert
        assert gr.args["chromosome"] is chrom
        assert gr.args["start"] is start
        assert gr.args["end"] is end

    def test___init___should_accept_all_args_when_optional_supplied(self):
        """Test GenomicRange instantiates with all optional args.

        Given:
            Required args plus optional strand and coord_system
        When:
            GenomicRange is instantiated
        Then:
            It should have all five args accessible
        """
        # Arrange
        chrom = exp.Literal.string("chr1")
        start = exp.Literal.number(1000)
        end = exp.Literal.number(2000)
        strand = exp.Literal.string("+")
        coord_system = exp.Literal.string("0-based")

        # Act
        gr = GenomicRange(
            chromosome=chrom,
            start=start,
            end=end,
            strand=strand,
            coord_system=coord_system,
        )

        # Assert
        assert gr.args["chromosome"] is chrom
        assert gr.args["start"] is start
        assert gr.args["end"] is end
        assert gr.args["strand"] is strand
        assert gr.args["coord_system"] is coord_system

    def test___init___should_default_optional_args_to_none_when_omitted(self):
        """Test GenomicRange defaults optional args to None.

        Given:
            Only required args provided
        When:
            GenomicRange is instantiated
        Then:
            It should leave strand and coord_system args as None
        """
        # Act
        gr = GenomicRange(
            chromosome=exp.Literal.string("chr1"),
            start=exp.Literal.number(1000),
            end=exp.Literal.number(2000),
        )

        # Assert
        assert gr.args.get("strand") is None
        assert gr.args.get("coord_system") is None


class TestSpatialPredicate:
    """Tests for SpatialPredicate subclasses."""

    def test___init___should_produce_spatial_predicate_and_binary_when_intersects(self):
        """Test Intersects inherits from SpatialPredicate and exp.Binary.

        Given:
            Two expression nodes (this, expression)
        When:
            Intersects is instantiated
        Then:
            It should produce an instance of SpatialPredicate and exp.Binary
        """
        # Arrange
        left = exp.Column(this=exp.Identifier(this="a"))
        right = exp.Column(this=exp.Identifier(this="b"))

        # Act
        node = Intersects(this=left, expression=right)

        # Assert
        assert isinstance(node, SpatialPredicate)
        assert isinstance(node, exp.Binary)

    def test___init___should_produce_spatial_predicate_and_binary_when_contains(self):
        """Test Contains inherits from SpatialPredicate and exp.Binary.

        Given:
            Two expression nodes
        When:
            Contains is instantiated
        Then:
            It should produce an instance of SpatialPredicate and exp.Binary
        """
        # Arrange
        left = exp.Column(this=exp.Identifier(this="a"))
        right = exp.Column(this=exp.Identifier(this="b"))

        # Act
        node = Contains(this=left, expression=right)

        # Assert
        assert isinstance(node, SpatialPredicate)
        assert isinstance(node, exp.Binary)

    def test___init___should_produce_spatial_predicate_and_binary_when_within(self):
        """Test Within inherits from SpatialPredicate and exp.Binary.

        Given:
            Two expression nodes
        When:
            Within is instantiated
        Then:
            It should produce an instance of SpatialPredicate and exp.Binary
        """
        # Arrange
        left = exp.Column(this=exp.Identifier(this="a"))
        right = exp.Column(this=exp.Identifier(this="b"))

        # Act
        node = Within(this=left, expression=right)

        # Assert
        assert isinstance(node, SpatialPredicate)
        assert isinstance(node, exp.Binary)


class TestSpatialSetPredicate:
    """Tests for SpatialSetPredicate expression node."""

    def test___init___should_set_all_args_when_required_args_supplied(self):
        """Test SpatialSetPredicate instantiates with all required args.

        Given:
            All required args (this, operator, quantifier, ranges)
        When:
            SpatialSetPredicate is instantiated
        Then:
            It should have all four args accessible
        """
        # Arrange
        this = exp.Column(this=exp.Identifier(this="interval"))
        operator = exp.Literal.string("INTERSECTS")
        quantifier = exp.Literal.string("ANY")
        ranges = exp.Array(
            expressions=[
                exp.Literal.string("chr1:1000-2000"),
                exp.Literal.string("chr1:5000-6000"),
            ]
        )

        # Act
        node = SpatialSetPredicate(
            this=this,
            operator=operator,
            quantifier=quantifier,
            ranges=ranges,
        )

        # Assert
        assert node.args["this"] is this
        assert node.args["operator"] is operator
        assert node.args["quantifier"] is quantifier
        assert node.args["ranges"] is ranges


class TestGIQLCluster:
    """Tests for GIQLCluster expression node parsing."""

    def test_parse_should_set_this_when_one_positional_arg(self):
        """Test CLUSTER parses with a single positional arg.

        Given:
            A CLUSTER expression with one positional arg (column)
        When:
            Parsed with GIQLDialect
        Then:
            It should produce a GIQLCluster instance with `this` set
        """
        # Act
        ast = parse_one(
            "SELECT CLUSTER(interval) FROM features",
            dialect=GIQLDialect,
        )

        # Assert
        nodes = list(ast.find_all(GIQLCluster))
        assert len(nodes) == 1
        assert nodes[0].args["this"] is not None

    def test_parse_should_set_distance_when_two_positional_args(self):
        """Test CLUSTER parses with column and distance positionals.

        Given:
            A CLUSTER expression with two positional args (column, distance)
        When:
            Parsed with GIQLDialect
        Then:
            It should produce a GIQLCluster instance with `this` and `distance` set
        """
        # Act
        ast = parse_one(
            "SELECT CLUSTER(interval, 1000) FROM features",
            dialect=GIQLDialect,
        )

        # Assert
        nodes = list(ast.find_all(GIQLCluster))
        assert len(nodes) == 1
        assert nodes[0].args["this"] is not None
        assert nodes[0].args["distance"].this == "1000"

    def test_parse_should_set_stranded_when_named_parameter_supplied(self):
        """Test CLUSTER parses with a stranded named parameter.

        Given:
            A CLUSTER expression with one positional and stranded := true
        When:
            Parsed with GIQLDialect
        Then:
            It should produce a GIQLCluster instance with `this` and `stranded` set
        """
        # Act
        ast = parse_one(
            "SELECT CLUSTER(interval, stranded := true) FROM features",
            dialect=GIQLDialect,
        )

        # Assert
        nodes = list(ast.find_all(GIQLCluster))
        assert len(nodes) == 1
        assert nodes[0].args["this"] is not None
        assert nodes[0].args["stranded"] is not None

    def test_parse_should_set_distance_and_stranded_when_both_supplied(self):
        """Test CLUSTER parses with both distance and stranded params.

        Given:
            A CLUSTER expression with two positionals and stranded := true
        When:
            Parsed with GIQLDialect
        Then:
            It should produce a GIQLCluster instance with `this`, `distance`, and `stranded` set
        """
        # Act
        ast = parse_one(
            "SELECT CLUSTER(interval, 1000, stranded := true) FROM features",
            dialect=GIQLDialect,
        )

        # Assert
        nodes = list(ast.find_all(GIQLCluster))
        assert len(nodes) == 1
        assert nodes[0].args["this"] is not None
        assert nodes[0].args["distance"].this == "1000"
        assert nodes[0].args["stranded"] is not None

    def test___init___should_leave_optional_args_absent_when_only_this_supplied(self):
        """Test GIQLCluster direct instantiation with just `this`.

        Given:
            Required arg `this` only
        When:
            GIQLCluster is instantiated directly
        Then:
            It should set `this` and leave `distance` and `stranded` absent
        """
        # Arrange
        col = exp.Column(this=exp.Identifier(this="interval"))

        # Act
        node = GIQLCluster(this=col)

        # Assert
        assert node.args["this"] is col
        assert node.args.get("distance") is None
        assert node.args.get("stranded") is None


class TestGIQLMerge:
    """Tests for GIQLMerge expression node parsing."""

    def test_parse_should_set_this_when_one_positional_arg(self):
        """Test MERGE parses with a single positional arg.

        Given:
            A MERGE expression with one positional arg (column)
        When:
            Parsed with GIQLDialect
        Then:
            It should produce a GIQLMerge instance with `this` set
        """
        # Act
        ast = parse_one(
            "SELECT MERGE(interval) FROM features",
            dialect=GIQLDialect,
        )

        # Assert
        nodes = list(ast.find_all(GIQLMerge))
        assert len(nodes) == 1
        assert nodes[0].args["this"] is not None

    def test_parse_should_set_distance_when_two_positional_args(self):
        """Test MERGE parses with column and distance positionals.

        Given:
            A MERGE expression with two positional args (column, distance)
        When:
            Parsed with GIQLDialect
        Then:
            It should produce a GIQLMerge instance with `this` and `distance` set
        """
        # Act
        ast = parse_one(
            "SELECT MERGE(interval, 1000) FROM features",
            dialect=GIQLDialect,
        )

        # Assert
        nodes = list(ast.find_all(GIQLMerge))
        assert len(nodes) == 1
        assert nodes[0].args["this"] is not None
        assert nodes[0].args["distance"].this == "1000"

    def test_parse_should_set_stranded_when_named_parameter_supplied(self):
        """Test MERGE parses with a stranded named parameter.

        Given:
            A MERGE expression with one positional and stranded := true
        When:
            Parsed with GIQLDialect
        Then:
            It should produce a GIQLMerge instance with `this` and `stranded` set
        """
        # Act
        ast = parse_one(
            "SELECT MERGE(interval, stranded := true) FROM features",
            dialect=GIQLDialect,
        )

        # Assert
        nodes = list(ast.find_all(GIQLMerge))
        assert len(nodes) == 1
        assert nodes[0].args["this"] is not None
        assert nodes[0].args["stranded"] is not None

    def test_parse_should_set_distance_and_stranded_when_both_supplied(self):
        """Test MERGE parses with both distance and stranded params.

        Given:
            A MERGE expression with two positionals and stranded := true
        When:
            Parsed with GIQLDialect
        Then:
            It should produce a GIQLMerge instance with `this`, `distance`, and `stranded` set
        """
        # Act
        ast = parse_one(
            "SELECT MERGE(interval, 1000, stranded := true) FROM features",
            dialect=GIQLDialect,
        )

        # Assert
        nodes = list(ast.find_all(GIQLMerge))
        assert len(nodes) == 1
        assert nodes[0].args["this"] is not None
        assert nodes[0].args["distance"].this == "1000"
        assert nodes[0].args["stranded"] is not None


class TestGIQLCoverage:
    """Tests for GIQLCoverage expression node parsing."""

    # ------------------------------------------------------------------
    # Example-based parsing (COV-001 to COV-007)
    # ------------------------------------------------------------------

    def test_from_arg_list_should_map_resolution_when_positional(self):
        """Test positional interval and resolution mapping.

        Given:
            A COVERAGE expression with positional interval and resolution
        When:
            Parsed with GIQLDialect
        Then:
            It should produce a GIQLCoverage node with resolution set and
            stat/target both None
        """
        # Act
        ast = parse_one(
            "SELECT COVERAGE(interval, 1000) FROM features",
            dialect=GIQLDialect,
        )

        # Assert
        coverage = list(ast.find_all(GIQLCoverage))
        assert len(coverage) == 1
        assert coverage[0].args["resolution"].this == "1000"
        assert coverage[0].args.get("stat") is None
        assert coverage[0].args.get("target") is None

    def test_from_arg_list_should_set_stat_when_walrus_syntax(self):
        """Test named stat parameter via := syntax.

        Given:
            A COVERAGE expression with := named stat parameter
        When:
            Parsed with GIQLDialect
        Then:
            It should produce a GIQLCoverage node with stat set to the given value
        """
        # Act
        ast = parse_one(
            "SELECT COVERAGE(interval, 500, stat := 'mean') FROM features",
            dialect=GIQLDialect,
        )

        # Assert
        coverage = list(ast.find_all(GIQLCoverage))
        assert len(coverage) == 1
        assert coverage[0].args["stat"].this == "mean"

    def test_from_arg_list_should_set_stat_when_arrow_syntax(self):
        """Test named stat parameter via => syntax.

        Given:
            A COVERAGE expression with => named stat parameter
        When:
            Parsed with GIQLDialect
        Then:
            It should produce a GIQLCoverage node with stat set to the given value
        """
        # Act
        ast = parse_one(
            "SELECT COVERAGE(interval, 500, stat => 'mean') FROM features",
            dialect=GIQLDialect,
        )

        # Assert
        coverage = list(ast.find_all(GIQLCoverage))
        assert len(coverage) == 1
        assert coverage[0].args["stat"].this == "mean"

    def test_from_arg_list_should_set_resolution_when_named(self):
        """Test named resolution parameter.

        Given:
            A COVERAGE expression with named resolution parameter
        When:
            Parsed with GIQLDialect
        Then:
            It should produce a GIQLCoverage node with resolution set via named param
        """
        # Act
        ast = parse_one(
            "SELECT COVERAGE(interval, resolution := 1000) FROM features",
            dialect=GIQLDialect,
        )

        # Assert
        coverage = list(ast.find_all(GIQLCoverage))
        assert len(coverage) == 1
        assert coverage[0].args["resolution"].this == "1000"

    def test_from_arg_list_should_set_target_when_walrus_syntax(self):
        """Test target parameter via := syntax.

        Given:
            A COVERAGE expression with := named target parameter
        When:
            Parsed with GIQLDialect
        Then:
            It should produce a GIQLCoverage node with target set
        """
        # Act
        ast = parse_one(
            "SELECT COVERAGE(interval, 1000, target := 'score') FROM features",
            dialect=GIQLDialect,
        )

        # Assert
        coverage = list(ast.find_all(GIQLCoverage))
        assert len(coverage) == 1
        assert coverage[0].args["target"].this == "score"

    def test_from_arg_list_should_set_target_when_arrow_syntax(self):
        """Test target parameter via => syntax.

        Given:
            A COVERAGE expression with => named target parameter
        When:
            Parsed with GIQLDialect
        Then:
            It should produce a GIQLCoverage node with target set
        """
        # Act
        ast = parse_one(
            "SELECT COVERAGE(interval, 1000, target => 'score') FROM features",
            dialect=GIQLDialect,
        )

        # Assert
        coverage = list(ast.find_all(GIQLCoverage))
        assert len(coverage) == 1
        assert coverage[0].args["target"].this == "score"

    def test_from_arg_list_should_set_all_params_when_all_named(self):
        """Test all parameters provided as named arguments.

        Given:
            A COVERAGE expression with stat, target, and resolution all named
        When:
            Parsed with GIQLDialect
        Then:
            It should produce a GIQLCoverage node with all three params set
        """
        # Act
        ast = parse_one(
            "SELECT COVERAGE(interval, resolution := 500, "
            "stat := 'mean', target := 'score') FROM features",
            dialect=GIQLDialect,
        )

        # Assert
        coverage = list(ast.find_all(GIQLCoverage))
        assert len(coverage) == 1
        assert coverage[0].args["resolution"].this == "500"
        assert coverage[0].args["stat"].this == "mean"
        assert coverage[0].args["target"].this == "score"

    # ------------------------------------------------------------------
    # Property-based parsing (PBT-001 to PBT-003)
    # ------------------------------------------------------------------

    @given(
        resolution=st.integers(min_value=1, max_value=10_000_000),
        stat=st.sampled_from(VALID_STATS),
        syntax=st.sampled_from([":=", "=>"]),
    )
    @settings(max_examples=50, suppress_health_check=[HealthCheck.function_scoped_fixture])
    def test_from_arg_list_should_parse_stat_and_resolution_when_varying_inputs(
        self, resolution, stat, syntax
    ):
        """Test stat and resolution parse correctly across input space.

        Given:
            Any valid resolution (1-10M), stat (sampled from valid values),
            and syntax (:= or =>)
        When:
            Parsed with GIQLDialect
        Then:
            It should produce a GIQLCoverage node with correct resolution and stat
        """
        # Act
        sql = (
            f"SELECT COVERAGE(interval, {resolution}, "
            f"stat {syntax} '{stat}') FROM features"
        )
        ast = parse_one(sql, dialect=GIQLDialect)

        # Assert
        coverage = list(ast.find_all(GIQLCoverage))
        assert len(coverage) == 1
        assert coverage[0].args["resolution"].this == str(resolution)
        assert coverage[0].args["stat"].this == stat

    @given(resolution=st.integers(min_value=1, max_value=10_000_000))
    @settings(max_examples=50, suppress_health_check=[HealthCheck.function_scoped_fixture])
    def test_from_arg_list_should_set_resolution_when_positional_only(self, resolution):
        """Test positional-only parsing across resolution range.

        Given:
            Any valid resolution (1-10M) with no stat or target
        When:
            Parsed with GIQLDialect
        Then:
            It should produce a GIQLCoverage node with resolution set and
            stat/target None
        """
        # Act
        ast = parse_one(
            f"SELECT COVERAGE(interval, {resolution}) FROM features",
            dialect=GIQLDialect,
        )

        # Assert
        coverage = list(ast.find_all(GIQLCoverage))
        assert len(coverage) == 1
        assert coverage[0].args["resolution"].this == str(resolution)
        assert coverage[0].args.get("stat") is None
        assert coverage[0].args.get("target") is None

    @given(syntax=st.sampled_from([":=", "=>"]))
    @settings(max_examples=50, suppress_health_check=[HealthCheck.function_scoped_fixture])
    def test_from_arg_list_should_set_target_when_varying_syntax(self, syntax):
        """Test target parameter parsing across syntax variants.

        Given:
            Either := or => syntax for target parameter
        When:
            Parsed with GIQLDialect
        Then:
            It should produce a GIQLCoverage node with target set
        """
        # Act
        ast = parse_one(
            f"SELECT COVERAGE(interval, 1000, target {syntax} 'score') FROM features",
            dialect=GIQLDialect,
        )

        # Assert
        coverage = list(ast.find_all(GIQLCoverage))
        assert len(coverage) == 1
        assert coverage[0].args["target"].this == "score"


class TestGIQLDistance:
    """Tests for GIQLDistance expression node parsing."""

    def test_parse_should_set_this_and_expression_when_two_positional_args(self):
        """Test DISTANCE parses with two positional interval args.

        Given:
            A DISTANCE expression with two positional args (interval_a, interval_b)
        When:
            Parsed with GIQLDialect
        Then:
            It should produce a GIQLDistance instance with `this` and `expression` set
        """
        # Act
        ast = parse_one(
            "SELECT DISTANCE(a.interval, b.interval) FROM a, b",
            dialect=GIQLDialect,
        )

        # Assert
        nodes = list(ast.find_all(GIQLDistance))
        assert len(nodes) == 1
        assert nodes[0].args["this"] is not None
        assert nodes[0].args["expression"] is not None

    def test_parse_should_set_stranded_and_signed_when_both_named_params(self):
        """Test DISTANCE parses with stranded and signed named params.

        Given:
            A DISTANCE expression with two positionals and stranded := true, signed := true
        When:
            Parsed with GIQLDialect
        Then:
            It should produce a GIQLDistance instance with `this`, `expression`, `stranded`, and `signed` set
        """
        # Act
        ast = parse_one(
            "SELECT DISTANCE(a.interval, b.interval, stranded := true, signed := true) FROM a, b",
            dialect=GIQLDialect,
        )

        # Assert
        nodes = list(ast.find_all(GIQLDistance))
        assert len(nodes) == 1
        assert nodes[0].args["this"] is not None
        assert nodes[0].args["expression"] is not None
        assert nodes[0].args["stranded"] is not None
        assert nodes[0].args["signed"] is not None

    def test_parse_should_leave_signed_absent_when_only_stranded_supplied(self):
        """Test DISTANCE parses with only stranded named param.

        Given:
            A DISTANCE expression with two positionals and only stranded := true
        When:
            Parsed with GIQLDialect
        Then:
            It should produce a GIQLDistance instance with `this`, `expression`, and `stranded` set; `signed` absent
        """
        # Act
        ast = parse_one(
            "SELECT DISTANCE(a.interval, b.interval, stranded := true) FROM a, b",
            dialect=GIQLDialect,
        )

        # Assert
        nodes = list(ast.find_all(GIQLDistance))
        assert len(nodes) == 1
        assert nodes[0].args["this"] is not None
        assert nodes[0].args["expression"] is not None
        assert nodes[0].args["stranded"] is not None
        assert nodes[0].args.get("signed") is None


class TestGIQLNearest:
    """Tests for GIQLNearest expression node parsing."""

    def test_parse_should_set_this_when_one_positional_arg(self):
        """Test NEAREST parses with a single positional table arg.

        Given:
            A NEAREST expression with one positional arg (table)
        When:
            Parsed with GIQLDialect
        Then:
            It should produce a GIQLNearest instance with `this` set
        """
        # Act
        ast = parse_one(
            "SELECT NEAREST(genes) FROM peaks",
            dialect=GIQLDialect,
        )

        # Assert
        nodes = list(ast.find_all(GIQLNearest))
        assert len(nodes) == 1
        assert nodes[0].args["this"] is not None

    def test_parse_should_set_k_when_named_parameter_supplied(self):
        """Test NEAREST parses with a k named parameter.

        Given:
            A NEAREST expression with one positional and k := 3
        When:
            Parsed with GIQLDialect
        Then:
            It should produce a GIQLNearest instance with `this` and `k` set
        """
        # Act
        ast = parse_one(
            "SELECT NEAREST(genes, k := 3) FROM peaks",
            dialect=GIQLDialect,
        )

        # Assert
        nodes = list(ast.find_all(GIQLNearest))
        assert len(nodes) == 1
        assert nodes[0].args["this"] is not None
        assert nodes[0].args["k"].this == "3"

    def test_parse_should_set_all_args_when_multiple_named_params(self):
        """Test NEAREST parses with multiple named params.

        Given:
            A NEAREST expression with one positional and multiple named params
        When:
            Parsed with GIQLDialect
        Then:
            It should produce a GIQLNearest instance with all provided args set
        """
        # Act
        ast = parse_one(
            "SELECT NEAREST(genes, k := 5, max_distance := 100000, stranded := true, signed := true) FROM peaks",
            dialect=GIQLDialect,
        )

        # Assert
        nodes = list(ast.find_all(GIQLNearest))
        assert len(nodes) == 1
        assert nodes[0].args["this"] is not None
        assert nodes[0].args["k"].this == "5"
        assert nodes[0].args["max_distance"].this == "100000"
        assert nodes[0].args["stranded"] is not None
        assert nodes[0].args["signed"] is not None
