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

    def test_instantiate_with_required_args(self):
        """GR-001: Instantiate with required args.

        Given:
            All required args (chromosome, start, end)
        When:
            GenomicRange is instantiated
        Then:
            Instance has correct chromosome, start, and end args
        """
        chrom = exp.Literal.string("chr1")
        start = exp.Literal.number(1000)
        end = exp.Literal.number(2000)

        gr = GenomicRange(chromosome=chrom, start=start, end=end)

        assert gr.args["chromosome"] is chrom
        assert gr.args["start"] is start
        assert gr.args["end"] is end

    def test_instantiate_with_all_args(self):
        """GR-002: Instantiate with all args including optional strand and coord_system.

        Given:
            Required args plus optional strand and coord_system
        When:
            GenomicRange is instantiated
        Then:
            Instance has all five args accessible
        """
        chrom = exp.Literal.string("chr1")
        start = exp.Literal.number(1000)
        end = exp.Literal.number(2000)
        strand = exp.Literal.string("+")
        coord_system = exp.Literal.string("0-based")

        gr = GenomicRange(
            chromosome=chrom,
            start=start,
            end=end,
            strand=strand,
            coord_system=coord_system,
        )

        assert gr.args["chromosome"] is chrom
        assert gr.args["start"] is start
        assert gr.args["end"] is end
        assert gr.args["strand"] is strand
        assert gr.args["coord_system"] is coord_system

    def test_optional_args_default_to_none(self):
        """GR-003: Optional args default to None.

        Given:
            Only required args provided
        When:
            GenomicRange is instantiated
        Then:
            strand and coord_system args are None
        """
        gr = GenomicRange(
            chromosome=exp.Literal.string("chr1"),
            start=exp.Literal.number(1000),
            end=exp.Literal.number(2000),
        )

        assert gr.args.get("strand") is None
        assert gr.args.get("coord_system") is None


class TestSpatialPredicate:
    """Tests for SpatialPredicate subclasses."""

    def test_intersects_is_spatial_predicate_and_binary(self):
        """SP-001: Intersects inheritance.

        Given:
            Two expression nodes (this, expression)
        When:
            Intersects is instantiated
        Then:
            Instance is a SpatialPredicate and exp.Binary
        """
        left = exp.Column(this=exp.Identifier(this="a"))
        right = exp.Column(this=exp.Identifier(this="b"))

        node = Intersects(this=left, expression=right)

        assert isinstance(node, SpatialPredicate)
        assert isinstance(node, exp.Binary)

    def test_contains_is_spatial_predicate_and_binary(self):
        """SP-002: Contains inheritance.

        Given:
            Two expression nodes
        When:
            Contains is instantiated
        Then:
            Instance is a SpatialPredicate and exp.Binary
        """
        left = exp.Column(this=exp.Identifier(this="a"))
        right = exp.Column(this=exp.Identifier(this="b"))

        node = Contains(this=left, expression=right)

        assert isinstance(node, SpatialPredicate)
        assert isinstance(node, exp.Binary)

    def test_within_is_spatial_predicate_and_binary(self):
        """SP-003: Within inheritance.

        Given:
            Two expression nodes
        When:
            Within is instantiated
        Then:
            Instance is a SpatialPredicate and exp.Binary
        """
        left = exp.Column(this=exp.Identifier(this="a"))
        right = exp.Column(this=exp.Identifier(this="b"))

        node = Within(this=left, expression=right)

        assert isinstance(node, SpatialPredicate)
        assert isinstance(node, exp.Binary)


class TestSpatialSetPredicate:
    """Tests for SpatialSetPredicate expression node."""

    def test_instantiate_with_all_required_args(self):
        """SSP-001: Instantiate with all required args.

        Given:
            All required args (this, operator, quantifier, ranges)
        When:
            SpatialSetPredicate is instantiated
        Then:
            Instance has all four args accessible
        """
        this = exp.Column(this=exp.Identifier(this="interval"))
        operator = exp.Literal.string("INTERSECTS")
        quantifier = exp.Literal.string("ANY")
        ranges = exp.Array(
            expressions=[
                exp.Literal.string("chr1:1000-2000"),
                exp.Literal.string("chr1:5000-6000"),
            ]
        )

        node = SpatialSetPredicate(
            this=this,
            operator=operator,
            quantifier=quantifier,
            ranges=ranges,
        )

        assert node.args["this"] is this
        assert node.args["operator"] is operator
        assert node.args["quantifier"] is quantifier
        assert node.args["ranges"] is ranges


class TestGIQLCluster:
    """Tests for GIQLCluster expression node parsing."""

    def test_parse_cluster_with_one_arg(self):
        """CL-001: Parse CLUSTER with one positional arg.

        Given:
            A CLUSTER expression with one positional arg (column)
        When:
            Parsed with GIQLDialect
        Then:
            GIQLCluster instance has `this` set
        """
        ast = parse_one(
            "SELECT CLUSTER(interval) FROM features",
            dialect=GIQLDialect,
        )

        nodes = list(ast.find_all(GIQLCluster))
        assert len(nodes) == 1
        assert nodes[0].args["this"] is not None

    def test_parse_cluster_with_distance(self):
        """CL-002: Parse CLUSTER with distance.

        Given:
            A CLUSTER expression with two positional args (column, distance)
        When:
            Parsed with GIQLDialect
        Then:
            GIQLCluster instance has `this` and `distance` set
        """
        ast = parse_one(
            "SELECT CLUSTER(interval, 1000) FROM features",
            dialect=GIQLDialect,
        )

        nodes = list(ast.find_all(GIQLCluster))
        assert len(nodes) == 1
        assert nodes[0].args["this"] is not None
        assert nodes[0].args["distance"].this == "1000"

    def test_parse_cluster_with_stranded(self):
        """CL-003: Parse CLUSTER with stranded parameter.

        Given:
            A CLUSTER expression with one positional and stranded := true
        When:
            Parsed with GIQLDialect
        Then:
            GIQLCluster instance has `this` and `stranded` set
        """
        ast = parse_one(
            "SELECT CLUSTER(interval, stranded := true) FROM features",
            dialect=GIQLDialect,
        )

        nodes = list(ast.find_all(GIQLCluster))
        assert len(nodes) == 1
        assert nodes[0].args["this"] is not None
        assert nodes[0].args["stranded"] is not None

    def test_parse_cluster_with_distance_and_stranded(self):
        """CL-004: Parse CLUSTER with distance and stranded.

        Given:
            A CLUSTER expression with two positionals and stranded := true
        When:
            Parsed with GIQLDialect
        Then:
            GIQLCluster instance has `this`, `distance`, and `stranded` set
        """
        ast = parse_one(
            "SELECT CLUSTER(interval, 1000, stranded := true) FROM features",
            dialect=GIQLDialect,
        )

        nodes = list(ast.find_all(GIQLCluster))
        assert len(nodes) == 1
        assert nodes[0].args["this"] is not None
        assert nodes[0].args["distance"].this == "1000"
        assert nodes[0].args["stranded"] is not None

    def test_direct_instantiation_minimal(self):
        """CL-005: Direct instantiation with just `this`.

        Given:
            Required arg `this` only
        When:
            GIQLCluster is instantiated directly
        Then:
            Instance has `this` set; `distance` and `stranded` are absent
        """
        col = exp.Column(this=exp.Identifier(this="interval"))

        node = GIQLCluster(this=col)

        assert node.args["this"] is col
        assert node.args.get("distance") is None
        assert node.args.get("stranded") is None


class TestGIQLMerge:
    """Tests for GIQLMerge expression node parsing."""

    def test_parse_merge_with_one_arg(self):
        """MG-001: Parse MERGE with one positional arg.

        Given:
            A MERGE expression with one positional arg (column)
        When:
            Parsed with GIQLDialect
        Then:
            GIQLMerge instance has `this` set
        """
        ast = parse_one(
            "SELECT MERGE(interval) FROM features",
            dialect=GIQLDialect,
        )

        nodes = list(ast.find_all(GIQLMerge))
        assert len(nodes) == 1
        assert nodes[0].args["this"] is not None

    def test_parse_merge_with_distance(self):
        """MG-002: Parse MERGE with distance.

        Given:
            A MERGE expression with two positional args (column, distance)
        When:
            Parsed with GIQLDialect
        Then:
            GIQLMerge instance has `this` and `distance` set
        """
        ast = parse_one(
            "SELECT MERGE(interval, 1000) FROM features",
            dialect=GIQLDialect,
        )

        nodes = list(ast.find_all(GIQLMerge))
        assert len(nodes) == 1
        assert nodes[0].args["this"] is not None
        assert nodes[0].args["distance"].this == "1000"

    def test_parse_merge_with_stranded(self):
        """MG-003: Parse MERGE with stranded parameter.

        Given:
            A MERGE expression with one positional and stranded := true
        When:
            Parsed with GIQLDialect
        Then:
            GIQLMerge instance has `this` and `stranded` set
        """
        ast = parse_one(
            "SELECT MERGE(interval, stranded := true) FROM features",
            dialect=GIQLDialect,
        )

        nodes = list(ast.find_all(GIQLMerge))
        assert len(nodes) == 1
        assert nodes[0].args["this"] is not None
        assert nodes[0].args["stranded"] is not None

    def test_parse_merge_with_distance_and_stranded(self):
        """MG-004: Parse MERGE with distance and stranded.

        Given:
            A MERGE expression with two positionals and stranded := true
        When:
            Parsed with GIQLDialect
        Then:
            GIQLMerge instance has `this`, `distance`, and `stranded` set
        """
        ast = parse_one(
            "SELECT MERGE(interval, 1000, stranded := true) FROM features",
            dialect=GIQLDialect,
        )

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

    def test_from_arg_list_with_positional_args(self):
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

    def test_from_arg_list_with_walrus_named_stat(self):
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

    def test_from_arg_list_with_arrow_named_stat(self):
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

    def test_from_arg_list_with_named_resolution(self):
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

    def test_from_arg_list_with_walrus_named_target(self):
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

    def test_from_arg_list_with_arrow_named_target(self):
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

    def test_from_arg_list_with_all_named_params(self):
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
    @settings(suppress_health_check=[HealthCheck.function_scoped_fixture])
    def test_from_arg_list_with_varying_stat_and_resolution(
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
    @settings(suppress_health_check=[HealthCheck.function_scoped_fixture])
    def test_from_arg_list_with_varying_positional_only(self, resolution):
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
    @settings(suppress_health_check=[HealthCheck.function_scoped_fixture])
    def test_from_arg_list_with_varying_target_syntax(self, syntax):
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

    def test_parse_distance_with_two_positional_args(self):
        """DI-001: Parse DISTANCE with two positional args.

        Given:
            A DISTANCE expression with two positional args (interval_a, interval_b)
        When:
            Parsed with GIQLDialect
        Then:
            GIQLDistance instance has `this` and `expression` set
        """
        ast = parse_one(
            "SELECT DISTANCE(a.interval, b.interval) FROM a, b",
            dialect=GIQLDialect,
        )

        nodes = list(ast.find_all(GIQLDistance))
        assert len(nodes) == 1
        assert nodes[0].args["this"] is not None
        assert nodes[0].args["expression"] is not None

    def test_parse_distance_with_stranded_and_signed(self):
        """DI-002: Parse DISTANCE with stranded and signed.

        Given:
            A DISTANCE expression with two positionals and stranded := true, signed := true
        When:
            Parsed with GIQLDialect
        Then:
            GIQLDistance instance has `this`, `expression`, `stranded`, and `signed` set
        """
        ast = parse_one(
            "SELECT DISTANCE(a.interval, b.interval, stranded := true, signed := true) FROM a, b",
            dialect=GIQLDialect,
        )

        nodes = list(ast.find_all(GIQLDistance))
        assert len(nodes) == 1
        assert nodes[0].args["this"] is not None
        assert nodes[0].args["expression"] is not None
        assert nodes[0].args["stranded"] is not None
        assert nodes[0].args["signed"] is not None

    def test_parse_distance_with_stranded_only(self):
        """DI-003: Parse DISTANCE with only stranded.

        Given:
            A DISTANCE expression with two positionals and only stranded := true
        When:
            Parsed with GIQLDialect
        Then:
            GIQLDistance instance has `this`, `expression`, and `stranded` set; `signed` absent
        """
        ast = parse_one(
            "SELECT DISTANCE(a.interval, b.interval, stranded := true) FROM a, b",
            dialect=GIQLDialect,
        )

        nodes = list(ast.find_all(GIQLDistance))
        assert len(nodes) == 1
        assert nodes[0].args["this"] is not None
        assert nodes[0].args["expression"] is not None
        assert nodes[0].args["stranded"] is not None
        assert nodes[0].args.get("signed") is None


class TestGIQLNearest:
    """Tests for GIQLNearest expression node parsing."""

    def test_parse_nearest_with_one_positional(self):
        """NR-001: Parse NEAREST with one positional arg.

        Given:
            A NEAREST expression with one positional arg (table)
        When:
            Parsed with GIQLDialect
        Then:
            GIQLNearest instance has `this` set
        """
        ast = parse_one(
            "SELECT NEAREST(genes) FROM peaks",
            dialect=GIQLDialect,
        )

        nodes = list(ast.find_all(GIQLNearest))
        assert len(nodes) == 1
        assert nodes[0].args["this"] is not None

    def test_parse_nearest_with_k(self):
        """NR-002: Parse NEAREST with k parameter.

        Given:
            A NEAREST expression with one positional and k := 3
        When:
            Parsed with GIQLDialect
        Then:
            GIQLNearest instance has `this` and `k` set
        """
        ast = parse_one(
            "SELECT NEAREST(genes, k := 3) FROM peaks",
            dialect=GIQLDialect,
        )

        nodes = list(ast.find_all(GIQLNearest))
        assert len(nodes) == 1
        assert nodes[0].args["this"] is not None
        assert nodes[0].args["k"].this == "3"

    def test_parse_nearest_with_multiple_named_params(self):
        """NR-003: Parse NEAREST with multiple named params.

        Given:
            A NEAREST expression with one positional and multiple named params
        When:
            Parsed with GIQLDialect
        Then:
            GIQLNearest instance has all provided args set
        """
        ast = parse_one(
            "SELECT NEAREST(genes, k := 5, max_distance := 100000, stranded := true, signed := true) FROM peaks",
            dialect=GIQLDialect,
        )

        nodes = list(ast.find_all(GIQLNearest))
        assert len(nodes) == 1
        assert nodes[0].args["this"] is not None
        assert nodes[0].args["k"].this == "5"
        assert nodes[0].args["max_distance"].this == "100000"
        assert nodes[0].args["stranded"] is not None
        assert nodes[0].args["signed"] is not None
