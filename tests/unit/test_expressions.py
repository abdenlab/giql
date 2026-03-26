"""Tests for custom AST expression nodes.

Test specification: specs/test_expressions.md
"""

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

    def test_parse_coverage_with_positional_args(self):
        """COV-001: Parse COVERAGE with positional args.

        Given:
            A COVERAGE expression with two positional args (column, resolution)
        When:
            Parsed with GIQLDialect
        Then:
            GIQLCoverage instance has `this` and `resolution` set
        """
        ast = parse_one(
            "SELECT COVERAGE(interval, 1000) FROM features",
            dialect=GIQLDialect,
        )

        nodes = list(ast.find_all(GIQLCoverage))
        assert len(nodes) == 1
        assert nodes[0].args["this"] is not None
        assert nodes[0].args["resolution"].this == "1000"
        assert nodes[0].args.get("stat") is None
        assert nodes[0].args.get("target") is None

    def test_parse_coverage_with_walrus_named_resolution(self):
        """COV-002: Parse COVERAGE with := named resolution.

        Given:
            A COVERAGE expression with one positional and resolution := 1000
        When:
            Parsed with GIQLDialect
        Then:
            GIQLCoverage instance has `this` and `resolution` set
        """
        ast = parse_one(
            "SELECT COVERAGE(interval, resolution := 1000) FROM features",
            dialect=GIQLDialect,
        )

        nodes = list(ast.find_all(GIQLCoverage))
        assert len(nodes) == 1
        assert nodes[0].args["this"] is not None
        assert nodes[0].args["resolution"].this == "1000"

    def test_parse_coverage_with_stat(self):
        """COV-003: Parse COVERAGE with stat parameter.

        Given:
            A COVERAGE expression with two positionals and stat := 'mean'
        When:
            Parsed with GIQLDialect
        Then:
            GIQLCoverage instance has `this`, `resolution`, and `stat` set
        """
        ast = parse_one(
            "SELECT COVERAGE(interval, 500, stat := 'mean') FROM features",
            dialect=GIQLDialect,
        )

        nodes = list(ast.find_all(GIQLCoverage))
        assert len(nodes) == 1
        assert nodes[0].args["resolution"].this == "500"
        assert nodes[0].args["stat"].this == "mean"

    def test_parse_coverage_with_stat_and_target(self):
        """COV-004: Parse COVERAGE with stat and target.

        Given:
            A COVERAGE expression with two positionals, stat := 'mean', and target := 'score'
        When:
            Parsed with GIQLDialect
        Then:
            GIQLCoverage instance has `this`, `resolution`, `stat`, and `target` set
        """
        ast = parse_one(
            "SELECT COVERAGE(interval, 1000, stat := 'mean', target := 'score') FROM features",
            dialect=GIQLDialect,
        )

        nodes = list(ast.find_all(GIQLCoverage))
        assert len(nodes) == 1
        assert nodes[0].args["resolution"].this == "1000"
        assert nodes[0].args["stat"].this == "mean"
        assert nodes[0].args["target"].this == "score"

    def test_parse_coverage_with_arrow_named_resolution(self):
        """COV-005: Parse COVERAGE with => named resolution.

        Given:
            A COVERAGE expression with one positional and resolution => 1000
        When:
            Parsed with GIQLDialect
        Then:
            GIQLCoverage instance has `this` and `resolution` set
        """
        ast = parse_one(
            "SELECT COVERAGE(interval, resolution => 1000) FROM features",
            dialect=GIQLDialect,
        )

        nodes = list(ast.find_all(GIQLCoverage))
        assert len(nodes) == 1
        assert nodes[0].args["this"] is not None
        assert nodes[0].args["resolution"].this == "1000"

    def test_parse_coverage_with_target_no_stat(self):
        """COV-006: Parse COVERAGE with target but no stat.

        Given:
            A COVERAGE expression with two positionals and target := 'score' only
        When:
            Parsed with GIQLDialect
        Then:
            GIQLCoverage instance has `this`, `resolution`, and `target` set; `stat` is absent
        """
        ast = parse_one(
            "SELECT COVERAGE(interval, 1000, target := 'score') FROM features",
            dialect=GIQLDialect,
        )

        nodes = list(ast.find_all(GIQLCoverage))
        assert len(nodes) == 1
        assert nodes[0].args["resolution"].this == "1000"
        assert nodes[0].args["target"].this == "score"
        assert nodes[0].args.get("stat") is None

    def test_direct_instantiation_minimal(self):
        """COV-007: Direct instantiation with required args only.

        Given:
            Required args `this` and `resolution` only
        When:
            GIQLCoverage is instantiated directly
        Then:
            Instance has `this` and `resolution` set; `stat` and `target` are absent
        """
        col = exp.Column(this=exp.Identifier(this="interval"))
        resolution = exp.Literal.number(1000)

        node = GIQLCoverage(this=col, resolution=resolution)

        assert node.args["this"] is col
        assert node.args["resolution"] is resolution
        assert node.args.get("stat") is None
        assert node.args.get("target") is None


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
