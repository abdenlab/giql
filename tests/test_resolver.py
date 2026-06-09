"""Tests for the ResolveOperatorRefs normalization pass.

Tests verify that pass 1 attaches well-formed resolution metadata to every GIQL
operator slot, that its resolution of the table-shaped reference slots mirrors
the generator's existing behavior, and that the closing validation boundary
accepts well-formed metadata while rejecting malformed metadata.
"""

import pytest
from sqlglot import exp
from sqlglot import parse_one

from giql.dialect import GIQLDialect
from giql.expressions import GIQLDisjoin
from giql.expressions import GIQLDistance
from giql.expressions import GIQLNearest
from giql.expressions import SpatialSetPredicate
from giql.resolver import META_KEY
from giql.resolver import OperatorResolution
from giql.resolver import ResolutionError
from giql.resolver import ResolvedRef
from giql.resolver import resolve_operator_refs
from giql.resolver import validate_operator_refs
from giql.table import Table
from giql.table import Tables

hypothesis = pytest.importorskip("hypothesis")
from hypothesis import given  # noqa: E402
from hypothesis import settings  # noqa: E402
from hypothesis import strategies as st  # noqa: E402

# The g_ prefix keeps generated identifiers clear of SQL keywords, which would
# otherwise break parsing when interpolated as a table name.
_identifiers = st.from_regex(r"g_[a-z0-9_]{0,8}", fullmatch=True)


def _tables(*names: str) -> Tables:
    """Build a Tables container registering each name with default columns."""
    container = Tables()
    for name in names:
        container.register(name, Table(name))
    return container


def _disjoin_node(ast: exp.Expression) -> GIQLDisjoin:
    """Return the single GIQLDisjoin node reachable from an annotated AST."""
    return next(n for n in ast.walk() if isinstance(n, GIQLDisjoin))


class TestResolveOperatorRefs:
    """Tests for the resolve_operator_refs pass."""

    def test_resolve_operator_refs_attaches_metadata_to_every_operator(self):
        """Test that every GIQL operator receives resolution metadata.

        Given:
            A query whose WHERE clause carries a spatial predicate.
        When:
            Running the resolve pass over the parsed AST.
        Then:
            It should attach an OperatorResolution under the giql meta key.
        """
        # Arrange
        ast = parse_one(
            "SELECT * FROM peaks WHERE interval INTERSECTS 'chr1:1000-2000'",
            dialect=GIQLDialect,
        )

        # Act
        resolve_operator_refs(ast, _tables("peaks"))

        # Assert
        predicate = next(n for n in ast.walk() if type(n).__name__ == "Intersects")
        assert isinstance(predicate.meta.get(META_KEY), OperatorResolution)

    def test_resolve_operator_refs_resolves_disjoin_target(self):
        """Test that a DISJOIN target resolves to its registered table.

        Given:
            A DISJOIN over a registered table with custom genomic columns.
        When:
            Running the resolve pass.
        Then:
            It should attach a registered-table ResolvedRef carrying the
            table's physical column names.
        """
        # Arrange
        tables = Tables()
        tables.register(
            "features",
            Table("features", chrom_col="chr", start_col="lo", end_col="hi"),
        )
        ast = parse_one("SELECT * FROM DISJOIN(features)", dialect=GIQLDialect)

        # Act
        resolve_operator_refs(ast, tables)

        # Assert
        ref = _disjoin_node(ast).meta[META_KEY].slot("this")
        assert ref == ResolvedRef(
            kind="registered_table",
            name="features",
            cols=("chr", "lo", "hi"),
            table=tables.get("features"),
            coverage_skippable=False,
        )

    def test_resolve_operator_refs_omitted_reference_is_self_reference(self):
        """Test that an omitted DISJOIN reference is a coverage-skippable self-reference.

        Given:
            A DISJOIN with no reference argument.
        When:
            Running the resolve pass.
        Then:
            It should resolve the reference slot to the target relation with
            coverage_skippable set.
        """
        # Arrange
        ast = parse_one("SELECT * FROM DISJOIN(features)", dialect=GIQLDialect)

        # Act
        resolve_operator_refs(ast, _tables("features"))

        # Assert
        ref = _disjoin_node(ast).meta[META_KEY].slot("reference")
        assert ref.kind == "registered_table"
        assert ref.name == "features"
        assert ref.coverage_skippable is True

    def test_resolve_operator_refs_distinct_registered_reference(self):
        """Test that a distinct registered DISJOIN reference is not coverage-skippable.

        Given:
            A DISJOIN whose reference is a different registered table.
        When:
            Running the resolve pass.
        Then:
            It should resolve to a registered-table reference that is not a
            self-reference.
        """
        # Arrange
        ast = parse_one(
            "SELECT * FROM DISJOIN(features, reference := other)",
            dialect=GIQLDialect,
        )

        # Act
        resolve_operator_refs(ast, _tables("features", "other"))

        # Assert
        ref = _disjoin_node(ast).meta[META_KEY].slot("reference")
        assert ref.kind == "registered_table"
        assert ref.name == "other"
        assert ref.coverage_skippable is False

    def test_resolve_operator_refs_cte_reference_shadows_table(self):
        """Test that an enclosing CTE reference resolves to a CTE.

        Given:
            A DISJOIN reference naming a CTE that shadows a registered table of
            the same name.
        When:
            Running the resolve pass.
        Then:
            It should resolve to a canonical CTE reference carrying no Table
            config.
        """
        # Arrange
        ast = parse_one(
            "WITH other AS (SELECT 1) SELECT * FROM DISJOIN(features, reference := other)",
            dialect=GIQLDialect,
        )

        # Act
        resolve_operator_refs(ast, _tables("features", "other"))

        # Assert
        ref = _disjoin_node(ast).meta[META_KEY].slot("reference")
        assert ref.kind == "cte"
        assert ref.name == "other"
        assert ref.table is None
        assert ref.cols == ("chrom", "start", "end")

    def test_resolve_operator_refs_subquery_reference(self):
        """Test that a subquery DISJOIN reference resolves to a subquery ref.

        Given:
            A DISJOIN whose reference is a parenthesized SELECT.
        When:
            Running the resolve pass.
        Then:
            It should resolve to an anonymous canonical subquery reference.
        """
        # Arrange
        ast = parse_one(
            "SELECT * FROM DISJOIN(features, reference := (SELECT * FROM other))",
            dialect=GIQLDialect,
        )

        # Act
        resolve_operator_refs(ast, _tables("features"))

        # Assert
        ref = _disjoin_node(ast).meta[META_KEY].slot("reference")
        assert ref.kind == "subquery"
        assert ref.name is None
        assert ref.table is None

    def test_resolve_operator_refs_unregistered_target_left_unresolved(self):
        """Test that an unregistered target leaves the slot unresolved.

        Given:
            A DISJOIN whose target table is not registered.
        When:
            Running the resolve pass.
        Then:
            It should attach an OperatorResolution with no resolved slots,
            deferring the error to the generator.
        """
        # Arrange
        ast = parse_one("SELECT * FROM DISJOIN(missing)", dialect=GIQLDialect)

        # Act
        resolve_operator_refs(ast, Tables())

        # Assert
        resolution = _disjoin_node(ast).meta[META_KEY]
        assert resolution.slots == {}

    def test_resolve_operator_refs_reserved_prefix_reference_left_unresolved(self):
        """Test that a reserved-prefix reference name leaves the slot unresolved.

        Given:
            A DISJOIN reference whose name uses the reserved __giql_dj_ prefix.
        When:
            Running the resolve pass.
        Then:
            It should leave the reference slot unresolved, deferring the error
            to the generator.
        """
        # Arrange
        ast = parse_one(
            "SELECT * FROM DISJOIN(features, reference := __giql_dj_ref)",
            dialect=GIQLDialect,
        )

        # Act
        resolve_operator_refs(ast, _tables("features"))

        # Assert
        resolution = _disjoin_node(ast).meta[META_KEY]
        assert resolution.slot("this") is not None
        assert resolution.slot("reference") is None

    def test_resolve_operator_refs_resolves_nearest_target(self):
        """Test that a NEAREST target resolves to its registered table.

        Given:
            A standalone NEAREST over a registered table.
        When:
            Running the resolve pass.
        Then:
            It should attach a registered-table ResolvedRef for the target slot.
        """
        # Arrange
        ast = parse_one(
            "SELECT * FROM NEAREST(genes, reference := 'chr1:1-2', k := 3)",
            dialect=GIQLDialect,
        )

        # Act
        resolve_operator_refs(ast, _tables("genes"))

        # Assert
        nearest = next(n for n in ast.walk() if isinstance(n, GIQLNearest))
        ref = nearest.meta[META_KEY].slot("this")
        assert ref.kind == "registered_table"
        assert ref.name == "genes"

    def test_resolve_operator_refs_returns_same_expression(self):
        """Test that the pass annotates and returns the same AST in place.

        Given:
            A parsed GIQL query.
        When:
            Running the resolve pass.
        Then:
            It should return the identical expression object it was given.
        """
        # Arrange
        ast = parse_one("SELECT * FROM DISJOIN(features)", dialect=GIQLDialect)

        # Act
        result = resolve_operator_refs(ast, _tables("features"))

        # Assert
        assert result is ast

    def test_resolve_operator_refs_outer_cte_visible_in_derived_table(self):
        """Test that an enclosing WITH's CTE is visible to a nested DISJOIN.

        Given:
            A DISJOIN nested in a derived table whose reference names a CTE
            defined by the enclosing outer WITH clause.
        When:
            Running the resolve pass.
        Then:
            It should resolve the reference to a CTE.
        """
        # Arrange
        ast = parse_one(
            "WITH mask AS (SELECT 1) "
            "SELECT * FROM (SELECT * FROM DISJOIN(features, reference := mask)) AS d",
            dialect=GIQLDialect,
        )

        # Act
        resolve_operator_refs(ast, _tables("features"))

        # Assert
        ref = _disjoin_node(ast).meta[META_KEY].slot("reference")
        assert ref.kind == "cte"
        assert ref.name == "mask"

    def test_resolve_operator_refs_metadata_survives_copy(self):
        """Test that attached resolution metadata survives a tree copy.

        Given:
            An AST annotated by the resolve pass.
        When:
            Copying the AST with Expression.copy().
        Then:
            It should carry an equal OperatorResolution on the copied operator.
        """
        # Arrange
        ast = parse_one("SELECT * FROM DISJOIN(features)", dialect=GIQLDialect)
        resolve_operator_refs(ast, _tables("features"))

        # Act
        copied = ast.copy()

        # Assert
        assert _disjoin_node(copied).meta[META_KEY] == _disjoin_node(ast).meta[META_KEY]

    def test_resolve_operator_refs_annotates_distance_and_set_predicate(self):
        """Test that non-reference-slot operators still receive metadata.

        Given:
            A query using DISTANCE in the projection and a spatial set
            predicate in the WHERE clause.
        When:
            Running the resolve pass.
        Then:
            It should attach an empty-slot OperatorResolution to both
            operators.
        """
        # Arrange
        ast = parse_one(
            "SELECT DISTANCE(a.interval, b.interval) FROM a, b "
            "WHERE a.interval INTERSECTS ANY('chr1:1-2', 'chr2:3-4')",
            dialect=GIQLDialect,
        )

        # Act
        resolve_operator_refs(ast, _tables("a", "b"))

        # Assert
        distance = next(n for n in ast.walk() if isinstance(n, GIQLDistance))
        predicate = next(n for n in ast.walk() if isinstance(n, SpatialSetPredicate))
        assert distance.meta[META_KEY].slots == {}
        assert predicate.meta[META_KEY].slots == {}

    def test_resolve_operator_refs_unregistered_nearest_target_left_unresolved(self):
        """Test that an unregistered NEAREST target leaves the slot unresolved.

        Given:
            A NEAREST whose target table is not registered.
        When:
            Running the resolve pass.
        Then:
            It should attach an OperatorResolution with no resolved slots,
            deferring the error to the generator.
        """
        # Arrange
        ast = parse_one(
            "SELECT * FROM NEAREST(missing, reference := 'chr1:1-2')",
            dialect=GIQLDialect,
        )

        # Act
        resolve_operator_refs(ast, Tables())

        # Assert
        nearest = next(n for n in ast.walk() if isinstance(n, GIQLNearest))
        assert nearest.meta[META_KEY].slots == {}

    @settings(max_examples=50)
    @given(names=st.lists(_identifiers, min_size=4, max_size=4, unique=True))
    def test_resolve_operator_refs_with_arbitrary_table_config(self, names):
        """Test resolution against arbitrary table and column names.

        Given:
            Any registered Table with generated identifier-safe table and
            chrom/start/end column names.
        When:
            Resolving a DISJOIN over that table.
        Then:
            The target ref should carry exactly the configured column names
            and validation should pass.
        """
        # Arrange
        table_name, chrom, start, end = names
        tables = Tables()
        tables.register(
            table_name,
            Table(table_name, chrom_col=chrom, start_col=start, end_col=end),
        )
        ast = parse_one(f"SELECT * FROM DISJOIN({table_name})", dialect=GIQLDialect)

        # Act
        resolve_operator_refs(ast, tables)

        # Assert
        ref = _disjoin_node(ast).meta[META_KEY].slot("this")
        assert ref.cols == (chrom, start, end)
        validate_operator_refs(ast)


class TestValidateOperatorRefs:
    """Tests for the validate_operator_refs validation boundary."""

    def test_validate_operator_refs_accepts_resolved_tree(self):
        """Test that validation accepts a well-formed annotated tree.

        Given:
            An AST already annotated by the resolve pass.
        When:
            Running validation over it.
        Then:
            It should not raise.
        """
        # Arrange
        ast = parse_one(
            "SELECT * FROM DISJOIN(features, reference := other)",
            dialect=GIQLDialect,
        )
        resolve_operator_refs(ast, _tables("features", "other"))

        # Act & assert
        validate_operator_refs(ast)

    def test_validate_operator_refs_rejects_missing_metadata(self):
        """Test that validation rejects an operator with no metadata.

        Given:
            A DISJOIN node that was never annotated by the resolve pass.
        When:
            Running validation over the tree.
        Then:
            It should raise a ResolutionError naming the missing metadata.
        """
        # Arrange
        ast = parse_one("SELECT * FROM DISJOIN(features)", dialect=GIQLDialect)

        # Act & assert
        with pytest.raises(ResolutionError, match="missing resolution metadata"):
            validate_operator_refs(ast)

    def test_validate_operator_refs_rejects_disallowed_kind(self):
        """Test that validation rejects a slot resolved to a disallowed kind.

        Given:
            A DISJOIN target slot annotated with a subquery ResolvedRef, a kind
            the target slot does not accept.
        When:
            Running validation over the tree.
        Then:
            It should raise a ResolutionError naming the rejected kind.
        """
        # Arrange
        ast = parse_one("SELECT * FROM DISJOIN(features)", dialect=GIQLDialect)
        node = _disjoin_node(ast)
        node.meta[META_KEY] = OperatorResolution(
            "GIQLDisjoin",
            {
                "this": ResolvedRef(
                    kind="subquery",
                    name=None,
                    cols=("chrom", "start", "end"),
                    table=None,
                    coverage_skippable=False,
                )
            },
        )

        # Act & assert
        with pytest.raises(ResolutionError, match="not accepted"):
            validate_operator_refs(ast)

    def test_validate_operator_refs_rejects_registered_table_without_config(self):
        """Test that validation rejects a registered-table ref with no Table.

        Given:
            A DISJOIN target annotated as a registered table but carrying no
            Table config.
        When:
            Running validation over the tree.
        Then:
            It should raise a ResolutionError naming the missing config.
        """
        # Arrange
        ast = parse_one("SELECT * FROM DISJOIN(features)", dialect=GIQLDialect)
        node = _disjoin_node(ast)
        node.meta[META_KEY] = OperatorResolution(
            "GIQLDisjoin",
            {
                "this": ResolvedRef(
                    kind="registered_table",
                    name="features",
                    cols=("chrom", "start", "end"),
                    table=None,
                    coverage_skippable=False,
                )
            },
        )

        # Act & assert
        with pytest.raises(ResolutionError, match="no Table config"):
            validate_operator_refs(ast)

    def test_validate_operator_refs_rejects_cte_with_table_config(self):
        """Test that validation rejects a CTE ref carrying a Table config.

        Given:
            A DISJOIN reference annotated as a CTE but carrying a Table
            config, contradicting the canonical assumption for CTEs.
        When:
            Running validation over the tree.
        Then:
            It should raise a ResolutionError naming the canonical assumption.
        """
        # Arrange
        ast = parse_one("SELECT * FROM DISJOIN(features)", dialect=GIQLDialect)
        node = _disjoin_node(ast)
        table = Table("features")
        node.meta[META_KEY] = OperatorResolution(
            "GIQLDisjoin",
            {
                "reference": ResolvedRef(
                    kind="cte",
                    name="features",
                    cols=("chrom", "start", "end"),
                    table=table,
                    coverage_skippable=False,
                )
            },
        )

        # Act & assert
        with pytest.raises(ResolutionError, match="assumed canonical"):
            validate_operator_refs(ast)

    def test_validate_operator_refs_rejects_malformed_cols(self):
        """Test that validation rejects a ref with malformed cols.

        Given:
            A DISJOIN target annotated with a ResolvedRef whose cols is not a
            3-tuple of column names.
        When:
            Running validation over the tree.
        Then:
            It should raise a ResolutionError naming the malformed cols.
        """
        # Arrange
        ast = parse_one("SELECT * FROM DISJOIN(features)", dialect=GIQLDialect)
        node = _disjoin_node(ast)
        table = Table("features")
        node.meta[META_KEY] = OperatorResolution(
            "GIQLDisjoin",
            {
                "this": ResolvedRef(
                    kind="registered_table",
                    name="features",
                    cols=("chrom", "start"),
                    table=table,
                    coverage_skippable=False,
                )
            },
        )

        # Act & assert
        with pytest.raises(ResolutionError, match="malformed cols"):
            validate_operator_refs(ast)

    def test_validate_operator_refs_rejects_non_resolvedref_slot_value(self):
        """Test that validation rejects a slot holding a non-ResolvedRef value.

        Given:
            A DISJOIN target slot whose resolution metadata holds a plain
            string instead of a ResolvedRef.
        When:
            Running validation over the tree.
        Then:
            It should raise a ResolutionError naming the expected type.
        """
        # Arrange
        ast = parse_one("SELECT * FROM DISJOIN(features)", dialect=GIQLDialect)
        node = _disjoin_node(ast)
        node.meta[META_KEY] = OperatorResolution(
            "GIQLDisjoin",
            {"this": "features"},
        )

        # Act & assert
        with pytest.raises(ResolutionError, match="expected ResolvedRef"):
            validate_operator_refs(ast)
