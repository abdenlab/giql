"""Tests for the operator-expander registry and the ``ExpandOperators`` pass.

These tests pin the dispatch infrastructure of epic #137, step 2:

* the :class:`~giql.expander.ExpanderRegistry` resolves an expander through the
  ``(target, op)`` → ``(generic, op)`` → ``None`` (legacy) fallback chain, and a
  later registration overrides an earlier one;
* the :func:`~giql.expander.register` decorator writes to the process-wide
  registry and returns the expander unchanged;
* the :class:`~giql.expander.ExpandOperators` pass dispatches an opted-in
  operator to its registered expander — including a *fake* custom target plus a
  fake operator, the public extension hook;
* with nothing flagged and an empty registry the pass is a strict no-op and the
  emitted SQL is byte-identical.
"""

from dataclasses import dataclass

import pytest
from sqlglot import exp
from sqlglot import parse_one

from giql.dialect import GIQLDialect
from giql.expander import EXPAND_ALIAS_PREFIX
from giql.expander import REGISTRY
from giql.expander import ExpanderRegistry
from giql.expander import ExpandOperators
from giql.expander import ExpansionContext
from giql.expander import OperatorExpander
from giql.expander import expand_operators
from giql.expander import register
from giql.expressions import GIQLDisjoin
from giql.resolver import resolve_operator_refs
from giql.table import Table
from giql.table import Tables
from giql.targets import Capabilities
from giql.targets import DataFusionTarget
from giql.targets import DuckDBTarget
from giql.targets import GenericTarget
from giql.targets import Target
from giql.transpile import transpile


def _record(label: str):
    """Build an expander function that replaces a node with a tagged literal."""

    def _expander(node: exp.Expression, ctx: ExpansionContext) -> exp.Expression:
        return exp.Literal.string(label)

    return _expander


@pytest.fixture
def clean_registry():
    """Isolate the process-wide REGISTRY, restoring its contents afterward."""
    saved = dict(REGISTRY._expanders)
    REGISTRY.clear()
    yield REGISTRY
    REGISTRY.clear()
    REGISTRY._expanders.update(saved)


class TestExpanderRegistry:
    """Tests for ExpanderRegistry registration and fallback resolution."""

    def test_resolve_returns_registered_target_expander(self):
        """Test that resolve returns the exact-target expander.

        Given:
            A registry with an expander registered for (DuckDBTarget, op).
        When:
            Resolving (DuckDBTarget(), op).
        Then:
            It should return that exact-target expander.
        """
        # Arrange
        registry = ExpanderRegistry()
        fn = _record("duckdb")
        registry.register(DuckDBTarget(), GIQLDisjoin, fn)

        # Act
        resolved = registry.resolve(DuckDBTarget(), GIQLDisjoin)

        # Assert
        assert resolved is fn

    def test_resolve_falls_back_to_generic_expander(self):
        """Test that resolve falls back from a target to the generic expander.

        Given:
            A registry with an expander registered only for (GenericTarget, op).
        When:
            Resolving (DuckDBTarget(), op), for which no exact entry exists.
        Then:
            It should fall back to the generic expander.
        """
        # Arrange
        registry = ExpanderRegistry()
        generic_fn = _record("generic")
        registry.register(GenericTarget(), GIQLDisjoin, generic_fn)

        # Act
        resolved = registry.resolve(DuckDBTarget(), GIQLDisjoin)

        # Assert
        assert resolved is generic_fn

    def test_resolve_prefers_target_over_generic(self):
        """Test that an exact-target entry shadows the generic fallback.

        Given:
            A registry with both a (DuckDBTarget, op) and a (GenericTarget, op)
            expander registered.
        When:
            Resolving (DuckDBTarget(), op).
        Then:
            It should return the exact-target expander, not the generic one.
        """
        # Arrange
        registry = ExpanderRegistry()
        duckdb_fn = _record("duckdb")
        generic_fn = _record("generic")
        registry.register(DuckDBTarget(), GIQLDisjoin, duckdb_fn)
        registry.register(GenericTarget(), GIQLDisjoin, generic_fn)

        # Act
        resolved = registry.resolve(DuckDBTarget(), GIQLDisjoin)

        # Assert
        assert resolved is duckdb_fn

    def test_resolve_returns_none_for_legacy_fallback(self):
        """Test that resolve returns None when no expander is registered.

        Given:
            An empty registry.
        When:
            Resolving any (target, op).
        Then:
            It should return None, signalling the legacy *_sql emitter fallback.
        """
        # Arrange
        registry = ExpanderRegistry()

        # Act
        resolved = registry.resolve(DuckDBTarget(), GIQLDisjoin)

        # Assert
        assert resolved is None

    def test_resolve_generic_target_does_not_self_fallback(self):
        """Test that a generic-target resolve does not double-check generic.

        Given:
            An empty registry.
        When:
            Resolving (GenericTarget(), op).
        Then:
            It should return None without raising (no redundant generic lookup).
        """
        # Arrange
        registry = ExpanderRegistry()

        # Act
        resolved = registry.resolve(GenericTarget(), GIQLDisjoin)

        # Assert
        assert resolved is None

    def test_register_overrides_existing_entry(self):
        """Test that a later registration overrides an earlier one.

        Given:
            A registry with an expander already registered for (target, op).
        When:
            Registering a second expander for the same key.
        Then:
            It should resolve to the second expander (last-write-wins).
        """
        # Arrange
        registry = ExpanderRegistry()
        first = _record("first")
        second = _record("second")
        registry.register(DuckDBTarget(), GIQLDisjoin, first)

        # Act
        registry.register(DuckDBTarget(), GIQLDisjoin, second)

        # Assert
        assert registry.resolve(DuckDBTarget(), GIQLDisjoin) is second

    def test_register_keys_on_target_value_not_identity(self):
        """Test that the registry keys on the frozen target value.

        Given:
            An expander registered against one DuckDBTarget() instance.
        When:
            Resolving against a distinct, freshly built DuckDBTarget() instance.
        Then:
            It should hit the same entry (frozen value-equal target keys).
        """
        # Arrange
        registry = ExpanderRegistry()
        fn = _record("duckdb")
        registry.register(DuckDBTarget(), GIQLDisjoin, fn)

        # Act
        resolved = registry.resolve(DuckDBTarget(), GIQLDisjoin)

        # Assert
        assert DuckDBTarget() is not DuckDBTarget()
        assert resolved is fn

    def test_register_accepts_expander_object(self):
        """Test that an OperatorExpander object registers and resolves.

        Given:
            An object with an expand method (an OperatorExpander).
        When:
            Registering it and resolving the key.
        Then:
            It should resolve to a callable that delegates to expand.
        """

        # Arrange
        class _Expander:
            def expand(self, node, ctx):
                return exp.Literal.string("object")

        registry = ExpanderRegistry()
        obj = _Expander()
        assert isinstance(obj, OperatorExpander)
        registry.register(GenericTarget(), GIQLDisjoin, obj)

        # Act
        resolved = registry.resolve(GenericTarget(), GIQLDisjoin)
        result = resolved(None, None)

        # Assert
        assert result == exp.Literal.string("object")

    def test_register_rejects_non_callable(self):
        """Test that registering a non-callable, non-expander raises.

        Given:
            A value that is neither callable nor an OperatorExpander.
        When:
            Registering it.
        Then:
            It should raise TypeError.
        """
        # Arrange
        registry = ExpanderRegistry()

        # Act & assert
        with pytest.raises(TypeError):
            registry.register(GenericTarget(), GIQLDisjoin, object())


class TestRegisterDecorator:
    """Tests for the @register extension-hook decorator."""

    def test_register_writes_to_process_registry(self, clean_registry):
        """Test that the decorator registers in the process-wide registry.

        Given:
            The decorator applied to an expander for (DuckDBTarget, op).
        When:
            Resolving that key against the process-wide registry.
        Then:
            It should resolve to the decorated expander.
        """

        # Arrange & act
        @register(DuckDBTarget, GIQLDisjoin)
        def _expander(node, ctx):
            return exp.Literal.string("decorated")

        # Assert
        assert clean_registry.resolve(DuckDBTarget(), GIQLDisjoin) is _expander

    def test_register_returns_expander_unchanged(self, clean_registry):
        """Test that the decorator returns its target unchanged.

        Given:
            An expander function.
        When:
            Decorating it with @register.
        Then:
            It should return the same function object (usable standalone).
        """

        # Arrange
        def _expander(node, ctx):
            return exp.Literal.string("x")

        # Act
        decorated = register(GenericTarget, GIQLDisjoin)(_expander)

        # Assert
        assert decorated is _expander

    def test_register_accepts_target_instance(self, clean_registry):
        """Test that the decorator accepts a target instance, not just a class.

        Given:
            The decorator applied with an explicit DuckDBTarget() instance.
        When:
            Resolving that key.
        Then:
            It should resolve to the decorated expander.
        """

        # Arrange & act
        @register(DuckDBTarget(), GIQLDisjoin)
        def _expander(node, ctx):
            return exp.Literal.string("instance")

        # Assert
        assert clean_registry.resolve(DuckDBTarget(), GIQLDisjoin) is _expander


class TestExpansionContext:
    """Tests for the ExpansionContext value object."""

    def test_capabilities_returns_target_capabilities(self):
        """Test that capabilities exposes the active target's capability set.

        Given:
            A context built with a DuckDBTarget.
        When:
            Reading its capabilities property.
        Then:
            It should be the target's Capabilities.
        """
        # Arrange
        target = DuckDBTarget()
        ctx = ExpansionContext(exp.true(), None, target, Tables())

        # Act
        caps = ctx.capabilities

        # Assert
        assert isinstance(caps, Capabilities)
        assert caps is target.capabilities

    def test_alias_mints_unique_prefixed_names(self):
        """Test that alias mints distinct names under the reserved prefix.

        Given:
            A single ExpansionContext.
        When:
            Minting two aliases.
        Then:
            They should differ and both carry the expander prefix.
        """
        # Arrange
        ctx = ExpansionContext(exp.true(), None, GenericTarget(), Tables())

        # Act
        first = ctx.alias()
        second = ctx.alias()

        # Assert
        assert first != second
        assert first.startswith(EXPAND_ALIAS_PREFIX)
        assert second.startswith(EXPAND_ALIAS_PREFIX)


# A fake custom target standing in for a user-defined engine: the extension-hook
# seam must dispatch to it exactly as to the built-in targets.
@dataclass(frozen=True)
class _FakeTarget(Target):
    name: str = "fake"
    sqlglot_dialect: str | None = None
    capabilities: Capabilities = Capabilities(
        supports_lateral=True,
        supports_star_replace=False,
        supports_qualify=False,
        range_join_strategy="binned",
    )


class TestExpandOperatorsPass:
    """Tests for the ExpandOperators dispatch pass."""

    def test_transform_dispatches_to_registered_expander(self, clean_registry):
        """Test that the pass replaces an opted-in node via its expander.

        Given:
            A DISJOIN AST whose operator is flagged GIQL_EXPAND and an expander
            registered for (GenericTarget, GIQLDisjoin).
        When:
            Running the ExpandOperators pass.
        Then:
            The DISJOIN node should be replaced by the expander's output.
        """
        # Arrange
        clean_registry.register(GenericTarget(), GIQLDisjoin, _record("expanded"))
        tables = _tables()
        ast = _prepare("SELECT * FROM DISJOIN(variants)", tables)
        pass_ = ExpandOperators(GenericTarget(), tables, clean_registry)

        # Act
        with _opted_in(GIQLDisjoin):
            result = pass_.transform(ast)

        # Assert
        assert not list(result.find_all(GIQLDisjoin))
        assert result.find(exp.Literal).this == "expanded"

    def test_transform_dispatches_to_fake_custom_target(self, clean_registry):
        """Test that the pass dispatches to a user-defined custom target.

        Given:
            A fake user target and an expander registered for it via @register,
            plus an opted-in DISJOIN operator.
        When:
            Running ExpandOperators with the fake target as active.
        Then:
            The pass should dispatch to the fake target's expander (the public
            extension hook works end to end).
        """

        # Arrange
        @register(_FakeTarget, GIQLDisjoin)
        def _fake_expander(node, ctx):
            assert ctx.target == _FakeTarget()
            return exp.Literal.string("fake-target")

        tables = _tables()
        ast = _prepare("SELECT * FROM DISJOIN(variants)", tables)
        pass_ = ExpandOperators(_FakeTarget(), tables, clean_registry)

        # Act
        with _opted_in(GIQLDisjoin):
            result = pass_.transform(ast)

        # Assert
        assert result.find(exp.Literal).this == "fake-target"

    def test_transform_passes_resolution_and_target_in_context(self, clean_registry):
        """Test that the expander receives the pass-1 metadata and active target.

        Given:
            An opted-in DISJOIN with a registered expander that captures its ctx.
        When:
            Running the pass with a DataFusionTarget.
        Then:
            The context should carry the operator's resolution metadata, the
            active target, and the registered tables.
        """
        # Arrange
        captured = {}

        def _expander(node, ctx):
            captured["ctx"] = ctx
            return exp.Literal.string("ok")

        clean_registry.register(DataFusionTarget(), GIQLDisjoin, _expander)
        tables = _tables()
        ast = _prepare("SELECT * FROM DISJOIN(variants)", tables)
        pass_ = ExpandOperators(DataFusionTarget(), tables, clean_registry)

        # Act
        with _opted_in(GIQLDisjoin):
            pass_.transform(ast)

        # Assert
        ctx = captured["ctx"]
        assert ctx.target == DataFusionTarget()
        assert ctx.tables is tables
        assert ctx.resolution is not None
        assert ctx.resolution.operator == "GIQLDisjoin"

    def test_transform_skips_unflagged_operator(self, clean_registry):
        """Test that an unflagged operator is left untouched even when registered.

        Given:
            An expander registered for (GenericTarget, GIQLDisjoin) but the
            operator's GIQL_EXPAND flag left at its default False.
        When:
            Running the pass.
        Then:
            The DISJOIN node should remain in the tree (gate requires both).
        """
        # Arrange
        clean_registry.register(GenericTarget(), GIQLDisjoin, _record("expanded"))
        tables = _tables()
        ast = _prepare("SELECT * FROM DISJOIN(variants)", tables)
        pass_ = ExpandOperators(GenericTarget(), tables, clean_registry)

        # Act (GIQL_EXPAND is False by default — no opt-in context)
        result = pass_.transform(ast)

        # Assert
        assert list(result.find_all(GIQLDisjoin))

    def test_transform_skips_flagged_operator_with_no_expander(self, clean_registry):
        """Test that a flagged operator with no expander is left untouched.

        Given:
            An opted-in DISJOIN but an empty registry (no expander resolves).
        When:
            Running the pass.
        Then:
            The DISJOIN node should remain (gate requires a registered expander).
        """
        # Arrange
        tables = _tables()
        ast = _prepare("SELECT * FROM DISJOIN(variants)", tables)
        pass_ = ExpandOperators(GenericTarget(), tables, clean_registry)

        # Act
        with _opted_in(GIQLDisjoin):
            result = pass_.transform(ast)

        # Assert
        assert list(result.find_all(GIQLDisjoin))


class TestNoOpWhenInert:
    """The pass touches nothing while no operator opts in and the registry is empty."""

    def test_expand_operators_is_identity_when_registry_empty(self):
        """Test that the standalone pass returns the AST unchanged when inert.

        Given:
            A DISJOIN AST, no operator flagged, and an empty registry.
        When:
            Running expand_operators directly.
        Then:
            The same AST is returned with the DISJOIN node intact.
        """
        # Arrange
        tables = _tables()
        ast = _prepare("SELECT * FROM DISJOIN(variants)", tables)

        # Act
        result = expand_operators(ast, GenericTarget(), tables, ExpanderRegistry())

        # Assert
        assert result is ast
        assert list(result.find_all(GIQLDisjoin))

    def test_transpile_sql_unchanged_with_pass_inert(self):
        """Test that transpile output is byte-identical with the pass inert.

        Given:
            A DISJOIN query, the default (empty) registry, and no operator flagged.
        When:
            Transpiling with the wired-in pass versus a pass-bypassed reference.
        Then:
            The SQL should match exactly and carry no expander alias prefix.
        """
        # Arrange
        query = "SELECT * FROM DISJOIN(variants)"
        tables = _tables()
        ast = _prepare(query, tables)
        from giql.generators import BaseGIQLGenerator

        expected = BaseGIQLGenerator(tables=tables).generate(ast)

        # Act
        actual = transpile(query, tables=[Table("variants")])

        # Assert
        assert actual == expected
        assert EXPAND_ALIAS_PREFIX not in actual


def _tables(names=("variants",)) -> Tables:
    """Build a Tables container registering each name with default encoding."""
    container = Tables()
    for name in names:
        container.register(name, Table(name))
    return container


def _prepare(query: str, tables: Tables) -> exp.Expression:
    """Parse *query* and run pass 1 (resolution) over the resulting AST."""
    ast = parse_one(query, dialect=GIQLDialect)
    return resolve_operator_refs(ast, tables)


class _opted_in:
    """Context manager opting an operator class into GIQL_EXPAND for a test."""

    def __init__(self, operator: type) -> None:
        self._operator = operator
        self._prior = operator.__dict__.get("GIQL_EXPAND", False)

    def __enter__(self):
        self._operator.GIQL_EXPAND = True
        return self._operator

    def __exit__(self, *exc):
        self._operator.GIQL_EXPAND = self._prior
        return False
