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
from giql.expressions import GIQLMerge
from giql.expressions import Intersects
from giql.generators import BaseGIQLGenerator
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


#: The registry contents at import — the built-in expanders registered by
#: ``giql.expanders`` for the already-migrated operators. The leak guards and
#: ``clean_registry`` treat this as the baseline rather than an empty registry,
#: so the real built-in registrations survive isolating fixtures and a leaking
#: test is still caught against the true baseline.
import giql.expanders  # noqa: F401, E402

_REGISTRY_BASELINE = REGISTRY.snapshot()


@pytest.fixture
def clean_registry():
    """Isolate the process-wide REGISTRY, restoring its baseline afterward.

    Saves the import-time baseline (the built-in expanders), empties the registry
    so a test sees only what it registers, and restores the baseline on the way
    out through the public ``snapshot()``/``restore()`` seam — so a test that
    registers a stand-in expander cannot leak it, and the built-in registrations
    survive this fixture's isolation.
    """
    saved = REGISTRY.snapshot()
    REGISTRY.clear()
    yield REGISTRY
    REGISTRY.restore(saved)


@pytest.fixture(autouse=True)
def _registry_leak_guard():
    """Assert the process-wide REGISTRY matches its baseline at each boundary.

    A leak guard: the registry holds the built-in expanders at import and must
    return to exactly that baseline after every test, so a test that registers
    without cleaning up (a leak that would silently change dispatch for a later
    test) fails loudly. Tests that mutate the process-wide REGISTRY do so through
    ``clean_registry``, which restores the baseline on the way out; this guard
    catches anything that bypasses it.
    """
    assert REGISTRY.snapshot() == _REGISTRY_BASELINE, (
        "REGISTRY differed from its baseline entering a test"
    )
    yield
    assert REGISTRY.snapshot() == _REGISTRY_BASELINE, (
        "a test leaked a registration into REGISTRY"
    )


@pytest.fixture(autouse=True)
def _expand_flag_leak_guard():
    """Assert every operator's GIQL_EXPAND is restored at each test boundary.

    The symmetric partner of the registry leak guard: each operator class ships
    a shipped GIQL_EXPAND default (``True`` for a migrated operator,
    ``False`` otherwise), and a test that flips one via ``_opted_in`` must restore
    it. A leaked flip would silently change the pass for a later test, so this
    catches anything that bypasses the exception-safe ``_opted_in`` manager by
    comparing against each operator's shipped default rather than a blanket False.
    """
    for op in _OPERATOR_CLASSES:
        assert op.__dict__.get("GIQL_EXPAND") is _SHIPPED_EXPAND_FLAGS[op], (
            f"{op.__name__}.GIQL_EXPAND leaked into a test from a prior one"
        )
    yield
    for op in _OPERATOR_CLASSES:
        assert op.__dict__.get("GIQL_EXPAND") is _SHIPPED_EXPAND_FLAGS[op], (
            f"a test leaked a GIQL_EXPAND opt-in on {op.__name__}"
        )


class _CountingExpander:
    """An expander recording each call so a walk's dispatch count can be pinned."""

    def __init__(self, label="counted"):
        self.calls = []
        self._label = label

    def __call__(self, node, ctx):
        self.calls.append((node, ctx))
        return exp.Literal.string(self._label)


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

    def test_register_rejects_object_with_non_callable_expand(self):
        """Test that an object whose expand attribute is not callable is rejected.

        Given:
            An object exposing an 'expand' attribute that is data, not a method,
            so it passes the runtime-checkable protocol's attribute check.
        When:
            Registering it.
        Then:
            It should raise TypeError (the registry requires a callable expand,
            not merely the attribute's presence).
        """

        # Arrange
        class _BadExpander:
            expand = "not callable"

        registry = ExpanderRegistry()

        # Act & assert
        with pytest.raises(TypeError):
            registry.register(GenericTarget(), GIQLDisjoin, _BadExpander())

    def test_snapshot_should_not_observe_later_registrations(self):
        """Test that a snapshot does not observe registrations made after it.

        Given:
            A registry with one entry, captured by snapshot.
        When:
            A second entry is registered after the snapshot is taken.
        Then:
            The snapshot should still hold only the first entry (it is a copy,
            not a live view).
        """
        # Arrange
        registry = ExpanderRegistry()
        registry.register(DuckDBTarget(), GIQLDisjoin, _record("first"))

        # Act
        saved = registry.snapshot()
        registry.register(GenericTarget(), Intersects, _record("second"))

        # Assert
        assert (DuckDBTarget(), GIQLDisjoin) in saved
        assert (GenericTarget(), Intersects) not in saved

    def test_restore_should_replace_entries_with_snapshot_contents(self):
        """Test that restore returns the registry to a captured snapshot.

        Given:
            A snapshot of a registry with one entry, after which the registry is
            cleared and a different entry registered.
        When:
            Restoring the snapshot.
        Then:
            The original entry should resolve again and the post-snapshot entry
            should be gone.
        """
        # Arrange
        registry = ExpanderRegistry()
        registry.register(DuckDBTarget(), GIQLDisjoin, _record("original"))
        saved = registry.snapshot()
        registry.clear()
        registry.register(GenericTarget(), Intersects, _record("transient"))

        # Act
        registry.restore(saved)

        # Assert
        assert (DuckDBTarget(), GIQLDisjoin) in registry
        assert (GenericTarget(), Intersects) not in registry


class TestExpanderRegistryFallbackGaps:
    """Edge cases of the registry fallback and op-scoped keying."""

    def test_non_generic_target_with_no_entries_resolves_none(self):
        """Test that a target with neither a specific nor generic entry yields None.

        Given:
            A registry with no entry for the target and no generic entry.
        When:
            Resolving (DataFusionTarget(), op).
        Then:
            It resolves to None (legacy emitter fallback), having exhausted the
            chain.
        """
        # Arrange
        registry = ExpanderRegistry()
        registry.register(DuckDBTarget(), GIQLDisjoin, _record("duckdb"))

        # Act
        resolved = registry.resolve(DataFusionTarget(), GIQLDisjoin)

        # Assert
        assert resolved is None

    def test_generic_entry_is_operator_scoped(self):
        """Test that a generic entry for one operator does not satisfy another.

        Given:
            A generic expander registered for operator X (GIQLDisjoin) only.
        When:
            Resolving operator Y (GIQLMerge) for a target.
        Then:
            It resolves to None — the generic fallback is keyed per operator, not
            a catch-all for the target.
        """
        # Arrange
        registry = ExpanderRegistry()
        registry.register(GenericTarget(), GIQLDisjoin, _record("generic"))

        # Act
        resolved = registry.resolve(DuckDBTarget(), GIQLMerge)

        # Assert
        assert resolved is None

    def test_contains_checks_exact_key_not_fallback(self):
        """Test that __contains__ tests the exact key, ignoring the generic chain.

        Given:
            A registry with a generic entry for op only.
        When:
            Testing membership of an exact (DuckDBTarget, op) key.
        Then:
            The exact target key reports absent while the generic key reports
            present (membership does not walk the fallback chain).
        """
        # Arrange
        registry = ExpanderRegistry()
        registry.register(GenericTarget(), GIQLDisjoin, _record("generic"))

        # Act & assert
        assert (GenericTarget(), GIQLDisjoin) in registry
        assert (DuckDBTarget(), GIQLDisjoin) not in registry

    def test_unregister_drops_entry_and_falls_through(self):
        """Test that unregister drops an entry so the key resolves to None after.

        Given:
            A registry with one (DuckDBTarget, op) entry.
        When:
            Unregistering that key.
        Then:
            Resolving it returns None and the key reports absent (public teardown
            seam).
        """
        # Arrange
        registry = ExpanderRegistry()
        registry.register(DuckDBTarget(), GIQLDisjoin, _record("duckdb"))

        # Act
        registry.unregister(DuckDBTarget(), GIQLDisjoin)

        # Assert
        assert registry.resolve(DuckDBTarget(), GIQLDisjoin) is None
        assert (DuckDBTarget(), GIQLDisjoin) not in registry

    def test_unregister_absent_key_is_noop(self):
        """Test that unregistering an absent key does not raise.

        Given:
            An empty registry.
        When:
            Unregistering a key that was never registered.
        Then:
            It returns without raising (idempotent teardown).
        """
        # Arrange
        registry = ExpanderRegistry()

        # Act & assert (no raise)
        registry.unregister(DuckDBTarget(), GIQLDisjoin)
        assert (DuckDBTarget(), GIQLDisjoin) not in registry

    def test_public_teardown_seam_resets_process_registry(self, clean_registry):
        """Test that the public seam tears a custom registration off REGISTRY.

        Given:
            A custom expander registered on the process-wide REGISTRY via the
            decorator (the public extension hook).
        When:
            Tearing it down through the public seam — unregister, then clear —
            without touching private state.
        Then:
            The key resolves to None and REGISTRY is empty again, so the leak
            guard's post-test assertion holds (the seam is a supported reset).
        """

        # Arrange
        @register(DuckDBTarget, GIQLDisjoin)
        def _expander(node, ctx):
            return exp.column("x")

        assert (DuckDBTarget(), GIQLDisjoin) in REGISTRY

        # Act
        REGISTRY.unregister(DuckDBTarget(), GIQLDisjoin)
        REGISTRY.clear()

        # Assert
        assert REGISTRY.resolve(DuckDBTarget(), GIQLDisjoin) is None
        assert (DuckDBTarget(), GIQLDisjoin) not in REGISTRY

    def test_has_override_should_return_true_when_exact_nongeneric_entry_registered(
        self,
    ):
        """Test that has_override is True for an exact non-generic entry.

        Given:
            A registry with an exact (DuckDBTarget, op) entry registered.
        When:
            Querying has_override for that exact key.
        Then:
            It should return True (a non-generic exact entry is an override).
        """
        # Arrange
        registry = ExpanderRegistry()
        registry.register(DuckDBTarget(), GIQLDisjoin, _record("duckdb"))

        # Act & assert
        assert registry.has_override(DuckDBTarget(), GIQLDisjoin) is True

    def test_has_override_should_return_false_for_generic_fallback_entry(self):
        """Test that has_override ignores the generic fallback entry.

        Given:
            A registry with only a (GenericTarget, op) entry registered.
        When:
            Querying has_override for a non-generic target's key.
        Then:
            It should return False (the portable generic fallback is not an
            override, even though resolve() would route to it).
        """
        # Arrange
        registry = ExpanderRegistry()
        registry.register(GenericTarget(), GIQLDisjoin, _record("generic"))

        # Act & assert
        assert registry.has_override(DuckDBTarget(), GIQLDisjoin) is False
        # And a generic target queried against its own entry is not an override.
        assert registry.has_override(GenericTarget(), GIQLDisjoin) is False

    def test_has_override_should_return_false_when_unregistered(self):
        """Test that has_override is False when nothing is registered.

        Given:
            An empty registry.
        When:
            Querying has_override for any key.
        Then:
            It should return False (no entry, so no override).
        """
        # Arrange
        registry = ExpanderRegistry()

        # Act & assert
        assert registry.has_override(DuckDBTarget(), GIQLDisjoin) is False


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

    def test_register_class_and_instance_land_on_same_entry(self, clean_registry):
        """Test that the class and an instance form key the same registry entry.

        Given:
            A registration via the class form @register(DuckDBTarget, op).
        When:
            Registering again via the instance form @register(DuckDBTarget(), op).
        Then:
            The second lands on the same entry and overrides the first (the
            decorator instantiates the class to a value-equal target key).
        """

        # Arrange
        @register(DuckDBTarget, GIQLDisjoin)
        def _first(node, ctx):
            return exp.Literal.string("first")

        # Act
        @register(DuckDBTarget(), GIQLDisjoin)
        def _second(node, ctx):
            return exp.Literal.string("second")

        # Assert
        assert clean_registry.resolve(DuckDBTarget(), GIQLDisjoin) is _second

    def test_register_stacks_one_expander_for_two_targets(self, clean_registry):
        """Test that stacking the decorator registers one expander for two targets.

        Given:
            An expander decorated with @register for both DuckDB and DataFusion.
        When:
            Resolving each target's key.
        Then:
            Both resolve to the same single expander object.
        """

        # Arrange & act
        @register(DuckDBTarget, GIQLDisjoin)
        @register(DataFusionTarget, GIQLDisjoin)
        def _shared(node, ctx):
            return exp.Literal.string("shared")

        # Assert
        assert clean_registry.resolve(DuckDBTarget(), GIQLDisjoin) is _shared
        assert clean_registry.resolve(DataFusionTarget(), GIQLDisjoin) is _shared

    def test_register_overrides_through_decorator(self, clean_registry):
        """Test that a second @register for one key overrides the first.

        Given:
            Two expanders decorated for the same (target, op) key in order.
        When:
            Resolving the key.
        Then:
            It resolves to the second (last-write-wins through the decorator).
        """

        # Arrange
        @register(GenericTarget, GIQLDisjoin)
        def _first(node, ctx):
            return exp.Literal.string("first")

        # Act
        @register(GenericTarget, GIQLDisjoin)
        def _second(node, ctx):
            return exp.Literal.string("second")

        # Assert
        assert clean_registry.resolve(GenericTarget(), GIQLDisjoin) is _second

    def test_register_accepts_object_form_through_decorator(self, clean_registry):
        """Test that the decorator accepts an OperatorExpander object.

        Given:
            An OperatorExpander object passed to the @register decorator.
        When:
            Resolving the key and invoking it.
        Then:
            It resolves to a callable delegating to the object's expand method.
        """

        # Arrange
        class _Obj:
            def expand(self, node, ctx):
                return exp.Literal.string("obj-form")

        obj = _Obj()

        # Act
        register(GenericTarget, GIQLDisjoin)(obj)
        resolved = clean_registry.resolve(GenericTarget(), GIQLDisjoin)

        # Assert
        assert resolved(None, None) == exp.Literal.string("obj-form")

    def test_register_rejects_non_callable_through_decorator(self, clean_registry):
        """Test that decorating a non-callable, non-expander raises TypeError.

        Given:
            A plain object that is neither callable nor an OperatorExpander.
        When:
            Passing it through the @register decorator.
        Then:
            It raises TypeError (the same guard as direct registration, proven
            through the public decorator).
        """
        # Arrange & act & assert
        with pytest.raises(TypeError):
            register(GenericTarget, GIQLDisjoin)(object())

    def test_register_does_not_leak_across_targets(self, clean_registry):
        """Test that registering for one target leaves another target unaffected.

        Given:
            An expander registered for (DuckDBTarget, op) only.
        When:
            Resolving (DataFusionTarget, op), with no generic fallback present.
        Then:
            It resolves to None (no cross-target leakage).
        """

        # Arrange & act
        @register(DuckDBTarget, GIQLDisjoin)
        def _expander(node, ctx):
            return exp.Literal.string("duck")

        # Assert
        assert clean_registry.resolve(DataFusionTarget(), GIQLDisjoin) is None

    def test_register_returns_callable_decorated_function(self, clean_registry):
        """Test that the bare decorated function stays callable standalone.

        Given:
            A function decorated with @register.
        When:
            Calling the decorated function directly (outside the pass).
        Then:
            It runs and returns its expansion, unchanged by decoration.
        """

        # Arrange & act
        @register(GenericTarget, GIQLDisjoin)
        def _expander(node, ctx):
            return exp.Literal.string("callable")

        # Assert
        assert _expander(None, None) == exp.Literal.string("callable")


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

    def test_resolution_none_is_tolerated(self):
        """Test that a None resolution is a tolerated first-class context state.

        Given:
            A context constructed with resolution=None.
        When:
            Reading its resolution attribute.
        Then:
            It is None without error (an un-annotated node is representable).
        """
        # Arrange
        ctx = ExpansionContext(exp.true(), None, GenericTarget(), Tables())

        # Act & assert
        assert ctx.resolution is None

    def test_tables_attribute_is_the_passed_container(self):
        """Test that the context exposes the tables container it was built with.

        Given:
            A context built with a specific Tables container.
        When:
            Reading its tables attribute.
        Then:
            It is that same container (identity), threaded for the expander.
        """
        # Arrange
        tables = _tables()
        ctx = ExpansionContext(exp.true(), None, GenericTarget(), tables)

        # Act & assert
        assert ctx.tables is tables

    def test_alias_uses_reserved_prefix_format(self):
        """Test that minted aliases carry the reserved expander prefix.

        Given:
            A fresh context.
        When:
            Minting an alias.
        Then:
            It starts with the EXPAND_ALIAS_PREFIX reserved namespace.
        """
        # Arrange
        ctx = ExpansionContext(exp.true(), None, GenericTarget(), Tables())

        # Act
        alias = ctx.alias()

        # Assert
        assert alias.startswith(EXPAND_ALIAS_PREFIX)

    def test_alias_is_monotonic_within_one_context(self):
        """Test that successive aliases in one context are all distinct.

        Given:
            A single context.
        When:
            Minting several aliases.
        Then:
            They are all distinct (a monotonic sequence, no repeats).
        """
        # Arrange
        ctx = ExpansionContext(exp.true(), None, GenericTarget(), Tables())

        # Act
        names = [ctx.alias() for _ in range(5)]

        # Assert
        assert len(set(names)) == 5

    def test_alias_continues_across_contexts_sharing_a_sequence(self):
        """Test that two contexts sharing one sequence mint non-colliding aliases.

        Given:
            Two contexts built with the same shared name_sequence callable.
        When:
            Each mints an alias.
        Then:
            The aliases differ (the sequence threads continuity across the
            per-node contexts of one run, preventing __giql_x_ collisions).
        """
        # Arrange
        from sqlglot.helper import name_sequence

        shared = name_sequence(EXPAND_ALIAS_PREFIX)
        ctx_a = ExpansionContext(exp.true(), None, GenericTarget(), Tables(), shared)
        ctx_b = ExpansionContext(exp.true(), None, GenericTarget(), Tables(), shared)

        # Act
        a = ctx_a.alias()
        b = ctx_b.alias()

        # Assert
        assert a != b

    def test_capabilities_delegates_on_second_target(self):
        """Test that capabilities delegates to a different target's capability set.

        Given:
            A context built with a DataFusionTarget.
        When:
            Reading its capabilities property.
        Then:
            It is the DataFusion target's Capabilities (delegation is per-target,
            not hardcoded to one engine).
        """
        # Arrange
        target = DataFusionTarget()
        ctx = ExpansionContext(exp.true(), None, target, Tables())

        # Act
        caps = ctx.capabilities

        # Assert
        assert caps is target.capabilities
        assert caps.supports_lateral is False


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
        """Test that an opted-out operator is left untouched even when registered.

        Given:
            A migrated operator with an expander registered for it, but its
            GIQL_EXPAND flag opted out for the test — so the *same* operator is
            registered, queried, and opted out, isolating the per-type gate.
        When:
            Running the pass.
        Then:
            The operator node should remain in the tree (the gate requires the
            flag, so opting it out alone holds expansion off).
        """
        # Arrange
        operator = _A_MIGRATED_OPERATOR
        clean_registry.register(GenericTarget(), operator, _record("expanded"))
        ast, tables = _prepare_operator(operator)
        pass_ = ExpandOperators(GenericTarget(), tables, clean_registry)

        # Act
        with _opted_out(operator):
            result = pass_.transform(ast)

        # Assert
        assert list(result.find_all(operator))

    def test_transform_expands_flagged_operator(self, clean_registry):
        """Test that the same operator expands once flagged (the paired positive).

        Given:
            The same migrated operator and registered expander as the opt-out
            control, but with its GIQL_EXPAND flag opted in.
        When:
            Running the pass.
        Then:
            The operator node should be replaced by the expander's output — so the
            contrast with the opt-out control pins the per-type gate as
            load-bearing, not vacuous.
        """
        # Arrange
        operator = _A_MIGRATED_OPERATOR
        clean_registry.register(GenericTarget(), operator, _record("expanded"))
        ast, tables = _prepare_operator(operator)
        pass_ = ExpandOperators(GenericTarget(), tables, clean_registry)

        # Act
        with _opted_in(operator):
            result = pass_.transform(ast)

        # Assert
        assert not list(result.find_all(operator))
        assert result.find(exp.Literal).this == "expanded"

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
        """Test that transpile output is byte-identical for an unmigrated operator.

        Given:
            A CLUSTER query (an operator not migrated onto the pass in any wave-3
            branch, so its GIQL_EXPAND is False and no expander resolves), with
            the default registry.
        When:
            Running the ExpandOperators pass (default REGISTRY) over the resolved
            AST and serializing both the original and the pass-run AST.
        Then:
            The pass leaves the operator node in place, the serialized SQL is
            byte-identical, and no expander alias prefix appears — the pass is
            inert for any operator that has not been migrated.
        """
        # Arrange
        query = "SELECT *, CLUSTER(interval) AS cluster_id FROM peaks"
        tables = _tables(("peaks",))
        ast = _prepare(query, tables)
        from giql.generators import BaseGIQLGenerator

        before = BaseGIQLGenerator(tables=tables).generate(ast)
        before_ops = len(list(ast.find_all(GIQLCluster)))

        # Act — the wired-in pass over the default REGISTRY must be a no-op here.
        result = expand_operators(ast, GenericTarget(), tables)
        after = BaseGIQLGenerator(tables=tables).generate(result)

        # Assert
        assert after == before
        assert len(list(result.find_all(GIQLCluster))) == before_ops
        assert EXPAND_ALIAS_PREFIX not in after


# The nine GIQL operator expression classes the ExpandOperators pass inspects.
# Each migrated operator ships opted in (GIQL_EXPAND=True) alongside its
# registered expander; the rest ship opted out (False) and fall through to the
# legacy emitter.
from giql.expressions import Contains  # noqa: E402
from giql.expressions import GIQLCluster  # noqa: E402
from giql.expressions import GIQLDistance  # noqa: E402
from giql.expressions import GIQLNearest  # noqa: E402
from giql.expressions import SpatialSetPredicate  # noqa: E402
from giql.expressions import Within  # noqa: E402

_OPERATOR_CLASSES = (
    Intersects,
    Contains,
    Within,
    SpatialSetPredicate,
    GIQLDistance,
    GIQLNearest,
    GIQLDisjoin,
    GIQLCluster,
    GIQLMerge,
)

#: Each operator's shipped GIQL_EXPAND default, captured from its own class dict
#: at import. A migrated operator ships ``True``; the rest ship ``False`` until
#: their migrations land. The flag leak guard restores to these shipped values
#: rather than a blanket ``False``.
# Pin the leak guard's assumption that every operator declares GIQL_EXPAND on its
# own class (so __dict__.get() reads the operator's shipped value rather than an
# inherited one). A future operator that omits its own flag would make the guard
# read ``None`` and silently lose its leak coverage; this one-time check closes
# that hole.
for _op in _OPERATOR_CLASSES:
    assert "GIQL_EXPAND" in _op.__dict__, (
        f"{_op.__name__} must declare GIQL_EXPAND on its own class"
    )
_SHIPPED_EXPAND_FLAGS = {op: op.__dict__.get("GIQL_EXPAND") for op in _OPERATOR_CLASSES}


#: Operators migrated onto the ExpandOperators pass — they ship GIQL_EXPAND=True.
_MIGRATED_OPERATORS = tuple(
    op for op in _OPERATOR_CLASSES if op.__dict__.get("GIQL_EXPAND") is True
)
assert _MIGRATED_OPERATORS, "expected at least one migrated operator"
#: An arbitrary migrated operator the operator-agnostic control tests target.
_A_MIGRATED_OPERATOR = _MIGRATED_OPERATORS[0]
#: Operators not yet migrated — they ship GIQL_EXPAND=False.
_UNMIGRATED_OPERATORS = tuple(
    op for op in _OPERATOR_CLASSES if op not in _MIGRATED_OPERATORS
)
assert _UNMIGRATED_OPERATORS, "expected at least one unmigrated operator"


#: A minimal GIQL query producing one node of each operator class, keyed by the
#: class. Lets the control tests build a node for *any* operator — whichever the
#: branch ships migrated/unmigrated — so they stay operator-agnostic rather than
#: hard-wiring a particular operator. The second element is the table names the
#: query references (registered before pass 1).
_OPERATOR_QUERIES: dict[type, tuple[str, tuple[str, ...]]] = {
    Intersects: (
        "SELECT * FROM variants WHERE interval INTERSECTS 'chr1:1000-2000'",
        ("variants",),
    ),
    Contains: (
        "SELECT * FROM variants WHERE interval CONTAINS 'chr1:1500-1600'",
        ("variants",),
    ),
    Within: (
        "SELECT * FROM variants WHERE interval WITHIN 'chr1:1000-5000'",
        ("variants",),
    ),
    SpatialSetPredicate: (
        "SELECT * FROM variants "
        "WHERE interval INTERSECTS ANY('chr1:1000-2000', 'chr1:5000-6000')",
        ("variants",),
    ),
    GIQLDistance: (
        "SELECT DISTANCE(a.interval, b.interval) AS d "
        "FROM features_a a CROSS JOIN features_b b",
        ("features_a", "features_b"),
    ),
    GIQLNearest: (
        "SELECT * FROM peaks CROSS JOIN LATERAL NEAREST(genes, k := 3)",
        ("peaks", "genes"),
    ),
    GIQLDisjoin: ("SELECT * FROM DISJOIN(variants)", ("variants",)),
    GIQLCluster: (
        "SELECT *, CLUSTER(interval) AS cluster_id FROM peaks",
        ("peaks",),
    ),
    GIQLMerge: ("SELECT MERGE(interval) AS m FROM peaks", ("peaks",)),
}


def _prepare_operator(operator: type) -> tuple[exp.Expression, Tables]:
    """Build a pass-1-resolved AST containing one node of *operator*.

    Operator-agnostic: looks the query up in :data:`_OPERATOR_QUERIES` so a
    control test can exercise whichever operator a branch ships migrated or
    unmigrated, rather than hard-wiring one operator class. Returns the resolved
    AST and the Tables container it was resolved against (the same container must
    be threaded into the pass).
    """
    query, names = _OPERATOR_QUERIES[operator]
    tables = _tables(names)
    return _prepare(query, tables), tables


class TestOperatorOptOut:
    """Migrated operators opt into the pass; the rest still ship opted out."""

    @pytest.mark.parametrize(
        "operator", _UNMIGRATED_OPERATORS, ids=lambda c: c.__name__
    )
    def test_operator_class_ships_expand_disabled(self, operator):
        """Test that each unmigrated operator class ships GIQL_EXPAND=False.

        Given:
            A GIQL operator expression class that has not been migrated onto the
            ExpandOperators pass.
        When:
            Reading its GIQL_EXPAND class attribute.
        Then:
            It should be False (the operator still uses the legacy emitter).
        """
        # Arrange & act
        flag = operator.GIQL_EXPAND

        # Assert
        assert flag is False

    @pytest.mark.parametrize(
        "operator", _MIGRATED_OPERATORS, ids=lambda c: c.__name__
    )
    def test_operator_class_ships_expand_enabled(self, operator):
        """Test that each migrated operator class ships GIQL_EXPAND=True.

        Given:
            A GIQL operator expression class migrated onto the ExpandOperators
            pass.
        When:
            Reading its GIQL_EXPAND class attribute.
        Then:
            It should be True (the operator expands through its registered
            expander instead of the deleted legacy emitter).
        """
        # Arrange & act
        flag = operator.GIQL_EXPAND

        # Assert
        assert flag is True


class TestMigratedOperatorsRegistered:
    """Every migrated operator resolves a built-in expander in the process REGISTRY."""

    @pytest.mark.parametrize(
        "operator", _MIGRATED_OPERATORS, ids=lambda c: c.__name__
    )
    def test_migrated_operator_resolves_in_process_registry(self, operator):
        """Test that each migrated operator resolves an expander in REGISTRY.

        Given:
            A GIQL operator class that ships GIQL_EXPAND=True (migrated onto the
            ExpandOperators pass).
        When:
            Resolving it against the import-populated process-wide REGISTRY for
            the generic target.
        Then:
            A built-in expander should resolve — a migrated operator always has a
            registered expander, so the pass never leaves it on a deleted emitter.
        """
        # Arrange & act
        resolved = REGISTRY.resolve(GenericTarget(), operator)

        # Assert
        assert resolved is not None


class TestOptedInRestoresFlag:
    """The _opted_in helper restores GIQL_EXPAND even when its body raises."""

    def test_opted_in_restores_flag_after_exception(self):
        """Test that _opted_in restores GIQL_EXPAND when the body raises.

        Given:
            An operator class at its default GIQL_EXPAND=False.
        When:
            Its _opted_in body raises an exception.
        Then:
            The flag should be restored to False (the manager is exception-safe,
            so a raising expansion test cannot leak an opt-in into a later test).
        """
        # Arrange
        assert GIQLMerge.GIQL_EXPAND is False

        # Act
        with pytest.raises(RuntimeError):
            with _opted_in(GIQLMerge):
                assert GIQLMerge.GIQL_EXPAND is True
                raise RuntimeError("boom")

        # Assert
        assert GIQLMerge.GIQL_EXPAND is False


class TestIEJoinRegistryDeferral:
    """The duckdb IEJoin path defers to a target-specific Intersects expander (#141).

    Resolves Finding 2: the IEJoin early return used to emit before the
    ExpandOperators pass, so a flagged operator on an IEJoin-eligible query was
    never expanded. Now a *target-specific* ``(DuckDBTarget, Intersects)``
    registry entry overrides the built-in join strategy entirely (the public
    extension hook), while the default duckdb path — with no such override —
    still emits the built-in IEJoin SQL.
    """

    def test_iejoin_query_expands_target_override_expander(self, clean_registry):
        """Test that a target-specific Intersects override fires on an IEJoin query.

        Given:
            A column-to-column INTERSECTS join eligible for the duckdb IEJoin
            path, with a (DuckDBTarget, Intersects) expander registered.
        When:
            Transpiling with dialect='duckdb'.
        Then:
            The override expander's sentinel should appear — the IEJoin path
            defers to the registry rather than short-circuiting expansion.
        """
        # Arrange
        clean_registry.register(
            DuckDBTarget(), Intersects, lambda n, c: exp.column("__giql_iejoin_sentinel")
        )
        query = (
            "SELECT a.start FROM peaks a "
            "JOIN genes b ON a.interval INTERSECTS b.interval"
        )

        # Act
        sql = transpile(query, tables=["peaks", "genes"], dialect="duckdb")

        # Assert
        assert "__giql_iejoin_sentinel" in sql
        assert "SET VARIABLE __giql_iejoin_" not in sql

    def test_iejoin_query_emits_builtin_iejoin_without_override(self):
        """Test that the default duckdb path emits the built-in IEJoin SQL.

        Given:
            The same IEJoin-eligible duckdb query and no target-specific
            Intersects override registered (only the built-in generic expander).
        When:
            Transpiling with dialect='duckdb'.
        Then:
            The built-in IEJoin SET VARIABLE SQL is emitted (the generic
            predicate expander does not disable the join strategy).
        """
        # Arrange
        query = (
            "SELECT a.start FROM peaks a "
            "JOIN genes b ON a.interval INTERSECTS b.interval"
        )

        # Act
        sql = transpile(query, tables=["peaks", "genes"], dialect="duckdb")

        # Assert
        assert "SET VARIABLE __giql_iejoin_" in sql


class TestTranspileExpanderDispatch:
    """Dispatch driven through the public transpile() entry point."""

    def test_target_entry_dispatches_for_matching_dialect(self, clean_registry):
        """Test that a (DuckDBTarget, op) entry dispatches for dialect='duckdb'.

        Given:
            An expander registered for (DuckDBTarget, GIQLDisjoin) and the
            operator flagged.
        When:
            Transpiling a DISJOIN query with dialect='duckdb'.
        Then:
            The expander's sentinel should reach the generated SQL.
        """
        # Arrange
        clean_registry.register(
            DuckDBTarget(), GIQLDisjoin, lambda n, c: exp.column("DUCK_SENTINEL")
        )

        # Act
        with _opted_in(GIQLDisjoin):
            sql = transpile(
                "SELECT * FROM DISJOIN(variants)", tables=["variants"], dialect="duckdb"
            )

        # Assert
        assert "DUCK_SENTINEL" in sql

    @pytest.mark.parametrize("dialect", [None, "datafusion"])
    def test_target_entry_falls_through_for_other_dialects(
        self, clean_registry, dialect
    ):
        """Test that a target-specific entry does not fire for other dialects.

        Given:
            An expander registered only for (DuckDBTarget, GIQLDisjoin), flagged.
        When:
            Transpiling with dialect=None or 'datafusion'.
        Then:
            The sentinel should be absent and the legacy SQL emitted (the
            target-specific entry does not resolve for a different target).
        """
        # Arrange
        clean_registry.register(
            DuckDBTarget(), GIQLDisjoin, lambda n, c: exp.column("DUCK_SENTINEL")
        )

        # Act
        with _opted_in(GIQLDisjoin):
            sql = transpile(
                "SELECT * FROM DISJOIN(variants)", tables=["variants"], dialect=dialect
            )

        # Assert
        assert "DUCK_SENTINEL" not in sql

    @pytest.mark.parametrize("dialect", [None, "duckdb", "datafusion"])
    def test_generic_entry_covers_all_dialects(self, clean_registry, dialect):
        """Test that a (GenericTarget, op) entry dispatches via fallback everywhere.

        Given:
            An expander registered only for (GenericTarget, GIQLDisjoin), flagged.
        When:
            Transpiling with each of the three dialects.
        Then:
            The generic expander's sentinel should reach the SQL in every case
            (the fallback chain routes every target to the generic entry).
        """
        # Arrange
        clean_registry.register(
            GenericTarget(), GIQLDisjoin, lambda n, c: exp.column("GENERIC_SENTINEL")
        )

        # Act
        with _opted_in(GIQLDisjoin):
            sql = transpile(
                "SELECT * FROM DISJOIN(variants)", tables=["variants"], dialect=dialect
            )

        # Assert
        assert "GENERIC_SENTINEL" in sql

    def test_target_entry_shadows_generic_per_dialect(self, clean_registry):
        """Test that a target entry shadows the generic entry for its dialect only.

        Given:
            Both a (DuckDBTarget, op) and a (GenericTarget, op) expander, flagged.
        When:
            Transpiling with dialect='duckdb' versus dialect=None.
        Then:
            The duckdb path uses the target sentinel; the generic path uses the
            generic sentinel (per-dialect shadowing through transpile()).
        """
        # Arrange
        clean_registry.register(
            DuckDBTarget(), GIQLDisjoin, lambda n, c: exp.column("DUCK_SENTINEL")
        )
        clean_registry.register(
            GenericTarget(), GIQLDisjoin, lambda n, c: exp.column("GENERIC_SENTINEL")
        )

        # Act
        with _opted_in(GIQLDisjoin):
            duck_sql = transpile(
                "SELECT * FROM DISJOIN(variants)", tables=["variants"], dialect="duckdb"
            )
            generic_sql = transpile(
                "SELECT * FROM DISJOIN(variants)", tables=["variants"], dialect=None
            )

        # Assert
        assert "DUCK_SENTINEL" in duck_sql and "GENERIC_SENTINEL" not in duck_sql
        assert "GENERIC_SENTINEL" in generic_sql and "DUCK_SENTINEL" not in generic_sql

    def test_context_target_matches_resolved_dialect(self, clean_registry):
        """Test that ctx.target equals resolve_target(dialect) at dispatch.

        Given:
            A flagged DISJOIN and an expander capturing its context.
        When:
            Transpiling with dialect='datafusion'.
        Then:
            The captured ctx.target should equal resolve_target('datafusion').
        """
        # Arrange
        from giql.targets import resolve_target

        captured = {}

        def _capture(node, ctx):
            captured["target"] = ctx.target
            return exp.column("ok")

        clean_registry.register(GenericTarget(), GIQLDisjoin, _capture)

        # Act
        with _opted_in(GIQLDisjoin):
            transpile(
                "SELECT * FROM DISJOIN(variants)",
                tables=["variants"],
                dialect="datafusion",
            )

        # Assert
        assert captured["target"] == resolve_target("datafusion")

    def test_throwing_expander_wraps_in_expansion_error(self, clean_registry):
        """Test that an unexpected expander error is wrapped as an Expansion error.

        Given:
            A flagged DISJOIN whose expander raises a non-ValueError.
        When:
            Transpiling.
        Then:
            transpile() should raise ValueError whose message starts with
            'Expansion error: ' (the stage boundary wraps unexpected errors).
        """

        # Arrange
        def _boom(node, ctx):
            raise RuntimeError("kaboom")

        clean_registry.register(GenericTarget(), GIQLDisjoin, _boom)

        # Act & assert
        with _opted_in(GIQLDisjoin):
            with pytest.raises(ValueError, match=r"^Expansion error: "):
                transpile("SELECT * FROM DISJOIN(variants)", tables=["variants"])

    def test_value_error_from_expander_passes_verbatim(self, clean_registry):
        """Test that a ValueError raised by an expander propagates unwrapped.

        Given:
            A flagged DISJOIN whose expander raises a ValueError with a targeted
            diagnostic.
        When:
            Transpiling.
        Then:
            The same ValueError message should propagate verbatim (the boundary
            lets user-facing ValueErrors through untouched).
        """

        # Arrange
        def _raise(node, ctx):
            raise ValueError("a targeted diagnostic")

        clean_registry.register(GenericTarget(), GIQLDisjoin, _raise)

        # Act & assert
        with _opted_in(GIQLDisjoin):
            with pytest.raises(ValueError, match=r"^a targeted diagnostic$"):
                transpile("SELECT * FROM DISJOIN(variants)", tables=["variants"])

    def test_expander_sees_canonicalized_operands(self, clean_registry):
        """Test that expansion runs after canonicalization (pass ordering).

        Given:
            A flagged DISJOIN with a non-canonical (1-based) target table, which
            pass 2 wraps in a __giql_canon_ CTE, and an expander capturing its
            node's resolution at dispatch.
        When:
            Transpiling.
        Then:
            The expander runs after canonicalization, so the surviving query
            already carries the canonical wrapper CTE (expansion sees pass-2
            output, not the raw AST).
        """
        # Arrange
        captured = {}

        def _capture(node, ctx):
            # Snapshot the node's SQL *at dispatch* — canonicalization ran before
            # this pass, so a non-canonical target is already wrapped in a
            # __giql_canon_ CTE the expander sees in its operand.
            captured["node_sql"] = ctx.node.sql(dialect=GIQLDialect)
            captured["resolution"] = ctx.resolution
            # Returning the node unchanged leaves the canonical CTE in the tree.
            return node

        clean_registry.register(GenericTarget(), GIQLDisjoin, _capture)
        one_based = Table("variants", coordinate_system="1based")

        # Act
        with _opted_in(GIQLDisjoin):
            transpile("SELECT * FROM DISJOIN(variants)", tables=[one_based])

        # Assert
        assert captured["resolution"] is not None
        assert "__giql_canon_" in captured["node_sql"]

    def test_replacement_reaches_generator(self, clean_registry):
        """Test that the expander's replacement node is what the generator renders.

        Given:
            A flagged DISJOIN and an expander returning a distinctive sentinel.
        When:
            Transpiling.
        Then:
            The sentinel appears in the final SQL and no DISJOIN remains, so the
            replacement reached the generator (end-to-end).
        """
        # Arrange
        clean_registry.register(
            GenericTarget(), GIQLDisjoin, lambda n, c: exp.column("REACHED_GENERATOR")
        )

        # Act
        with _opted_in(GIQLDisjoin):
            sql = transpile("SELECT * FROM DISJOIN(variants)", tables=["variants"])

        # Assert
        assert "REACHED_GENERATOR" in sql
        assert "DISJOIN" not in sql.upper()

    def test_unknown_dialect_raises_before_dispatch(self, clean_registry):
        """Test that an unknown dialect raises before any expander dispatches.

        Given:
            A flagged DISJOIN and a generic expander that would record a call.
        When:
            Transpiling with an unrecognized dialect.
        Then:
            transpile() raises ValueError for the unknown dialect and the
            expander never runs (dialect resolution precedes dispatch).
        """
        # Arrange
        counting = _CountingExpander()
        clean_registry.register(GenericTarget(), GIQLDisjoin, counting)

        # Act & assert
        with _opted_in(GIQLDisjoin):
            with pytest.raises(ValueError, match="Unknown dialect"):
                transpile(
                    "SELECT * FROM DISJOIN(variants)",
                    tables=["variants"],
                    dialect="postgres",
                )
        assert counting.calls == []


class TestExpandOperatorsWalk:
    """Real-AST walk behavior of the ExpandOperators pass."""

    def test_walk_visits_every_sibling_operator(self, clean_registry):
        """Test that the pass dispatches once per sibling operator and clears them.

        Given:
            A query with two sibling DISJOIN operators, flagged, and a counting
            expander.
        When:
            Running the pass.
        Then:
            The expander is called exactly twice and no DISJOIN node remains.
        """
        # Arrange
        counting = _CountingExpander()
        clean_registry.register(GenericTarget(), GIQLDisjoin, counting)
        tables = _tables(("variants", "genes"))
        ast = _prepare(
            "SELECT * FROM DISJOIN(variants) UNION ALL SELECT * FROM DISJOIN(genes)",
            tables,
        )
        pass_ = ExpandOperators(GenericTarget(), tables, clean_registry)

        # Act
        with _opted_in(GIQLDisjoin):
            result = pass_.transform(ast)

        # Assert
        assert len(counting.calls) == 2
        assert not list(result.find_all(GIQLDisjoin))

    def test_walk_reaches_nested_operator(self, clean_registry):
        """Test that the walk reaches an operator nested in a derived table.

        Given:
            A DISJOIN nested inside a derived table, flagged, with a counting
            expander.
        When:
            Running the pass.
        Then:
            The nested operator is dispatched exactly once.
        """
        # Arrange
        counting = _CountingExpander()
        clean_registry.register(GenericTarget(), GIQLDisjoin, counting)
        tables = _tables()
        ast = _prepare("SELECT * FROM (SELECT * FROM DISJOIN(variants)) t", tables)
        pass_ = ExpandOperators(GenericTarget(), tables, clean_registry)

        # Act
        with _opted_in(GIQLDisjoin):
            result = pass_.transform(ast)

        # Assert
        assert len(counting.calls) == 1
        assert not list(result.find_all(GIQLDisjoin))

    @pytest.mark.parametrize(
        "query",
        [
            "SELECT * FROM peaks WHERE interval INTERSECTS 'chr1:1-100'",
            "SELECT * FROM peaks a JOIN genes b ON a.interval INTERSECTS 'chr1:1-100'",
            "SELECT * FROM (SELECT * FROM peaks WHERE interval INTERSECTS 'chr1:1-100') t",
            "WITH d AS (SELECT * FROM peaks WHERE interval INTERSECTS 'chr1:1-100') "
            "SELECT * FROM d",
        ],
        ids=["where", "join_on", "subquery", "cte"],
    )
    def test_walk_reaches_operator_in_each_location(self, clean_registry, query):
        """Test that the walk reaches an operator in any clause location.

        Given:
            An INTERSECTS predicate placed in a WHERE, JOIN-ON, subquery, or CTE,
            flagged, with a counting expander.
        When:
            Running the pass.
        Then:
            The operator is dispatched exactly once regardless of location.
        """
        # Arrange
        counting = _CountingExpander()
        clean_registry.register(GenericTarget(), Intersects, counting)
        tables = _tables(("peaks", "genes"))
        ast = _prepare(query, tables)
        pass_ = ExpandOperators(GenericTarget(), tables, clean_registry)

        # Act
        with _opted_in(Intersects):
            result = pass_.transform(ast)

        # Assert
        assert len(counting.calls) == 1
        assert not list(result.find_all(Intersects))

    def test_walk_replacement_serializes_through_generator(self, clean_registry):
        """Test that the pass's replacement serializes cleanly to SQL.

        Given:
            A flagged DISJOIN and an expander returning a column sentinel.
        When:
            Running the pass and serializing the result with BaseGIQLGenerator.
        Then:
            The sentinel appears in the serialized SQL (replacement integrity).
        """
        # Arrange
        clean_registry.register(
            GenericTarget(), GIQLDisjoin, lambda n, c: exp.column("WALK_SENTINEL")
        )
        tables = _tables()
        ast = _prepare("SELECT * FROM DISJOIN(variants)", tables)
        pass_ = ExpandOperators(GenericTarget(), tables, clean_registry)

        # Act
        with _opted_in(GIQLDisjoin):
            result = pass_.transform(ast)
        sql = BaseGIQLGenerator(tables=tables).generate(result)

        # Assert
        assert "WALK_SENTINEL" in sql

    def test_walk_identity_return_leaves_node_in_place(self, clean_registry):
        """Test that an expander returning the node itself triggers no replacement.

        Given:
            A flagged DISJOIN whose expander returns the node unchanged.
        When:
            Running the pass.
        Then:
            The DISJOIN node remains in the tree (the identity-return guard skips
            the replace call).
        """
        # Arrange
        clean_registry.register(GenericTarget(), GIQLDisjoin, lambda n, c: n)
        tables = _tables()
        ast = _prepare("SELECT * FROM DISJOIN(variants)", tables)
        pass_ = ExpandOperators(GenericTarget(), tables, clean_registry)

        # Act
        with _opted_in(GIQLDisjoin):
            result = pass_.transform(ast)

        # Assert
        assert list(result.find_all(GIQLDisjoin))

    def test_walk_dispatches_mixed_operator_types_per_type(self, clean_registry):
        """Test that the walk dispatches each operator type to its own expander.

        Given:
            A query containing both a DISJOIN and an INTERSECTS, both flagged,
            each with its own counting expander.
        When:
            Running the pass.
        Then:
            Each expander is called once, dispatched by operator type.
        """
        # Arrange
        disjoin_counter = _CountingExpander("disjoin")
        intersects_counter = _CountingExpander("intersects")
        clean_registry.register(GenericTarget(), GIQLDisjoin, disjoin_counter)
        clean_registry.register(GenericTarget(), Intersects, intersects_counter)
        tables = _tables(("variants", "peaks"))
        ast = _prepare(
            "SELECT * FROM DISJOIN(variants) "
            "WHERE EXISTS (SELECT * FROM peaks WHERE interval INTERSECTS 'chr1:1-100')",
            tables,
        )
        pass_ = ExpandOperators(GenericTarget(), tables, clean_registry)

        # Act
        with _opted_in(GIQLDisjoin):
            with _opted_in(Intersects):
                pass_.transform(ast)

        # Assert
        assert len(disjoin_counter.calls) == 1
        assert len(intersects_counter.calls) == 1

    def test_walk_partial_opt_in_replaces_only_flagged_type(self, clean_registry):
        """Test that only the flagged operator type is replaced when both registered.

        Given:
            A genuinely-unmigrated operator (GIQLCluster, shipping
            GIQL_EXPAND=False in every wave-3 branch) and an INTERSECTS, both with
            registered expanders, but only INTERSECTS flagged GIQL_EXPAND for the
            test.
        When:
            Running the pass.
        Then:
            The INTERSECTS is replaced while the unmigrated operator node remains
            on its own shipped ``False`` flag — no opt-out ceremony needed (the
            gate is per-type).
        """
        # Arrange — the held-off subject is genuinely unmigrated: it survives on
        # its own shipped GIQL_EXPAND=False, not on a test opt-out.
        held_off = GIQLCluster
        assert held_off in _UNMIGRATED_OPERATORS
        assert held_off.GIQL_EXPAND is False
        clean_registry.register(
            GenericTarget(), held_off, lambda n, c: exp.column("CL")
        )
        clean_registry.register(
            GenericTarget(), Intersects, lambda n, c: exp.column("IX")
        )
        tables = _tables(("peaks",))
        ast = _prepare(
            "SELECT *, CLUSTER(interval) AS cluster_id FROM peaks "
            "WHERE EXISTS (SELECT * FROM peaks WHERE interval INTERSECTS 'chr1:1-100')",
            tables,
        )
        pass_ = ExpandOperators(GenericTarget(), tables, clean_registry)

        # Act — only INTERSECTS is opted in; the unmigrated operator stays off.
        with _opted_in(Intersects):
            result = pass_.transform(ast)

        # Assert
        assert list(result.find_all(held_off))
        assert not list(result.find_all(Intersects))

    def test_walk_shares_alias_sequence_across_sibling_expanders(self, clean_registry):
        """Test that sibling expanders draw from one alias sequence (no collision).

        Given:
            Two sibling DISJOIN operators, flagged, with an expander that mints an
            alias per call and records it.
        When:
            Running the pass.
        Then:
            The two minted aliases differ (the run threads a single
            name_sequence across every per-node context, so no __giql_x_
            collision).
        """
        # Arrange
        seen = []

        def _mint(node, ctx):
            seen.append(ctx.alias())
            return exp.column("c")

        clean_registry.register(GenericTarget(), GIQLDisjoin, _mint)
        tables = _tables(("variants", "genes"))
        ast = _prepare(
            "SELECT * FROM DISJOIN(variants) UNION ALL SELECT * FROM DISJOIN(genes)",
            tables,
        )
        pass_ = ExpandOperators(GenericTarget(), tables, clean_registry)

        # Act
        with _opted_in(GIQLDisjoin):
            pass_.transform(ast)

        # Assert
        assert len(seen) == 2
        assert seen[0] != seen[1]
        assert all(a.startswith(EXPAND_ALIAS_PREFIX) for a in seen)

    def test_walk_expands_inner_before_outer_replaces_subtree(self, clean_registry):
        """Test that a nested operator expands into the live tree before its ancestor.

        Given:
            An outer DISJOIN whose reference subquery contains an inner DISJOIN,
            both flagged. The outer's expander returns a fresh node that does NOT
            re-attach the inner subtree, and inspects its own subtree so the
            recorded order reveals whether the inner had already expanded when the
            outer ran.
        When:
            Running the pass.
        Then:
            BOTH operators expand and the outer sees the inner already gone — the
            deepest-first order expands the inner into the live tree before the
            outer replaces it. Under the former outer-first order
            the outer replacement detached the still-pending inner, whose later
            replace() landed in a discarded subtree, so the inner never reached
            the live tree (#154). The recorded order pins the discriminator: the
            inner must expand first.
        """
        # Arrange
        order = []
        tables = _tables(("variants", "genes"))
        ast = _prepare(
            "SELECT * FROM "
            "DISJOIN(variants, reference := (SELECT * FROM DISJOIN(genes)))",
            tables,
        )
        # Identify the outer node up front: it is the one carrying a descendant
        # DISJOIN. Tag by identity so the label survives the inner's replacement.
        outer = next(
            n for n in ast.find_all(GIQLDisjoin) if list(n.find_all(GIQLDisjoin))[1:]
        )

        def _expander(node, ctx):
            if node is outer:
                # When the outer runs, the inner must already be gone from the
                # outer's subtree (deepest-first expanded and replaced it). Under
                # the old outer-first order an un-expanded inner DISJOIN would
                # still be present here.
                order.append(
                    "outer_sees_no_inner"
                    if not list(node.find_all(GIQLDisjoin))[1:]
                    else "outer_sees_inner"
                )
                return exp.column("EXPANDED_outer")
            order.append("inner")
            return exp.column("EXPANDED_inner")

        clean_registry.register(GenericTarget(), GIQLDisjoin, _expander)
        pass_ = ExpandOperators(GenericTarget(), tables, clean_registry)

        # Act
        with _opted_in(GIQLDisjoin):
            result = pass_.transform(ast)

        # Assert
        # The inner expands first, into the live tree, so the outer sees it gone.
        assert order == ["inner", "outer_sees_no_inner"]
        assert not list(result.find_all(GIQLDisjoin))

    def test_transform_raises_on_node_missing_resolution(self, clean_registry):
        """Test that the pass raises when an operator lacks resolution metadata.

        Given:
            A flagged DISJOIN whose pass-1 resolution metadata has been stripped
            (an internal invariant violation), with a registered expander.
        When:
            Running the pass.
        Then:
            It raises ResolutionError rather than dispatching with a None
            resolution (pass 1 must annotate every operator node).
        """
        # Arrange
        from giql.resolver import META_KEY
        from giql.resolver import ResolutionError

        clean_registry.register(GenericTarget(), GIQLDisjoin, _record("x"))
        tables = _tables()
        ast = _prepare("SELECT * FROM DISJOIN(variants)", tables)
        for node in ast.find_all(GIQLDisjoin):
            node.meta.pop(META_KEY, None)
        pass_ = ExpandOperators(GenericTarget(), tables, clean_registry)

        # Act & assert
        with _opted_in(GIQLDisjoin):
            with pytest.raises(ResolutionError):
                pass_.transform(ast)

    def test_lazy_import_has_no_cycle_in_fresh_process(self):
        """Test that importing the expander module standalone raises no cycle.

        Given:
            A fresh Python process with nothing pre-imported.
        When:
            Importing giql.expander and running its public expand_operators pass
            over a parsed operator AST.
        Then:
            It imports and runs the inert pass without an import cycle, leaving at
            least one operator node intact (the empty registry is a no-op), so the
            lazy operator import resolves through a public entry point.
        """
        # Arrange
        import os
        import subprocess
        import sys

        code = (
            "from giql.expander import expand_operators, ExpanderRegistry; "
            "from giql.expressions import GIQLDisjoin; "
            "from giql.dialect import GIQLDialect; "
            "from giql.table import Table, Tables; "
            "from giql.targets import GenericTarget; "
            "from sqlglot import parse_one; "
            "t = Tables(); t.register('variants', Table('variants')); "
            "ast = parse_one('SELECT * FROM DISJOIN(variants)', dialect=GIQLDialect); "
            "out = expand_operators(ast, GenericTarget(), t, ExpanderRegistry()); "
            "ops = list(out.find_all(GIQLDisjoin)); "
            "assert len(ops) >= 1, ops; "
            "print('ok')"
        )
        # Strip coverage's subprocess auto-start hooks so the child is a genuinely
        # clean interpreter exercising only the import (not the parent run's
        # coverage instrumentation, which is irrelevant to the cycle check).
        env = {
            k: v
            for k, v in os.environ.items()
            if not k.startswith("COV_CORE") and k != "COVERAGE_PROCESS_START"
        }

        # Act
        result = subprocess.run(
            [sys.executable, "-c", code],
            capture_output=True,
            text=True,
            env=env,
        )

        # Assert
        assert result.returncode == 0, result.stderr
        assert result.stdout.strip() == "ok"


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


class _opted_out:
    """Context manager opting an operator class out of GIQL_EXPAND for a test.

    The complement of :class:`_opted_in`: used by a control test that needs a
    *migrated* operator (one shipping GIQL_EXPAND=True) to behave as if
    unflagged, so the test can prove the pass gates per-type without the
    operator's shipped opt-in interfering. Restores the prior flag on exit.
    """

    def __init__(self, operator: type) -> None:
        self._operator = operator
        self._prior = operator.__dict__.get("GIQL_EXPAND", False)

    def __enter__(self):
        self._operator.GIQL_EXPAND = False
        return self._operator

    def __exit__(self, *exc):
        self._operator.GIQL_EXPAND = self._prior
        return False
