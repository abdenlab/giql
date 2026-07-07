"""Tests for the target-engine model."""

import re
from dataclasses import FrozenInstanceError
from dataclasses import dataclass

import pytest
from hypothesis import given
from hypothesis import strategies as st
from sqlglot import parse_one

from giql import REGISTRY
from giql import ExpansionContext
from giql import register
from giql import transpile
from giql.expressions import Within
from giql.targets import Capabilities
from giql.targets import DataFusionTarget
from giql.targets import DuckDBTarget
from giql.targets import GenericTarget
from giql.targets import Target
from giql.targets import resolve_target


class TestCapabilities:
    """Tests for the Capabilities dataclass."""

    def test___init___stores_all_descriptors(self):
        """Test that Capabilities records each descriptor.

        Given:
            A full set of capability descriptor values.
        When:
            A Capabilities instance is constructed.
        Then:
            It should expose each value on the corresponding attribute.
        """
        # Act
        caps = Capabilities(
            supports_lateral=True,
            supports_star_replace=False,
            supports_qualify=True,
            range_join_strategy="naive",
        )

        # Assert
        assert caps.supports_lateral is True
        assert caps.supports_star_replace is False
        assert caps.supports_qualify is True
        assert caps.range_join_strategy == "naive"

    def test___eq___with_identical_values(self):
        """Test value equality of Capabilities.

        Given:
            Two Capabilities built from identical descriptor values.
        When:
            They are compared for equality.
        Then:
            It should treat them as equal (frozen value semantics).
        """
        # Arrange
        first = Capabilities(
            supports_lateral=True,
            supports_star_replace=True,
            supports_qualify=True,
            range_join_strategy="iejoin",
        )
        second = Capabilities(
            supports_lateral=True,
            supports_star_replace=True,
            supports_qualify=True,
            range_join_strategy="iejoin",
        )

        # Act & assert
        assert first == second

    def test___eq___with_one_differing_field(self):
        """Test value inequality of Capabilities.

        Given:
            Two Capabilities identical in every descriptor except one.
        When:
            They are compared for equality.
        Then:
            It should treat them as unequal (per-field value semantics).
        """
        # Arrange
        base = Capabilities(
            supports_lateral=True,
            supports_star_replace=True,
            supports_qualify=True,
            range_join_strategy="iejoin",
        )
        differing = Capabilities(
            supports_lateral=True,
            supports_star_replace=True,
            supports_qualify=True,
            range_join_strategy="naive",
        )

        # Act & assert
        assert base != differing

    def test___setattr___is_frozen(self):
        """Test that Capabilities is immutable.

        Given:
            A constructed Capabilities instance.
        When:
            An attribute is reassigned.
        Then:
            It should raise FrozenInstanceError.
        """
        # Arrange
        caps = Capabilities(
            supports_lateral=True,
            supports_star_replace=False,
            supports_qualify=False,
            range_join_strategy="naive",
        )

        # Act & assert
        with pytest.raises(FrozenInstanceError):
            caps.supports_lateral = False  # type: ignore[misc]


class TestGenericTarget:
    """Tests for the GenericTarget."""

    def test___init___exposes_portable_baseline(self):
        """Test the generic target's identity and capabilities.

        Given:
            The GenericTarget class.
        When:
            An instance is constructed.
        Then:
            It should carry the portable SQL-92 baseline: no engine dialect,
            lateral supported, no star-REPLACE, no QUALIFY, naive-predicate joins.
        """
        # Act
        target = GenericTarget()

        # Assert
        assert target.name == "generic"
        assert target.sqlglot_dialect is None
        assert target.capabilities.supports_lateral is True
        assert target.capabilities.supports_star_replace is False
        assert target.capabilities.supports_qualify is False
        assert target.capabilities.range_join_strategy == "naive"


class TestDuckDBTarget:
    """Tests for the DuckDBTarget."""

    def test___init___exposes_duckdb_features(self):
        """Test the DuckDB target's identity and capabilities.

        Given:
            The DuckDBTarget class.
        When:
            An instance is constructed.
        Then:
            It should serialize through the duckdb dialect and enable
            lateral, star-REPLACE, QUALIFY, and the IEJoin strategy.
        """
        # Act
        target = DuckDBTarget()

        # Assert
        assert target.name == "duckdb"
        assert target.sqlglot_dialect == "duckdb"
        assert target.capabilities.supports_lateral is True
        assert target.capabilities.supports_star_replace is True
        assert target.capabilities.supports_qualify is True
        assert target.capabilities.range_join_strategy == "iejoin"


class TestDataFusionTarget:
    """Tests for the DataFusionTarget."""

    def test___init___exposes_conservative_capabilities(self):
        """Test the DataFusion target's identity and capabilities.

        Given:
            The DataFusionTarget class.
        When:
            An instance is constructed.
        Then:
            It should fall back to generic serialization (no sqlglot
            DataFusion dialect), disable star-REPLACE and QUALIFY, and use
            the naive-predicate join strategy.
        """
        # Act
        target = DataFusionTarget()

        # Assert
        assert target.name == "datafusion"
        assert target.sqlglot_dialect is None
        assert target.capabilities.supports_lateral is False
        assert target.capabilities.supports_star_replace is False
        assert target.capabilities.supports_qualify is False
        assert target.capabilities.range_join_strategy == "naive"


class TestTarget:
    """Tests for the shared Target value semantics."""

    def test___eq___same_target_type(self):
        """Test that two instances of the same target are value-equal.

        Given:
            Two separately constructed DuckDBTarget instances.
        When:
            They are compared for equality.
        Then:
            It should treat them as equal (frozen value semantics), so a
            resolved target is a stable registry key.
        """
        # Act & assert
        assert DuckDBTarget() == DuckDBTarget()

    def test___eq___different_target_types(self):
        """Test that distinct target classes are never equal.

        Given:
            A GenericTarget and a DataFusionTarget (whose fields partly
            overlap).
        When:
            They are compared for equality.
        Then:
            It should treat them as unequal, since equality is class-scoped.
        """
        # Act & assert
        assert GenericTarget() != DataFusionTarget()

    def test___hash___supports_set_and_dict_keys(self):
        """Test that targets are hashable and usable as keys.

        Given:
            Value-equal and value-distinct target instances.
        When:
            They are placed in a set and used as dict keys.
        Then:
            It should dedup equal instances and resolve a key built from a
            separately constructed equal instance — the contract the
            operator-expander registry (#138) relies on.
        """
        # Arrange
        registry = {(DuckDBTarget(), "Intersects"): "expander"}

        # Act & assert
        assert len({DuckDBTarget(), DuckDBTarget(), GenericTarget()}) == 2
        assert registry[(DuckDBTarget(), "Intersects")] == "expander"


def test_resolve_target_with_none_returns_generic():
    """Test resolution of the default dialect.

    Given:
        A None dialect argument.
    When:
        resolve_target is called.
    Then:
        It should return a GenericTarget.
    """
    # Act
    target = resolve_target(None)

    # Assert
    assert isinstance(target, GenericTarget)


def test_resolve_target_with_duckdb_returns_duckdb():
    """Test resolution of the duckdb dialect.

    Given:
        The dialect name "duckdb".
    When:
        resolve_target is called.
    Then:
        It should return a DuckDBTarget.
    """
    # Act
    target = resolve_target("duckdb")

    # Assert
    assert isinstance(target, DuckDBTarget)


def test_resolve_target_with_datafusion_returns_datafusion():
    """Test resolution of the datafusion dialect.

    Given:
        The dialect name "datafusion".
    When:
        resolve_target is called.
    Then:
        It should return a DataFusionTarget.
    """
    # Act
    target = resolve_target("datafusion")

    # Assert
    assert isinstance(target, DataFusionTarget)


@pytest.mark.parametrize(
    "dialect",
    [
        "postgres",  # plain unknown engine
        "sqlite",  # second unknown — guards against a hardcoded message
        "DuckDB",  # case-variant of a real target — lookup is case-sensitive
        "",  # empty string is distinct from None
        "generic",  # the internal target name — None is the only public path
    ],
)
def test_resolve_target_with_unsupported_dialect_raises(dialect):
    """Test resolution of dialect names that are not public targets.

    Given:
        A dialect string that is not one of the supported public names
        (an unknown engine, a case-variant, the empty string, or the
        internal "generic" name).
    When:
        resolve_target is called.
    Then:
        It should raise ValueError naming the offending value and the
        supported targets.
    """
    # Act & assert
    pattern = (
        re.escape(f"Unknown dialect: {dialect!r}.")
        + r".*"
        + re.escape("'duckdb', 'datafusion', None,")
    )
    with pytest.raises(ValueError, match=pattern):
        resolve_target(dialect)


@given(st.text().filter(lambda s: s not in ("duckdb", "datafusion")))
def test_resolve_target_with_arbitrary_unsupported_string_raises(dialect):
    """Test resolution over the open domain of unsupported strings.

    Given:
        Any string that is not a registered public dialect name (the only
        public string names are "duckdb" and "datafusion").
    When:
        resolve_target is called.
    Then:
        It should raise ValueError rather than returning a target.
    """
    # Act & assert
    with pytest.raises(ValueError, match="Unknown dialect"):
        resolve_target(dialect)


@dataclass(frozen=True)
class _PostgresTarget(Target):
    """A capability-only custom target used to exercise the plugin hub."""

    name: str = "postgres"
    sqlglot_dialect: str | None = "postgres"
    capabilities: Capabilities = Capabilities(
        supports_lateral=True,
        supports_star_replace=False,
        supports_qualify=True,
        range_join_strategy="naive",
    )


#: The process-wide REGISTRY state at import (built-in expanders). The custom
#: target tests below mutate the global registry, so — mirroring test_expander.py
#: — an autouse guard asserts the baseline is restored at every test boundary,
#: catching any leaked target/expander (e.g. the "postgres" name this file both
#: registers and, elsewhere, expects to be unresolved) deterministically.
_REGISTRY_BASELINE = REGISTRY.snapshot()


@pytest.fixture(autouse=True)
def _registry_leak_guard():
    assert REGISTRY.snapshot() == _REGISTRY_BASELINE, (
        "REGISTRY differed from its baseline entering a test"
    )
    yield
    assert REGISTRY.snapshot() == _REGISTRY_BASELINE, (
        "a test leaked a registration into REGISTRY"
    )


@pytest.fixture
def clean_registry():
    saved = REGISTRY.snapshot()
    try:
        yield
    finally:
        REGISTRY.restore(saved)


class TestCustomTargetInjection:
    """Registering a custom target makes it selectable via ``transpile``."""

    def test_register_target_makes_custom_target_resolvable(self, clean_registry):
        """Test that a declared custom target resolves by name.

        Given:
            A custom Target declared through REGISTRY.register_target.
        When:
            resolve_target is called with its name.
        Then:
            It should return that exact target instance.
        """
        # Arrange
        target = _PostgresTarget()
        REGISTRY.register_target(target)

        # Act & assert
        assert resolve_target("postgres") is target

    def test_transpile_resolves_capability_only_custom_target(self, clean_registry):
        """Test end-to-end transpilation against a capability-only custom target.

        Given:
            A custom target that overrides no operators, only its capabilities.
        When:
            A WITHIN query is transpiled with dialect set to its name.
        Then:
            It should transpile through the generic expanders (its
            ``supports_star_replace=False`` selects the portable form) rather
            than raising ``Unknown dialect``.
        """
        # Arrange
        REGISTRY.register_target(_PostgresTarget())
        query = "SELECT * FROM peaks WHERE interval WITHIN 'chr1:1000-5000'"

        # Act
        output = transpile(query, tables=["peaks"], dialect="postgres")

        # Assert
        assert output == transpile(query, tables=["peaks"])  # generic form

    def test_register_auto_declares_its_target(self, clean_registry):
        """Test that registering an expander also declares its target by name.

        Given:
            An operator expander registered against a custom target via
            @register, with no separate register_target call.
        When:
            resolve_target is called with the target's name.
        Then:
            It should resolve the custom target (declared as a side effect).
        """

        # Arrange
        @register(_PostgresTarget, Within)
        def _noop(node, ctx: ExpansionContext):
            return node

        # Act & assert
        assert resolve_target("postgres") == _PostgresTarget()

    def test_operator_override_on_custom_target_changes_output(self, clean_registry):
        """Test that a per-target expander override reshapes the emitted SQL.

        Given:
            A custom target with a WITHIN expander registered that emits a
            BETWEEN predicate instead of the generic two-sided comparison.
        When:
            A WITHIN query is transpiled against that target.
        Then:
            The emitted SQL should carry the override's BETWEEN form.
        """

        # Arrange
        @register(_PostgresTarget, Within)
        def _within_between(node, ctx: ExpansionContext):
            col = ctx.resolution.column("this")
            return parse_one(f"{col.start} BETWEEN 1000 AND 5000")

        query = "SELECT * FROM peaks WHERE interval WITHIN 'chr1:1000-5000'"

        # Act
        output = transpile(query, tables=["peaks"], dialect="postgres")

        # Assert
        assert "BETWEEN 1000 AND 5000" in output

    def test_register_target_skips_generic(self, clean_registry):
        """Test that GenericTarget is never selectable by the name "generic".

        Given:
            An attempt to declare GenericTarget on the registry.
        When:
            resolve_target is called with "generic".
        Then:
            It should still raise (None is the sole way to select the generic
            target), because register_target is a no-op for GenericTarget.
        """
        # Arrange
        REGISTRY.register_target(GenericTarget())

        # Act & assert
        assert REGISTRY.target("generic") is None
        with pytest.raises(ValueError, match="Unknown dialect"):
            resolve_target("generic")

    def test_resolve_target_builtin_name_shadows_registry(self, clean_registry):
        """Test that a built-in dialect name resolves the built-in, not a shadow.

        Given:
            A custom Target registered under a built-in name ("duckdb").
        When:
            resolve_target is called with that name.
        Then:
            It should return the built-in DuckDBTarget — ``_TARGETS_BY_NAME`` is
            consulted before the plugin registry, so a custom target cannot
            hijack a reserved name.
        """

        # Arrange
        @dataclass(frozen=True)
        class _ShadowTarget(Target):
            name: str = "duckdb"
            sqlglot_dialect: str | None = "postgres"
            capabilities: Capabilities = Capabilities(
                supports_lateral=True,
                supports_star_replace=False,
                supports_qualify=False,
                range_join_strategy="naive",
            )

        REGISTRY.register_target(_ShadowTarget())

        # Act
        resolved = resolve_target("duckdb")

        # Assert
        assert type(resolved) is DuckDBTarget

    def test_resolve_target_error_names_custom_registration_path(self):
        """Test that the unknown-dialect error points at the plugin registry.

        Given:
            A dialect name that is neither built-in nor registered.
        When:
            resolve_target is called with it.
        Then:
            The ValueError message should mention register_target, guiding the
            user toward the custom-target extension path.
        """
        # Act & assert
        with pytest.raises(ValueError, match="register_target"):
            resolve_target("no-such-engine")


class TestTargetDrivesSerialization:
    """The active target's ``sqlglot_dialect`` selects the stock serializer."""

    def test_duckdb_dialect_emits_engine_specific_ordering(self):
        """Test that the duckdb target serializes through the duckdb dialect.

        Given:
            A DISJOIN query whose expansion carries a window ORDER BY.
        When:
            It is transpiled for the duckdb target versus the generic target.
        Then:
            The duckdb output should make DuckDB's default null ordering
            explicit (``NULLS FIRST``) while the generic (dialect-less) output
            should not — confirming ``sqlglot_dialect`` drives serialization.
        """
        # Arrange
        query = "SELECT * FROM DISJOIN(peaks)"

        # Act
        generic_sql = transpile(query, tables=["peaks"])
        duckdb_sql = transpile(query, tables=["peaks"], dialect="duckdb")

        # Assert
        assert "NULLS FIRST" in duckdb_sql
        assert "NULLS FIRST" not in generic_sql

    def test_custom_target_sqlglot_dialect_drives_serialization(self, clean_registry):
        """Test that a registered custom target's sqlglot_dialect reaches ast.sql.

        Given:
            A capability-only custom target whose ``sqlglot_dialect="duckdb"``.
        When:
            A window-carrying operator (DISJOIN) is transpiled under its name
            versus the generic (dialect-less) target.
        Then:
            The custom output should carry the duckdb-dialect ``NULLS FIRST``
            that the generic output lacks — proving a *registered* target's
            ``sqlglot_dialect`` is threaded into serialization, not only a
            built-in's.
        """

        # Arrange
        @dataclass(frozen=True)
        class _DuckLikeTarget(Target):
            name: str = "ducklike"
            sqlglot_dialect: str | None = "duckdb"
            capabilities: Capabilities = Capabilities(
                supports_lateral=True,
                supports_star_replace=False,
                supports_qualify=False,
                range_join_strategy="naive",
            )

        REGISTRY.register_target(_DuckLikeTarget())
        query = "SELECT * FROM DISJOIN(peaks)"

        # Act
        generic_sql = transpile(query, tables=["peaks"])
        custom_sql = transpile(query, tables=["peaks"], dialect="ducklike")

        # Assert
        assert "NULLS FIRST" in custom_sql
        assert "NULLS FIRST" not in generic_sql
