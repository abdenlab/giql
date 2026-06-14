"""Tests for the target-engine model."""

import re
from dataclasses import FrozenInstanceError

import pytest
from hypothesis import given
from hypothesis import strategies as st

from giql.targets import Capabilities
from giql.targets import DataFusionTarget
from giql.targets import DuckDBTarget
from giql.targets import GenericTarget
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
            range_join_strategy="binned",
        )

        # Assert
        assert caps.supports_lateral is True
        assert caps.supports_star_replace is False
        assert caps.supports_qualify is True
        assert caps.range_join_strategy == "binned"

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
            range_join_strategy="binned",
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
            range_join_strategy="binned",
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
            lateral supported, no star-REPLACE, no QUALIFY, binned joins.
        """
        # Act
        target = GenericTarget()

        # Assert
        assert target.name == "generic"
        assert target.sqlglot_dialect is None
        assert target.capabilities.supports_lateral is True
        assert target.capabilities.supports_star_replace is False
        assert target.capabilities.supports_qualify is False
        assert target.capabilities.range_join_strategy == "binned"


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
            the binned join strategy.
        """
        # Act
        target = DataFusionTarget()

        # Assert
        assert target.name == "datafusion"
        assert target.sqlglot_dialect is None
        assert target.capabilities.supports_lateral is False
        assert target.capabilities.supports_star_replace is False
        assert target.capabilities.supports_qualify is False
        assert target.capabilities.range_join_strategy == "binned"


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
        + re.escape("'duckdb', 'datafusion', or None.")
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
