"""Engine-free unit tests for the cross-target oracle internals (#139, T1).

These prove the oracle *itself* is sound before any engine runs: that
:func:`normalize` / :func:`scalar` collapse engine quirks without hiding real
differences, that :func:`assert_cross_target` genuinely raises on a wrong
``expected`` and on engine disagreement, and that :func:`resolve_routing` maps
and guards targets/engines. No DuckDB / DataFusion / pyarrow needed, so they
run on every machine and localize oracle bugs away from engine bugs.
"""

import math

import pytest

from tests.integration._oracle import assert_cross_target
from tests.integration._oracle import normalize
from tests.integration._oracle import resolve_routing
from tests.integration._oracle import scalar


class _FakeNumpyScalar:
    """A scalar that only reveals its value (and NaN-ness) after ``.item()``."""

    def __init__(self, value):
        self._value = value

    def item(self):
        return self._value


class TestScalar:
    """`scalar` coercion of engine-specific values."""

    def test_scalar_none_returns_none(self):
        """Test None passes through as None.

        Given:
            A None value.
        When:
            scalar() coerces it.
        Then:
            It should return None.
        """
        # Arrange / Act / Assert
        assert scalar(None) is None

    def test_scalar_plain_float_nan_returns_none(self):
        """Test a plain float NaN collapses to None.

        Given:
            A bare float NaN as a DuckDB/DataFusion null surrogate.
        When:
            scalar() coerces it.
        Then:
            It should return None so it compares equal to a real NULL.
        """
        # Arrange / Act / Assert
        assert scalar(float("nan")) is None

    def test_scalar_numpy_nan_after_item_returns_none(self):
        """Test a NaN that only surfaces after .item() collapses to None.

        Given:
            A numpy-like scalar whose .item() unwraps to a float NaN.
        When:
            scalar() coerces it and re-checks NaN AFTER .item() (Finding 6).
        Then:
            It should return None, never a bare NaN sort key.
        """
        # Arrange
        value = _FakeNumpyScalar(float("nan"))

        # Act
        result = scalar(value)

        # Assert
        assert result is None

    def test_scalar_numpy_int_returns_plain_int(self):
        """Test a numpy-like int scalar unwraps to a plain Python int.

        Given:
            A numpy-like scalar wrapping an int.
        When:
            scalar() coerces it.
        Then:
            It should return the plain int value.
        """
        # Arrange / Act
        result = scalar(_FakeNumpyScalar(42))

        # Assert
        assert result == 42
        assert isinstance(result, int)


class TestNormalize:
    """`normalize` multiset construction."""

    def test_normalize_sorts_rows_order_free(self):
        """Test normalize returns a sorted, order-free row list.

        Given:
            Two rows supplied in descending order.
        When:
            normalize() builds the multiset.
        Then:
            The result should be sorted ascending regardless of input order.
        """
        # Arrange / Act
        result = normalize([("chr2", 5, 6), ("chr1", 1, 2)])

        # Assert
        assert result == [("chr1", 1, 2), ("chr2", 5, 6)]

    def test_normalize_preserves_multiplicity(self):
        """Test normalize keeps duplicate rows (list, not set).

        Given:
            Two identical rows.
        When:
            normalize() builds the multiset.
        Then:
            Both copies should survive so duplicate-row bugs surface.
        """
        # Arrange / Act
        result = normalize([("chr1", 1, 2), ("chr1", 1, 2)])

        # Assert
        assert result == [("chr1", 1, 2), ("chr1", 1, 2)]

    def test_normalize_is_determinor_under_shuffle(self):
        """Test normalize yields an identical result under row shuffling.

        Given:
            A row multiset and a reversed copy of it.
        When:
            normalize() is applied to both orderings.
        Then:
            The two normalized results should be identical (determinism
            self-test), including for None-bearing rows.
        """
        # Arrange
        rows = [
            ("chr1", 1, 2),
            ("chr1", None, 9),
            ("chr2", 3, 4),
            ("chr1", 1, 2),
        ]

        # Act
        forward = normalize(rows)
        reverse = normalize(list(reversed(rows)))

        # Assert
        assert forward == reverse

    def test_normalize_treats_nan_and_none_alike(self):
        """Test a NaN and a None sort into the same normalized position.

        Given:
            One row with a None and one with a float NaN in the same column.
        When:
            normalize() coerces both.
        Then:
            They should compare equal so engine null surrogates do not diverge.
        """
        # Arrange / Act
        with_none = normalize([("chr1", None, 2)])
        with_nan = normalize([("chr1", math.nan, 2)])

        # Assert
        assert with_none == with_nan


class TestAssertCrossTarget:
    """`assert_cross_target` differential + expectation phases."""

    def test_passes_when_targets_agree_and_match_expected(self):
        """Test no error when every target agrees and equals expected.

        Given:
            Two targets returning the same normalized rows and a matching
            expected multiset.
        When:
            assert_cross_target() runs.
        Then:
            It should not raise.
        """
        # Arrange
        rows = normalize([("chr1", 1, 2)])
        results = {"generic": rows, "duckdb": rows}

        # Act / Assert
        assert_cross_target(results, normalize([("chr1", 1, 2)]))

    def test_raises_on_value_difference_from_expected(self):
        """Test a value diff from expected raises in the expectation phase.

        Given:
            Agreeing targets whose value differs from expected.
        When:
            assert_cross_target() runs.
        Then:
            It should raise AssertionError naming the expectation mismatch.
        """
        # Arrange
        rows = normalize([("chr1", 1, 2)])
        results = {"generic": rows, "duckdb": rows}

        # Act / Assert
        with pytest.raises(AssertionError):
            assert_cross_target(results, normalize([("chr1", 1, 99)]))

    def test_raises_on_extra_row_versus_expected(self):
        """Test an extra row over expected raises.

        Given:
            Agreeing targets returning one more row than expected.
        When:
            assert_cross_target() runs.
        Then:
            It should raise AssertionError.
        """
        # Arrange
        rows = normalize([("chr1", 1, 2), ("chr1", 3, 4)])
        results = {"generic": rows, "duckdb": rows}

        # Act / Assert
        with pytest.raises(AssertionError):
            assert_cross_target(results, normalize([("chr1", 1, 2)]))

    def test_raises_on_missing_row_versus_expected(self):
        """Test a missing row versus expected raises.

        Given:
            Agreeing targets returning one fewer row than expected.
        When:
            assert_cross_target() runs.
        Then:
            It should raise AssertionError.
        """
        # Arrange
        rows = normalize([("chr1", 1, 2)])
        results = {"generic": rows, "duckdb": rows}

        # Act / Assert
        with pytest.raises(AssertionError):
            assert_cross_target(results, normalize([("chr1", 1, 2), ("chr1", 3, 4)]))

    def test_raises_on_empty_versus_nonempty_expected(self):
        """Test an empty result against a non-empty expectation raises.

        Given:
            Agreeing targets returning zero rows but a non-empty expected set.
        When:
            assert_cross_target() runs.
        Then:
            It should raise AssertionError.
        """
        # Arrange
        results = {"generic": [], "duckdb": []}

        # Act / Assert
        with pytest.raises(AssertionError):
            assert_cross_target(results, normalize([("chr1", 1, 2)]))

    def test_raises_on_nonempty_versus_empty_expected(self):
        """Test a non-empty result against an empty expectation raises.

        Given:
            Agreeing targets returning a row but an empty expected set.
        When:
            assert_cross_target() runs.
        Then:
            It should raise AssertionError.
        """
        # Arrange
        rows = normalize([("chr1", 1, 2)])
        results = {"generic": rows, "duckdb": rows}

        # Act / Assert
        with pytest.raises(AssertionError):
            assert_cross_target(results, [])

    def test_raises_on_multiset_multiplicity_divergence(self):
        """Test differing duplicate counts raise (list-not-set semantics).

        Given:
            Agreeing targets returning a row twice but expected only once.
        When:
            assert_cross_target() runs.
        Then:
            It should raise because multiplicity is preserved, not collapsed.
        """
        # Arrange
        rows = normalize([("chr1", 1, 2), ("chr1", 1, 2)])
        results = {"generic": rows, "duckdb": rows}

        # Act / Assert
        with pytest.raises(AssertionError):
            assert_cross_target(results, normalize([("chr1", 1, 2)]))

    def test_raises_on_none_versus_value_difference(self):
        """Test a None where expected holds a value raises.

        Given:
            Agreeing targets with a None column but expected with a real value.
        When:
            assert_cross_target() runs.
        Then:
            It should raise: None and a value are genuinely distinct.
        """
        # Arrange
        rows = normalize([("chr1", None, 2)])
        results = {"generic": rows, "duckdb": rows}

        # Act / Assert
        with pytest.raises(AssertionError):
            assert_cross_target(results, normalize([("chr1", 5, 2)]))

    def test_coercion_normalizes_type_but_preserves_value_diff(self):
        """Test type coercion does not mask a genuine value difference.

        Given:
            Targets whose values are numpy-like scalars equal to expected ints
            in one case and differing in another.
        When:
            assert_cross_target() compares the coerced multisets.
        Then:
            The equal case passes and the differing case raises — coercion
            normalizes type without hiding a real difference.
        """
        # Arrange
        equal_rows = normalize([("chr1", _FakeNumpyScalar(1), _FakeNumpyScalar(2))])
        diff_rows = normalize([("chr1", _FakeNumpyScalar(1), _FakeNumpyScalar(2))])

        # Act / Assert
        assert_cross_target(
            {"generic": equal_rows, "duckdb": equal_rows}, normalize([("chr1", 1, 2)])
        )
        with pytest.raises(AssertionError):
            assert_cross_target(
                {"generic": diff_rows, "duckdb": diff_rows},
                normalize([("chr1", 1, 99)]),
            )

    def test_raises_with_engines_disagree_message(self):
        """Test two disagreeing engines trigger the differential assertion.

        Given:
            A results dict where the duckdb target returns a different row from
            the generic target, both matching neither each other.
        When:
            assert_cross_target() runs (now reachable after Finding 1).
        Then:
            It should raise AssertionError whose message says the engines
            disagree and names both targets — independent of expected.
        """
        # Arrange
        generic_rows = normalize([("chr1", 1, 2)])
        duckdb_rows = normalize([("chr1", 1, 3)])
        results = {"generic": generic_rows, "duckdb": duckdb_rows}

        # Act
        with pytest.raises(AssertionError) as excinfo:
            assert_cross_target(results, generic_rows)

        # Assert
        message = str(excinfo.value)
        assert "engines disagree" in message
        assert "generic" in message and "duckdb" in message

    def test_disagreement_checked_before_expectation(self):
        """Test the differential phase fires even when one target matches expected.

        Given:
            A generic target equal to expected and a duckdb target that differs,
            so a mid-loop ``== expected`` design would have passed generic and
            never compared the engines.
        When:
            assert_cross_target() runs.
        Then:
            It should raise the engines-disagree error first, proving Finding 1
            made the cross-target assertion genuinely reachable.
        """
        # Arrange
        expected = normalize([("chr1", 1, 2)])
        results = {"generic": expected, "duckdb": normalize([("chr1", 9, 9)])}

        # Act / Assert
        with pytest.raises(AssertionError, match="engines disagree"):
            assert_cross_target(results, expected)


class TestResolveRouting:
    """`resolve_routing` default map, overrides, and guards."""

    def test_default_routing_maps_each_target_to_its_engine(self):
        """Test the default routing maps targets to their intended engines.

        Given:
            All three default targets and no overrides.
        When:
            resolve_routing() resolves them.
        Then:
            generic/datafusion should route to DataFusion and duckdb to DuckDB,
            each carrying its transpile dialect.
        """
        # Arrange / Act
        routing = resolve_routing(("generic", "datafusion", "duckdb"))

        # Assert
        assert routing == {
            "generic": ("datafusion", None),
            "datafusion": ("datafusion", "datafusion"),
            "duckdb": ("duckdb", "duckdb"),
        }

    def test_engines_override_reroutes_target(self):
        """Test an engines override reroutes a target's engine but not its dialect.

        Given:
            The generic and duckdb targets with generic routed to duckdb.
        When:
            resolve_routing() applies the override.
        Then:
            generic should execute on duckdb while keeping its None dialect.
        """
        # Arrange / Act
        routing = resolve_routing(("generic", "duckdb"), engines={"generic": "duckdb"})

        # Assert
        assert routing["generic"] == ("duckdb", None)
        assert routing["duckdb"] == ("duckdb", "duckdb")

    def test_unknown_target_raises_value_error(self):
        """Test an unknown target name raises a clear ValueError (Finding 7).

        Given:
            A target name that is not a registered target.
        When:
            resolve_routing() resolves it.
        Then:
            It should raise ValueError, not a bare KeyError.
        """
        # Arrange / Act / Assert
        with pytest.raises(ValueError, match="Unknown target"):
            resolve_routing(("postgres",))

    def test_unknown_engine_override_raises_value_error(self):
        """Test an override to an unknown engine raises a clear ValueError.

        Given:
            A valid target rerouted to an engine the oracle cannot run.
        When:
            resolve_routing() resolves it.
        Then:
            It should raise ValueError naming the unknown engine.
        """
        # Arrange / Act / Assert
        with pytest.raises(ValueError, match="Unknown engine"):
            resolve_routing(("generic",), engines={"generic": "spark"})

    def test_unknown_target_in_override_raises_value_error(self):
        """Test an override keyed on an unknown target raises a clear ValueError.

        Given:
            An engines override naming a target that does not exist.
        When:
            resolve_routing() resolves it.
        Then:
            It should raise ValueError naming the unknown override target.
        """
        # Arrange / Act / Assert
        with pytest.raises(ValueError, match="Unknown target"):
            resolve_routing(("generic",), engines={"nope": "duckdb"})
