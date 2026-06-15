"""Capability-truthing and lane/CI guards for the cross-target oracle (#139, T5).

These pin the invariants the oracle's routing and xfails depend on: the
DataFusion target's declared capabilities (no LATERAL, no ``* REPLACE``), that
the engine-free oracle internals are importable and reusable from a non-
DataFusion lane, the skip-path behaviour when an engine is absent, and a
version-floor guard for the optional engines so a too-old DuckDB / DataFusion
fails with a clear message rather than a cryptic runtime error.
"""

import importlib.util

import pytest

from giql.targets import DataFusionTarget
from giql.targets import DuckDBTarget
from tests.integration._oracle import assert_cross_target
from tests.integration._oracle import normalize
from tests.integration._oracle import resolve_routing


class TestDataFusionCapabilityTruthing:
    """The DataFusion capability flags the oracle's xfails rely on."""

    def test_datafusion_does_not_support_lateral(self):
        """Test DataFusionTarget declares no LATERAL support.

        Given:
            A DataFusionTarget instance.
        When:
            Its capabilities are inspected.
        Then:
            supports_lateral should be False, matching the NEAREST xfail (#142).
        """
        # Arrange / Act / Assert
        assert DataFusionTarget().capabilities.supports_lateral is False

    def test_datafusion_does_not_support_star_replace(self):
        """Test DataFusionTarget declares no ``* REPLACE`` support.

        Given:
            A DataFusionTarget instance.
        When:
            Its capabilities are inspected.
        Then:
            supports_star_replace should be False (DataFusion has EXCEPT/EXCLUDE
            only).
        """
        # Arrange / Act / Assert
        assert DataFusionTarget().capabilities.supports_star_replace is False

    def test_duckdb_supports_star_replace(self):
        """Test DuckDBTarget declares ``* REPLACE`` support as a contrast.

        Given:
            A DuckDBTarget instance.
        When:
            Its capabilities are inspected.
        Then:
            supports_star_replace should be True, confirming the asymmetry the
            DISJOIN note describes is about duplicate output names, not a
            ``* REPLACE`` capability gap.
        """
        # Arrange / Act / Assert
        assert DuckDBTarget().capabilities.supports_star_replace is True


class TestSharedOracleReuse:
    """The engine-free internals are importable from any lane."""

    def test_oracle_internals_reusable_without_engines(self):
        """Test the comparison core works imported from a non-DataFusion lane.

        Given:
            Two agreeing normalized results assembled with no engine present.
        When:
            assert_cross_target compares them against a matching expectation.
        Then:
            It should pass, proving the shared module is reusable outside the
            DataFusion lane (e.g. for the bedtools / coordinate_space lanes).
        """
        # Arrange
        rows = normalize([("chr1", 1, 2)])

        # Act / Assert
        assert_cross_target({"a": rows, "b": rows}, normalize([("chr1", 1, 2)]))

    def test_routing_resolves_without_engines(self):
        """Test resolve_routing is a pure function needing no engine import.

        Given:
            The default targets.
        When:
            resolve_routing resolves them with no engine installed or imported.
        Then:
            It should return the full routing map, confirming routing is a pure
            decision (Finding 7).
        """
        # Arrange / Act
        routing = resolve_routing(("generic", "datafusion", "duckdb"))

        # Assert
        assert set(routing) == {"generic", "datafusion", "duckdb"}


class TestEngineSkipPaths:
    """Skip / availability behaviour for the optional engines."""

    def test_duckdb_available_or_skipped(self):
        """Test the DuckDB skip-path: present engine imports, absent one skips.

        Given:
            The optional DuckDB dependency, which may or may not be installed.
        When:
            The lane probes for it via importorskip.
        Then:
            It should import when present (and otherwise skip cleanly), mirroring
            the lane's module-level guard.
        """
        # Arrange / Act
        duckdb = pytest.importorskip("duckdb")

        # Assert
        assert hasattr(duckdb, "connect")

    def test_datafusion_available_or_skipped(self):
        """Test the DataFusion skip-path: present engine imports, absent skips.

        Given:
            The optional DataFusion and pyarrow dependencies.
        When:
            The lane probes for them via importorskip.
        Then:
            They should import when present (and otherwise skip cleanly).
        """
        # Arrange / Act
        datafusion = pytest.importorskip("datafusion")
        pytest.importorskip("pyarrow")

        # Assert
        assert hasattr(datafusion, "SessionContext")


class TestEngineVersionFloor:
    """Version-floor guards for the optional engines."""

    def test_duckdb_meets_version_floor(self):
        """Test the installed DuckDB meets the lane's minimum version.

        Given:
            The optional DuckDB dependency.
        When:
            Its version tuple is compared to the floor the IEJoin SQL needs.
        Then:
            It should be at least the supported floor, failing loudly otherwise.
        """
        # Arrange
        duckdb = pytest.importorskip("duckdb")
        floor = (0, 9)

        # Act
        version = tuple(int(p) for p in duckdb.__version__.split(".")[:2])

        # Assert
        assert version >= floor, f"DuckDB {duckdb.__version__} below floor {floor}"

    def test_datafusion_meets_version_floor(self):
        """Test the installed DataFusion meets the lane's minimum version.

        Given:
            The optional DataFusion dependency.
        When:
            Its version tuple is compared to the floor the lane was validated on.
        Then:
            It should be at least the supported floor, failing loudly otherwise.
        """
        # Arrange
        datafusion = pytest.importorskip("datafusion")
        floor = (40, 0)

        # Act
        version = tuple(int(p) for p in datafusion.__version__.split(".")[:2])

        # Assert
        assert version >= floor, (
            f"DataFusion {datafusion.__version__} below floor {floor}"
        )


class TestLaneCollection:
    """A smoke check that the lane and oracle modules collect."""

    def test_oracle_modules_are_importable(self):
        """Test the oracle's importable modules resolve by spec.

        Given:
            The non-test ``_oracle`` module and the conftest-backed fixture.
        When:
            Their import specs are looked up.
        Then:
            The ``_oracle`` module should be findable, anchoring the lane's
            shared-helper layout.
        """
        # Arrange / Act
        spec = importlib.util.find_spec("tests.integration._oracle")

        # Assert
        assert spec is not None
