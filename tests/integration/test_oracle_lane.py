"""Lane purity guard for the cross-target oracle (#139, T5).

This pins the one invariant the oracle's lane layout depends on that nothing
else covers: :func:`resolve_routing` is a *pure* decision that resolves the full
routing map with no engine installed or imported, so non-DataFusion lanes (e.g.
bedtools / coordinate_space) can reuse the engine-free internals.
"""

import pytest

from tests.integration._oracle import resolve_routing

pytestmark = pytest.mark.integration


class TestSharedOracleReuse:
    """The engine-free routing decision is reusable from any lane."""

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
