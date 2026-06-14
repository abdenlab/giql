"""Cross-target result-oracle fixture for the operator-expansion epic (#137).

The oracle is the verification backbone every later operator migration
(#140-#144) consumes: it transpiles one GIQL query for each registered target
(``Generic`` via ``dialect=None``, ``DuckDB`` via ``dialect="duckdb"``,
``DataFusion`` via ``dialect="datafusion"``), executes each target's SQL on the
engine that target is meant for, and asserts that every target's result set is
identical to the others and equal to an expected set.

It compares ROWS, not SQL strings, on purpose: the DuckDB IEJoin SQL
(``SET VARIABLE`` / ``getvariable``) is DuckDB-only and will not even parse on
DataFusion. Result identity is the contract that survives that divergence.

Engine routing per target:

* ``"generic"`` (``dialect=None``) emits portable SQL -> executed on DataFusion.
* ``"datafusion"`` emits the DataFusion target SQL -> executed on DataFusion.
* ``"duckdb"`` emits the DuckDB-specific IEJoin SQL -> executed on DuckDB.

Adding a new operator case is trivial: write one ``oracle(...)`` call with the
query, the table data, and the expected rows. Cases where a target cannot run
an operator yet (e.g. correlated ``LATERAL`` on DataFusion) restrict the run to
the supported targets via the ``targets`` argument.

The engine-free internals live in :mod:`tests.integration._oracle` (an
importable, non-test module) so they can be unit-tested without engines; this
fixture is a thin wrapper around them.

This lane is deliberately kept an **explicit-case fixture**. It does NOT use
the SDLC guide's pairwise / Scenario-model / Hypothesis apparatus: each case is
a hand-written ``(query, data, expected)`` triple so a regression points at a
named operator scenario, and the cross-engine differential oracle -- not a
generator -- is what makes the cases load-bearing. The coordinate-encoding
sweep (``coordinate_space`` / #145) and strand (the schema here is
``chrom``/``start``/``end`` only) are out of scope by design.
"""

from __future__ import annotations

import pytest

from giql import Table
from giql import transpile

from ._oracle import ENGINE_RUNNERS
from ._oracle import assert_cross_target
from ._oracle import normalize
from ._oracle import resolve_routing

# Default per-target schema: GIQL's default interval columns. Reserved words
# (``end``) are quoted by the loaders below, mirroring ``test_binned_join.py``.
DEFAULT_COLUMNS: tuple[tuple[str, str], ...] = (
    ("chrom", "utf8"),
    ("start", "int64"),
    ("end", "int64"),
)

# Default target spread: every target executed on its intended engine.
DEFAULT_TARGETS: tuple[str, ...] = ("generic", "datafusion", "duckdb")


@pytest.fixture
def cross_target_oracle():
    """Return a callable asserting cross-target result identity for a GIQL query.

    The callable signature::

        oracle(
            query,
            *,
            expected,
            tables=None,
            columns=DEFAULT_COLUMNS,
            targets=DEFAULT_TARGETS,
            engines=None,
            **table_data,
        )

    ``query``
        The GIQL query string.
    ``expected``
        The expected result as an iterable of row tuples (order-free; compared
        as a sorted multiset).
    ``tables``
        Optional list of table names or :class:`giql.Table` specs passed to
        :func:`giql.transpile`. Defaults to the keys of ``table_data`` using the
        default column mapping.
    ``columns``
        The ``(name, kind)`` schema (``kind`` in ``{"utf8", "int64"}``) used to
        register ``table_data`` into both engines. Defaults to the GIQL default
        ``chrom``/``start``/``end`` interval schema, where ``chrom`` is ``utf8``
        and ``start``/``end`` are ``int64``.
    ``targets``
        The target names to exercise. Defaults to all three. Restrict this when
        a target cannot run the operator on its engine yet (e.g. NEAREST's
        correlated LATERAL has no DataFusion physical plan).
    ``engines``
        Optional ``{target: engine}`` overrides for the default routing
        (``generic``/``datafusion`` -> DataFusion, ``duckdb`` -> DuckDB). Use
        this when a target's *portable* SQL must be executed on a different
        engine for the case at hand -- e.g. routing ``generic`` to DuckDB so a
        LATERAL-based operator can still be compared against the duckdb target.
    ``**table_data``
        ``name=[row, ...]`` table contents, where each row is a tuple matching
        ``columns``.

    Returns the per-target normalized result so callers may make extra
    assertions.
    """

    def _oracle(
        query: str,
        *,
        expected,
        tables=None,
        columns=DEFAULT_COLUMNS,
        targets=DEFAULT_TARGETS,
        engines=None,
        **table_data,
    ) -> dict[str, list[tuple]]:
        if tables is None:
            tables = [_default_table(name, columns) for name in table_data]

        routing = resolve_routing(targets, engines)

        expected_rows = normalize(expected)
        results: dict[str, list[tuple]] = {}
        sql_by_target: dict[str, str] = {}

        # Collect ALL target results first; defer every assertion to
        # assert_cross_target so the first divergent target can no longer
        # short-circuit the cross-target comparison (Finding 1).
        for target in targets:
            engine, dialect = routing[target]
            pytest.importorskip("duckdb" if engine == "duckdb" else "datafusion")
            if engine == "datafusion":
                pytest.importorskip("pyarrow")

            sql = transpile(query, tables=tables, dialect=dialect)
            sql_by_target[target] = sql
            results[target] = ENGINE_RUNNERS[engine](sql, table_data, columns)

        assert_cross_target(results, expected_rows, sql_by_target)
        return results

    return _oracle


def _default_table(name, columns) -> Table:
    """Build a :class:`giql.Table` mapping the canonical interval columns.

    The default chrom/start/end column names match :class:`giql.Table`'s own
    defaults; any extra columns ride along as plain data columns the query can
    still project. Custom column schemas should pass an explicit ``tables`` list
    to the oracle instead.
    """
    names = {col for col, _ in columns}
    kwargs = {}
    if "chrom" in names:
        kwargs["chrom_col"] = "chrom"
    if "start" in names:
        kwargs["start_col"] = "start"
    if "end" in names:
        kwargs["end_col"] = "end"
    return Table(name, **kwargs)
