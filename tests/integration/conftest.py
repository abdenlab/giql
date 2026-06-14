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
"""

from __future__ import annotations

import pytest

from giql import Table
from giql import transpile

# Default per-target schema: GIQL's default interval columns. Reserved words
# (``end``) are quoted by the loaders below, mirroring ``test_binned_join.py``.
DEFAULT_COLUMNS: tuple[tuple[str, str], ...] = (
    ("chrom", "utf8"),
    ("start", "int64"),
    ("end", "int64"),
)

# Target -> the engine the target's SQL is meant to execute on.
TARGET_ENGINE: dict[str, str] = {
    "generic": "datafusion",
    "datafusion": "datafusion",
    "duckdb": "duckdb",
}

# Default target spread: every target executed on its intended engine.
DEFAULT_TARGETS: tuple[str, ...] = ("generic", "datafusion", "duckdb")

# transpile() dialect argument per target name.
_TARGET_DIALECT: dict[str, str | None] = {
    "generic": None,
    "datafusion": "datafusion",
    "duckdb": "duckdb",
}


def _normalize(rows) -> list[tuple]:
    """Return a row multiset as a sorted list of tuples for order-free compare.

    A list (not a set) preserves multiplicity so duplicate-row bugs surface;
    sorting strips engine result ordering. Values are coerced to plain Python
    scalars so DuckDB and DataFusion (pandas) rows compare equal.
    """
    out = []
    for row in rows:
        out.append(tuple(_scalar(v) for v in row))
    return sorted(out, key=lambda r: tuple((v is None, v) for v in r))


def _scalar(value):
    """Coerce engine-specific scalars to comparable plain Python values."""
    if value is None:
        return None
    # pandas / numpy nullable integers and NaN surface from DataFusion.
    try:
        import math

        if isinstance(value, float) and math.isnan(value):
            return None
    except (TypeError, ValueError):
        pass
    if hasattr(value, "item"):
        try:
            return value.item()
        except (ValueError, AttributeError):
            return value
    return value


def _run_duckdb(sql: str, table_data: dict, columns) -> list[tuple]:
    """Register tables in an in-memory DuckDB and return normalized result rows."""
    import duckdb

    conn = duckdb.connect(":memory:")
    try:
        for name, rows in table_data.items():
            cols_ddl = ", ".join(
                f'"{col}" {"VARCHAR" if kind == "utf8" else "BIGINT"}'
                for col, kind in columns
            )
            conn.execute(f"CREATE TABLE {name} ({cols_ddl})")
            if rows:
                placeholders = ", ".join("?" for _ in columns)
                conn.executemany(
                    f"INSERT INTO {name} VALUES ({placeholders})",
                    [tuple(r) for r in rows],
                )
        return _normalize(conn.execute(sql).fetchall())
    finally:
        conn.close()


def _run_datafusion(sql: str, table_data: dict, columns) -> list[tuple]:
    """Register tables in a DataFusion context and return normalized result rows."""
    import pyarrow as pa
    from datafusion import SessionContext

    pa_fields = [
        (col, pa.utf8() if kind == "utf8" else pa.int64()) for col, kind in columns
    ]
    schema = pa.schema(pa_fields)

    ctx = SessionContext()
    for name, rows in table_data.items():
        arrays = {
            col: [r[idx] for r in rows] for idx, (col, _kind) in enumerate(columns)
        }
        ctx.register_record_batches(name, [pa.table(arrays, schema=schema).to_batches()])
    return _normalize(ctx.sql(sql).to_pandas().itertuples(index=False, name=None))


_ENGINE_RUNNERS = {
    "duckdb": _run_duckdb,
    "datafusion": _run_datafusion,
}


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
        ``chrom``/``start``/``end`` interval schema.
    ``targets``
        The target names to exercise. Defaults to all three. Restrict this when
        a target cannot run the operator on its engine yet (e.g. NEAREST's
        correlated LATERAL has no DataFusion physical plan).
    ``engines``
        Optional ``{target: engine}`` overrides for the default routing
        (``generic``/``datafusion`` -> DataFusion, ``duckdb`` -> DuckDB). Use
        this when a target's *portable* SQL must be executed on a different
        engine for the case at hand — e.g. routing ``generic`` to DuckDB so a
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

        engine_for = dict(TARGET_ENGINE)
        if engines:
            engine_for.update(engines)

        expected_rows = _normalize(expected)
        results: dict[str, list[tuple]] = {}

        for target in targets:
            engine = engine_for[target]
            importorskip = pytest.importorskip
            importorskip("duckdb" if engine == "duckdb" else "datafusion")
            if engine == "datafusion":
                importorskip("pyarrow")

            sql = transpile(query, tables=tables, dialect=_TARGET_DIALECT[target])
            rows = _ENGINE_RUNNERS[engine](sql, table_data, columns)
            results[target] = rows

            assert rows == expected_rows, (
                f"target {target!r} on engine {engine!r} returned {rows!r}, "
                f"expected {expected_rows!r}\nSQL: {sql}"
            )

        # Cross-target identity: every executed target agrees row-for-row.
        run_targets = list(results)
        if len(run_targets) > 1:
            reference = run_targets[0]
            for other in run_targets[1:]:
                assert results[other] == results[reference], (
                    f"targets {reference!r} and {other!r} disagree: "
                    f"{results[reference]!r} != {results[other]!r}"
                )

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
