"""Engine-free internals for the cross-target result oracle (epic #137, #139).

This is an importable, *non-test* module: it holds the pure helpers the
``cross_target_oracle`` fixture (``tests/integration/conftest.py``) wraps, so
they can be unit-tested directly without spinning up DuckDB or DataFusion.

Layout:

* :func:`normalize` / :func:`scalar` -- coerce engine rows into a sorted,
  order-free, type-normalized multiset for comparison.
* :func:`resolve_routing` -- the *pure* routing decision (target -> engine +
  transpile dialect), pulled out of the run loop so it is unit-testable and so
  unknown targets / engines raise a clear error rather than a bare ``KeyError``
  (Finding 7).
* :func:`assert_cross_target` -- the comparison core (Finding 1): assert the
  engines agree *with each other* first, then assert the agreed result equals
  ``expected``. Both assertions are genuinely reachable.
* :func:`register_record_batches` -- the single pyarrow loader dance reused by
  :func:`run_datafusion` and the #132 ``datafusion_ctx`` fixture (Finding 8).
* :func:`run_duckdb` / :func:`run_datafusion` -- the per-engine runners.
"""

from __future__ import annotations

import math

# Target -> the engine the target's SQL is meant to execute on.
TARGET_ENGINE: dict[str, str] = {
    "generic": "datafusion",
    "datafusion": "datafusion",
    "duckdb": "duckdb",
}

# transpile() dialect argument per target name.
TARGET_DIALECT: dict[str, str | None] = {
    "generic": None,
    "datafusion": "datafusion",
    "duckdb": "duckdb",
}

# The engines a target may be routed onto.
KNOWN_ENGINES: frozenset[str] = frozenset({"duckdb", "datafusion"})


def scalar(value):
    """Coerce an engine-specific scalar to a comparable plain Python value.

    ``None`` and any NaN (whether a ``float`` subclass or a numpy scalar that
    only reveals its NaN-ness after ``.item()``) both collapse to ``None`` so
    nullable integers from DataFusion compare equal to DuckDB's ``NULL`` and
    can never reach the ``(is_none, value)`` sort key as a bare NaN.
    """
    if value is None:
        return None
    # Fast path: a plain float NaN.
    if isinstance(value, float) and math.isnan(value):
        return None
    if hasattr(value, "item"):
        try:
            value = value.item()
        except (ValueError, AttributeError):
            return value
        # Re-check NaN AFTER .item(): a numpy NaN unwraps to a float NaN that
        # must not survive as a sort key (Finding 6).
        if value is None:
            return None
        if isinstance(value, float) and math.isnan(value):
            return None
    return value


def normalize(rows) -> list[tuple]:
    """Return a row multiset as a sorted list of tuples for order-free compare.

    A list (not a set) preserves multiplicity so duplicate-row bugs surface;
    sorting strips engine result ordering. Values are coerced to plain Python
    scalars so DuckDB and DataFusion (pandas) rows compare equal.
    """
    out = [tuple(scalar(v) for v in row) for row in rows]
    return sorted(out, key=lambda r: tuple((v is None, v) for v in r))


def resolve_routing(targets, engines=None) -> dict[str, tuple[str, str | None]]:
    """Resolve each target to its ``(engine, transpile_dialect)`` pair.

    This is the pure routing decision lifted out of the oracle run loop so it
    can be unit-tested without engines (Finding 7).

    Parameters
    ----------
    targets : iterable of str
        Target names to route (e.g. ``"generic"``, ``"datafusion"``,
        ``"duckdb"``).
    engines : dict[str, str] | None
        Optional ``{target: engine}`` overrides for the default routing. Use
        this to execute a target's portable SQL on a different engine (e.g.
        route ``generic`` to DuckDB so a LATERAL operator can still be
        compared).

    Returns
    -------
    dict[str, tuple[str, str | None]]
        ``{target: (engine, transpile_dialect)}`` for each requested target.

    Raises
    ------
    ValueError
        If a target name is not a known target, or an override routes to an
        unknown engine.
    """
    engine_for = dict(TARGET_ENGINE)
    if engines:
        for target, engine in engines.items():
            if target not in TARGET_ENGINE:
                raise ValueError(
                    f"Unknown target {target!r} in engines override. "
                    f"Known targets: {sorted(TARGET_ENGINE)}."
                )
            engine_for[target] = engine

    routing: dict[str, tuple[str, str | None]] = {}
    for target in targets:
        if target not in TARGET_ENGINE:
            raise ValueError(
                f"Unknown target {target!r}. Known targets: {sorted(TARGET_ENGINE)}."
            )
        engine = engine_for[target]
        if engine not in KNOWN_ENGINES:
            raise ValueError(
                f"Unknown engine {engine!r} for target {target!r}. "
                f"Known engines: {sorted(KNOWN_ENGINES)}."
            )
        routing[target] = (engine, TARGET_DIALECT[target])
    return routing


def assert_cross_target(results, expected, sql_by_target=None) -> None:
    """Assert every target agrees with the others and with ``expected``.

    The comparison runs in two genuinely-reachable phases (Finding 1):

    1. **Differential oracle** -- assert the engines agree *with each other*,
       independent of ``expected``. With every engine's output already
       collected, a disagreement yields a clear "engines disagree" message
       naming both targets and showing both result sets.
    2. **Expectation check** -- assert the (now agreed) result equals the
       expected multiset.

    Parameters
    ----------
    results : dict[str, list[tuple]]
        ``{target: normalized_rows}`` -- every target's already-normalized
        result. Collected in full *before* any assertion so the first failing
        target no longer short-circuits the cross-target comparison.
    expected : list[tuple]
        The normalized expected multiset.
    sql_by_target : dict[str, str] | None
        Optional ``{target: sql}`` for richer diagnostics.
    """
    sql_by_target = sql_by_target or {}
    run_targets = list(results)

    # Phase 1: engines must agree with each other (the differential oracle).
    if len(run_targets) > 1:
        reference = run_targets[0]
        ref_rows = results[reference]
        for other in run_targets[1:]:
            assert results[other] == ref_rows, (
                f"engines disagree: target {reference!r} returned "
                f"{ref_rows!r} but target {other!r} returned "
                f"{results[other]!r}\n"
                f"{reference} SQL: {sql_by_target.get(reference, '<n/a>')}\n"
                f"{other} SQL: {sql_by_target.get(other, '<n/a>')}"
            )

    # Phase 2: the agreed result must equal the expectation.
    for target in run_targets:
        assert results[target] == expected, (
            f"target {target!r} returned {results[target]!r}, "
            f"expected {expected!r}\n"
            f"SQL: {sql_by_target.get(target, '<n/a>')}"
        )


def register_record_batches(ctx, name, rows, schema, columns) -> None:
    """Register ``rows`` as a record-batch table on a DataFusion context.

    The single pyarrow loader dance shared by :func:`run_datafusion` and the
    #132 ``datafusion_ctx`` fixture (Finding 8). ``columns`` is the
    ``(name, kind)`` schema; ``schema`` is the matching :class:`pyarrow.Schema`.

    Empty ``rows`` need care: ``pa.table(...).to_batches()`` yields an EMPTY
    batch *list* for zero rows, and ``SessionContext.register_record_batches``
    PANICS on an empty partition (``index out of bounds`` in DataFusion's Rust
    core). DuckDB sidesteps this with an ``if rows`` guard in :func:`run_duckdb`;
    DataFusion has no such guard, so this helper synthesizes a single empty
    record batch with the declared schema instead. This keeps empty-input cases
    runnable on DataFusion (T4 pins both the panic and this guard).
    """
    import pyarrow as pa

    arrays = {col: [r[idx] for r in rows] for idx, (col, _kind) in enumerate(columns)}
    batches = pa.table(arrays, schema=schema).to_batches()
    if not batches:
        batches = [pa.RecordBatch.from_pylist([], schema=schema)]
    ctx.register_record_batches(name, [batches])


def arrow_schema(columns):
    """Build a :class:`pyarrow.Schema` from a ``(name, kind)`` column spec."""
    import pyarrow as pa

    fields = [
        (col, pa.utf8() if kind == "utf8" else pa.int64()) for col, kind in columns
    ]
    return pa.schema(fields)


def run_duckdb(sql: str, table_data: dict, columns) -> list[tuple]:
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
        return normalize(conn.execute(sql).fetchall())
    finally:
        conn.close()


def run_datafusion(sql: str, table_data: dict, columns) -> list[tuple]:
    """Register tables in a DataFusion context and return normalized result rows."""
    from datafusion import SessionContext

    schema = arrow_schema(columns)
    ctx = SessionContext()
    for name, rows in table_data.items():
        register_record_batches(ctx, name, rows, schema, columns)
    return normalize(ctx.sql(sql).to_pandas().itertuples(index=False, name=None))


ENGINE_RUNNERS = {
    "duckdb": run_duckdb,
    "datafusion": run_datafusion,
}
