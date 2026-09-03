"""Transpile GIQL queries to SQL.

This module provides the main entry point for transpiling GIQL queries
to standard SQL.
"""

from contextlib import contextmanager
from typing import Iterator

from sqlglot import parse_one

import giql.expanders  # noqa: F401
from giql.canonicalizer import canonicalize_coordinates
from giql.dialect import GIQLDialect
from giql.expander import ExpandOperators
from giql.resolver import resolve_operator_refs
from giql.table import Table
from giql.table import Tables
from giql.targets import DialectName
from giql.targets import Target
from giql.targets import resolve_target


def transpile(
    giql: str,
    tables: list[str | Table] | None = None,
    *,
    dialect: DialectName | str | Target | None = None,
) -> str:
    """Transpile a GIQL query to SQL.

    Parses the GIQL syntax and converts it to portable SQL for the active target
    (uses LATERAL joins where needed for operations like NEAREST). The output is
    not strictly SQL-92: depending on the target it may use engine extensions
    such as ``LATERAL`` or ``SELECT * EXCEPT`` (see the ``dialect`` parameter).

    Parameters
    ----------
    giql : str
        The GIQL query string containing genomic extensions like
        INTERSECTS, CONTAINS, WITHIN, CLUSTER, MERGE, NEAREST, or DISJOIN.
    tables : list[str | :class:`Table`] | None
        Table configurations. Strings use default column mappings
        (chrom, start, end, strand). :class:`Table` objects provide
        custom column name mappings.
    dialect : DialectName | str | Target | None
        Optional target engine, as either a :class:`giql.targets.Target` or the
        name of one. A Target carries the engine's capability set; ``None``
        selects the generic portable target, ``"duckdb"`` / ``"datafusion"`` the
        built-in targets, and any other name a custom
        :class:`giql.targets.Target` registered on the plugin hub (see
        :func:`giql.expander.register` /
        :meth:`giql.expander.ExpanderRegistry.register_target`). Passing the
        instance directly is the way to use a one-off Target without registering
        it under a name, since registration is process-global and outlives the
        call. :class:`giql.targets.GenericTarget` is accepted as an instance but
        has no selectable name -- ``None`` is the one public spelling for it. When set to
        ``"duckdb"``, column-to-column ``INTERSECTS`` joins (INNER, SEMI, ANTI,
        and LEFT / RIGHT outer) are transpiled into a per-chromosome dynamic-SQL
        pattern (``SET VARIABLE`` + ``query(getvariable(...))``) that DuckDB
        plans through its range-join family (``IE_JOIN`` /
        ``PIECEWISE_MERGE_JOIN``). A LEFT or RIGHT outer join is decomposed into
        an INNER half, a ``NOT EXISTS`` half carrying the preserved side's
        unmatched rows, and a third branch carrying its NULL-chromosome rows,
        combined with ``UNION ALL``; the first two reach the fast operator and
        the third is a filtered scan. Because the per-chromosome form emits one
        branch per distinct chromosome, it degrades to a single non-partitioned
        query above a partition-count ceiling, decided at execution time: the
        partitioned form wins by orders of magnitude on a primary assembly and
        loses by orders of magnitude on a scaffold-level one. Shapes the IEJoin
        path declines fall through to the generic naive overlap predicate; the
        performance guide enumerates them, and that enumeration is deliberately
        not repeated here, because a second copy drifts. ``"datafusion"`` and
        ``None`` always emit that naive predicate — a plain
        ``ON a.chrom = b.chrom AND a.start < b.end AND b.start < a.end``
        condition the engine's own optimizer plans as a range join — for both
        inner and outer column-to-column INTERSECTS joins. For the INNER /
        SEMI / ANTI shapes, a projection the rewrite cannot attribute to a side
        raises ``ValueError`` at transpile time. The outer-join form instead
        declines silently to the naive predicate for every projection it cannot
        split, so ``dialect="duckdb"`` never changes an outer join's result —
        only its plan. See the performance guide for the full enumeration.
        This dialect also names SELECT-list items the other targets leave to
        the engine: an unaliased scalar expression or aggregate is emitted as
        ``__giql_expr_<n>`` / ``__giql_agg_<n>``, because the rebuilt text
        DuckDB would derive a name from refers to inner columns the caller
        never wrote. Alias such an item to pin its result-column name across
        targets.

        The target's capabilities also choose the
        coordinate-canonicalization emit form for a non-canonically-encoded
        table: ``"duckdb"`` emits ``SELECT * REPLACE (...)``, while the generic
        (``None``) and ``"datafusion"`` targets emit the portable
        ``SELECT * EXCEPT (...), <start>, <end>`` form, which runs on
        ``* EXCEPT``-capable engines (the DataFusion family) but **not** on
        DuckDB — pass ``dialect="duckdb"`` for DuckDB-runnable output. Tables in
        the canonical 0-based half-open encoding are unaffected (they emit
        portable SQL on every target).

    Returns
    -------
    str
        The transpiled SQL. Every target but ``"duckdb"`` returns a single
        statement, as does ``"duckdb"`` for any shape its IEJoin path declines.
        An accelerated column-to-column ``INTERSECTS`` join instead returns a
        multi-statement script: one ``SET VARIABLE`` statement, or two for a
        decomposed outer join, followed by the final ``SELECT``. Execute it on a
        single connection -- a driver that splits statements, or forwards only
        the last, drops the variable the ``SELECT`` reads and yields empty or
        NULL results. Variable names are derived from the content they hold, so
        re-running one query shape rebinds its own names rather than leaving a
        fresh pair behind per call.

    Raises
    ------
    ValueError
        If the query cannot be parsed or transpiled, or if ``dialect`` is
        unknown.
    TypeError
        If an entry in ``tables`` is neither a ``str`` nor a
        :class:`giql.table.Table`.

    Examples
    --------
    Basic usage with default column mappings::

        sql = transpile(
            "SELECT * FROM peaks WHERE interval INTERSECTS 'chr1:1000-2000'",
            tables=["peaks"],
        )

    Custom :class:`Table` configuration::

        sql = transpile(
            "SELECT * FROM peaks WHERE interval INTERSECTS 'chr1:1000-2000'",
            tables=[
                Table(
                    "peaks",
                    genomic_col="interval",
                    chrom_col="chrom",
                    start_col="start",
                    end_col="end",
                )
            ],
        )

    Column-to-column INTERSECTS join (naive overlap predicate; inner or
    outer, planned as a range join by the target engine)::

        sql = transpile(
            "SELECT a.*, b.* FROM peaks a JOIN genes b "
            "ON a.interval INTERSECTS b.interval",
            tables=["peaks", "genes"],
        )

    DuckDB IEJoin dialect (column-to-column INNER/SEMI/ANTI/LEFT/RIGHT JOIN).
    INNER, SEMI and ANTI require qualified projections and raise otherwise; the
    outer-join forms decline to the naive predicate instead::

        sql = transpile(
            "SELECT a.chrom, a.start, b.start "
            "FROM peaks a JOIN genes b ON a.interval INTERSECTS b.interval",
            tables=["peaks", "genes"],
            dialect="duckdb",
        )
    """
    target = resolve_target(dialect)

    tables_container = _build_tables(tables)

    with _reraise_as_value_error("Parse error", query=giql):
        ast = parse_one(giql, dialect=GIQLDialect)

    # Every INTERSECTS join strategy is now registry-driven (epic #137, #169). A
    # column-to-column INTERSECTS *join* flows untouched through the pipeline into
    # pass 3, where ``ExpandOperators`` dispatches it to its registered expander:
    # the ``(GenericTarget, Intersects)`` expander renders the naive overlap
    # predicate (a plain ``ON`` condition the engine plans as a range join) on
    # ``None`` / ``"datafusion"``, and the ``(DuckDBTarget, Intersects)`` override
    # (giql.expanders.intersects_duckdb) rewrites it into the per-chromosome IEJoin
    # plan on ``"duckdb"``, deferring to that same naive predicate for the shapes it
    # declines. Literal-range and residual INTERSECTS predicates flow the same way.
    # There is no capability-gated pre-pass anymore — the former DuckDB IEJoin
    # pre-pass (and the CLUSTER / MERGE pre-passes, #144) were all relocated into
    # the operator-expander registry, and the generic binned equi-join was dropped
    # in favor of the naive predicate (#167).

    # Pass 1 of the normalization pipeline (epic #114): attach resolution
    # metadata to every GIQL operator slot ahead of generation. Every migrated
    # operator's expander consumes this metadata in pass 3 (CLUSTER/MERGE carry an
    # empty resolution, deriving their columns from the FROM table instead).
    with _reraise_as_value_error("Resolution error"):
        ast = resolve_operator_refs(ast, tables_container)

    # Pass 2 of the normalization pipeline (epic #114): synthesize canonical
    # __giql_canon_* wrapper CTEs for non-canonical interval operands of operators
    # that opt in via GIQL_CANONICALIZE; those operators are rewritten here, and
    # operators that do not opt in are left untouched. The active target's
    # capabilities choose the wrapper's emit strategy (* REPLACE vs the portable
    # * EXCEPT form — epic #137 / #145).
    with _reraise_as_value_error("Canonicalization error"):
        ast = canonicalize_coordinates(ast, target.capabilities)

    # Pass 3 of the normalization pipeline (epic #137): replace every GIQL operator
    # node with the AST its registered expander produces for the active target.
    # With every operator migrated, this pass fully consumes the GIQL dialect —
    # nothing GIQL-specific survives into serialization.
    expand_pass = ExpandOperators(target, tables_container)
    with _reraise_as_value_error("Expansion error"):
        ast = expand_pass.transform(ast)

    # Serialize the now-standard AST with the stock sqlglot serializer for the
    # active target (epic #137, #146). The target's ``sqlglot_dialect`` selects the
    # engine's serialization (``None`` is sqlglot's portable default); there is no
    # custom GIQL generator anymore.
    with _reraise_as_value_error("Transpilation error"):
        sql = ast.sql(dialect=target.sqlglot_dialect)

    return sql


def _build_tables(tables: list[str | Table] | None) -> Tables:
    """Build a :class:`Tables` container from table specifications.

    Parameters
    ----------
    tables : list[str | :class:`Table`] | None
        Table specifications. Strings use default column mappings.
        :class:`Table` objects provide custom column mappings.

    Returns
    -------
    Tables
        Container with all tables registered.
    """
    container = Tables()

    if tables is None:
        return container

    for item in tables:
        if isinstance(item, str):
            container.register(item, Table(item))
        elif isinstance(item, Table):
            container.register(item.name, item)
        else:
            # Reject here rather than duck-typing ``item.name``. An object that
            # merely looks like a Table reaches pass 1 and fails there with an
            # AttributeError naming an internal column attribute, which tells a
            # caller nothing about which argument was wrong.
            raise TypeError(
                f"tables entries must be str or Table, got {type(item).__name__}"
            )

    return container


@contextmanager
def _reraise_as_value_error(prefix: str, query: str | None = None) -> Iterator[None]:
    """Re-raise non-:class:`ValueError` exceptions as :class:`ValueError` with *prefix*.

    Lets user-facing :class:`ValueError`\\s from the parser, the resolution /
    canonicalization / operator-expansion passes, and the stock serializer
    propagate verbatim (so the dialect's targeted error messages survive the
    boundary) while wrapping unexpected
    exceptions in a uniform
    :class:`ValueError` prefixed with the stage name. When *query* is supplied,
    the original input is appended to the message so parse errors retain the
    offending text.
    """
    try:
        yield
    except ValueError:
        raise
    except Exception as e:
        msg = f"{prefix}: {e}"
        if query is not None:
            msg += f"\nQuery: {query}"
        raise ValueError(msg) from e
