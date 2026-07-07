"""Transpile GIQL queries to SQL.

This module provides the main entry point for transpiling GIQL queries
to standard SQL.
"""

from contextlib import contextmanager
from typing import Iterator
from typing import Literal
from typing import overload

from sqlglot import parse_one

import giql.expanders  # noqa: F401
from giql.canonicalizer import canonicalize_coordinates
from giql.dialect import GIQLDialect
from giql.expander import REGISTRY
from giql.expander import ExpandOperators
from giql.expressions import Intersects
from giql.resolver import resolve_operator_refs
from giql.table import Table
from giql.table import Tables
from giql.targets import resolve_target
from giql.transformer import IntersectsDuckDBIEJoinTransformer


@overload
def transpile(
    giql: str,
    tables: list[str | Table] | None = None,
    *,
    dialect: Literal["datafusion"] | None = None,
) -> str: ...


@overload
def transpile(
    giql: str,
    tables: list[str | Table] | None = None,
    *,
    dialect: Literal["duckdb"],
) -> str: ...


# Widened overload for arbitrary custom-target dialect names (the registry
# plugin hub). It subsumes the ``Literal["duckdb"]`` overload above, which is
# kept deliberately so editors still autocomplete the ``"duckdb"`` literal
# (symmetric with the ``Literal["datafusion"]`` overload).
@overload
def transpile(
    giql: str,
    tables: list[str | Table] | None = None,
    *,
    dialect: str,
) -> str: ...


def transpile(
    giql: str,
    tables: list[str | Table] | None = None,
    *,
    dialect: str | None = None,
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
    dialect : str | None
        Optional target engine. Resolves to a :class:`giql.targets.Target`
        carrying the engine's capability set; ``None`` selects the generic
        portable target, ``"duckdb"`` / ``"datafusion"`` the built-in targets,
        and any other name a custom :class:`giql.targets.Target` registered on
        the plugin hub (see :func:`giql.expander.register` /
        :meth:`giql.expander.ExpanderRegistry.register_target`). When set to
        ``"duckdb"``, column-to-column
        ``INTERSECTS`` joins (INNER, SEMI, or ANTI) are transpiled into a
        per-chromosome dynamic-SQL pattern (``SET VARIABLE`` +
        ``query(getvariable(...))``) that DuckDB plans through its
        range-join family (``IE_JOIN`` / ``PIECEWISE_MERGE_JOIN``). Shapes
        the IEJoin path declines (LEFT/RIGHT/FULL joins, self-joins, multiple
        INTERSECTS, extra predicates, non-base operands, 3+ tables) fall
        through to the generic naive overlap predicate. ``"datafusion"`` and
        ``None`` always emit that naive predicate — a plain
        ``ON a.chrom = b.chrom AND a.start < b.end AND b.start < a.end``
        condition the engine's own optimizer plans as a range join — for both
        inner and outer column-to-column INTERSECTS joins. Hard-error
        projection shapes raise ``ValueError`` at transpile time; see the
        performance guide for the full enumeration. The target's capabilities
        also choose the
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
        The transpiled SQL query.

    Raises
    ------
    ValueError
        If the query cannot be parsed or transpiled, or if ``dialect`` is
        unknown.

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

    DuckDB IEJoin dialect (column-to-column INNER/SEMI/ANTI JOIN only;
    projections must be qualified)::

        sql = transpile(
            "SELECT a.chrom, a.start, b.start "
            "FROM peaks a JOIN genes b ON a.interval INTERSECTS b.interval",
            tables=["peaks", "genes"],
            dialect="duckdb",
        )
    """
    target = resolve_target(dialect)
    uses_iejoin = target.capabilities.range_join_strategy == "iejoin"

    # The INTERSECTS join-rewrite override.
    #
    # A *target-specific* ``(target, Intersects)`` registry entry — the public
    # extension hook, matched by ``ExpanderRegistry.has_override`` — takes over the
    # INTERSECTS join rewrite entirely. ``has_override`` deliberately matches only
    # an *exact non-generic* entry: the built-in ``(GenericTarget(), Intersects)``
    # predicate expander is NOT a join-strategy override — it renders the naive
    # overlap predicate that IS the generic join plan, so it must not disable the
    # IEJoin short-circuit. When an override is present, the DuckDB IEJoin
    # short-circuit below is skipped (the registry-deferral the IEJoin early-return
    # used to preclude, #141) so the INTERSECTS node flows untouched into
    # ExpandOperators, which dispatches it to that expander.
    target_overrides_intersects = REGISTRY.has_override(target, Intersects)

    tables_container = _build_tables(tables)

    with _reraise_as_value_error("Parse error", query=giql):
        ast = parse_one(giql, dialect=GIQLDialect)

    # The DuckDB IEJoin plan is the one remaining capability-gated pre-pass
    # transformer (epic #137, #141): ``range_join_strategy == "iejoin"`` selects
    # it. It runs on the raw parsed AST (before resolution, which rewrites the
    # genomic column name) and consumes a column-to-column INTERSECTS *join* so it
    # never reaches the predicate expander. Every other target — and every shape
    # the IEJoin path declines — leaves the column-to-column INTERSECTS in place
    # for pass 3, where the ``(GenericTarget, Intersects)`` expander renders it as
    # the naive overlap predicate (a plain ``ON`` condition the engine plans as a
    # range join). Literal-range and residual INTERSECTS predicates flow the same
    # way. ``target_overrides_intersects`` (computed above) gates whether the
    # IEJoin path runs — see its definition for the override rationale.
    #
    # The IEJoin transformer emits a whole-query string, so when it produces output
    # it must short-circuit the AST pipeline. This never skips an INTERSECTS that
    # pass 3 would expand: ``_classify_extras`` declines (returning None here) for
    # any query carrying a residual INTERSECTS alongside the join, so a query the
    # IEJoin transformer accepts has no residual INTERSECTS left for the expander.
    # (A residual CONTAINS/WITHIN/ANY beside an IEJoin is a pre-existing IEJoin
    # limitation that errors identically on main.)
    #
    # CLUSTER and MERGE used to be rewritten as pre-pass transformers too; they are
    # now expanded in pass 3 (ExpandOperators) by giql.expanders.cluster / .merge
    # (#144), and the generic binned equi-join was dropped in favor of the naive
    # predicate (#167).
    if uses_iejoin and not target_overrides_intersects:
        duckdb_transformer = IntersectsDuckDBIEJoinTransformer(tables_container)
        with _reraise_as_value_error("Transformation error"):
            duckdb_sql = duckdb_transformer.transform_to_sql(ast)
        if duckdb_sql is not None:
            return duckdb_sql

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
        else:
            container.register(item.name, item)

    return container


@contextmanager
def _reraise_as_value_error(prefix: str, query: str | None = None) -> Iterator[None]:
    """Re-raise non-:class:`ValueError` exceptions as :class:`ValueError` with *prefix*.

    Lets user-facing :class:`ValueError`\\s from the parser, transformer chain,
    expander pass, and stock serializer propagate verbatim (so the dialect's
    targeted error messages survive the boundary) while wrapping unexpected
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
