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
from giql.generators import BaseGIQLGenerator
from giql.resolver import resolve_operator_refs
from giql.table import Table
from giql.table import Tables
from giql.targets import resolve_target
from giql.transformer import ClusterTransformer
from giql.transformer import IntersectsBinnedJoinTransformer
from giql.transformer import IntersectsDuckDBIEJoinTransformer
from giql.transformer import MergeTransformer


@overload
def transpile(
    giql: str,
    tables: list[str | Table] | None = None,
    *,
    dialect: Literal["datafusion"] | None = None,
    intersects_bin_size: int | None = None,
) -> str: ...


@overload
def transpile(
    giql: str,
    tables: list[str | Table] | None = None,
    *,
    dialect: Literal["duckdb"],
    intersects_bin_size: None = None,
) -> str: ...


def transpile(
    giql: str,
    tables: list[str | Table] | None = None,
    *,
    dialect: Literal["duckdb", "datafusion"] | None = None,
    intersects_bin_size: int | None = None,
) -> str:
    """Transpile a GIQL query to SQL.

    Parses the GIQL syntax and converts it to standard SQL-92 compatible
    output (uses LATERAL joins where needed for operations like NEAREST).

    Parameters
    ----------
    giql : str
        The GIQL query string containing genomic extensions like
        INTERSECTS, CONTAINS, WITHIN, CLUSTER, MERGE, NEAREST, or DISJOIN.
    tables : list[str | :class:`Table`] | None
        Table configurations. Strings use default column mappings
        (chrom, start, end, strand). :class:`Table` objects provide
        custom column name mappings.
    dialect : Literal["duckdb", "datafusion"] | None
        Optional target engine. Resolves to a :class:`giql.targets.Target`
        carrying the engine's capability set; ``None`` selects the generic
        portable target. When set to ``"duckdb"``, column-to-column
        ``INTERSECTS`` joins (INNER, SEMI, or ANTI) are transpiled into a
        per-chromosome dynamic-SQL pattern (``SET VARIABLE`` +
        ``query(getvariable(...))``) that DuckDB plans through its
        range-join family (``IE_JOIN`` / ``PIECEWISE_MERGE_JOIN``); this
        IEJoin plan is mutually exclusive with ``intersects_bin_size``.
        ``"datafusion"`` and ``None`` use the generic binned equi-join path
        and accept ``intersects_bin_size``. Hard-error projection shapes
        raise ``ValueError`` at transpile time; see the performance guide
        for the full enumeration.
    intersects_bin_size : int | None
        Bin size for INTERSECTS equi-join optimization. When a query
        contains a full-table column-to-column INTERSECTS join, the
        transpiler rewrites it as a binned equi-join for performance.
        Defaults to 10,000 if not specified. Cannot be combined with
        ``dialect="duckdb"``.

    Returns
    -------
    str
        The transpiled SQL query.

    Raises
    ------
    ValueError
        If the query cannot be parsed or transpiled, if ``dialect`` is
        unknown, or if ``dialect="duckdb"`` and ``intersects_bin_size``
        are both set.

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

    Binned equi-join with custom bin size::

        sql = transpile(
            "SELECT a.*, b.* FROM peaks a JOIN genes b "
            "ON a.interval INTERSECTS b.interval",
            tables=["peaks", "genes"],
            intersects_bin_size=100000,
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
    if uses_iejoin and intersects_bin_size is not None:
        raise ValueError(
            f"intersects_bin_size has no effect with dialect={target.name!r}; "
            f"the {target.name} target uses an IEJoin per-partition plan instead "
            "of the binned equi-join. Pass one or the other, not both."
        )

    # The INTERSECTS join-rewrite override (governs the three uses below).
    #
    # A *target-specific* ``(target, Intersects)`` registry entry — the public
    # extension hook, matched by ``ExpanderRegistry.has_override`` — takes over the
    # INTERSECTS join rewrite entirely. ``has_override`` deliberately matches only
    # an *exact non-generic* entry: the built-in ``(GenericTarget(), Intersects)``
    # predicate expander is NOT a join-strategy override (it only renders the
    # residual / literal-range predicates the join transformers leave behind), so
    # it must not disable the join rewrite. When an override is present, all three
    # built-in join paths below are bypassed for that target so the INTERSECTS node
    # flows untouched into ExpandOperators, which dispatches it to that expander:
    #   1. ``intersects_bin_size`` is rejected — it only configures the built-in
    #      binned transformer the override supersedes (rejected here, parallel to
    #      the iejoin rejection above, rather than silently dropped);
    #   2. the DuckDB IEJoin short-circuit is skipped (the registry-deferral the
    #      IEJoin early-return used to preclude, #141);
    #   3. the binned-join transformer is skipped.
    target_overrides_intersects = REGISTRY.has_override(target, Intersects)
    if target_overrides_intersects and intersects_bin_size is not None:
        raise ValueError(
            "intersects_bin_size has no effect when a target-specific "
            f"(target={target.name!r}, Intersects) expander is registered; that "
            "expander supersedes the built-in binned join transformer the bin "
            "size configures. Pass one or the other, not both."
        )

    tables_container = _build_tables(tables)

    with _reraise_as_value_error("Parse error", query=giql):
        ast = parse_one(giql, dialect=GIQLDialect)

    # The column-to-column INTERSECTS *join* rewrites are capability-gated
    # pre-pass transformers (epic #137, #141): the target's
    # ``range_join_strategy`` selects the DuckDB IEJoin plan or the generic
    # binned equi-join. They run on the raw parsed AST (before resolution, which
    # rewrites the genomic column name) and consume a column-to-column INTERSECTS
    # *join* so it never reaches the predicate expander; a literal-range or
    # residual column-to-column INTERSECTS *predicate* survives to pass 3.
    # ``target_overrides_intersects`` (computed above) gates whether these
    # built-in join paths run — see its definition for the override rationale.

    # Falls back to the binned plan for unsupported shapes — see
    # IntersectsDuckDBIEJoinTransformer.transform_to_sql for the complete
    # fallback set. The IEJoin transformer emits a whole-query string, so when it
    # produces output it must short-circuit the AST pipeline. This never skips an
    # INTERSECTS that pass 3 would expand: ``_classify_extras`` forces the binned
    # fallback (returning None here) for any query carrying a residual INTERSECTS
    # alongside the join, so a query the IEJoin transformer accepts has no residual
    # INTERSECTS left for the expander. (A residual CONTAINS/WITHIN/ANY beside an
    # IEJoin is a pre-existing IEJoin limitation that errors identically on main.)
    if uses_iejoin and not target_overrides_intersects:
        duckdb_transformer = IntersectsDuckDBIEJoinTransformer(tables_container)
        with _reraise_as_value_error("Transformation error"):
            duckdb_sql = duckdb_transformer.transform_to_sql(ast)
        if duckdb_sql is not None:
            return duckdb_sql

    merge_transformer = MergeTransformer(tables_container)
    cluster_transformer = ClusterTransformer(tables_container)
    generator = BaseGIQLGenerator(tables=tables_container)

    with _reraise_as_value_error("Transformation error"):
        # Reaching here with an iejoin target means the IEJoin transformer
        # declined the query (returned None) and fell back to the binned plan,
        # exactly as before. ``intersects_bin_size`` is rejected up front for
        # iejoin targets, so the binned transformer always sees its default there.
        if not target_overrides_intersects:
            intersects_transformer = IntersectsBinnedJoinTransformer(
                tables_container,
                bin_size=intersects_bin_size,
            )
            ast = intersects_transformer.transform(ast)
        ast = merge_transformer.transform(ast)
        ast = cluster_transformer.transform(ast)

    # Pass 1 of the normalization pipeline (epic #114): attach resolution
    # metadata to every GIQL operator slot ahead of generation. DISJOIN's
    # emitter consumes this metadata (step 2); the remaining operators still
    # use the generator's legacy resolver paths until their ports land.
    with _reraise_as_value_error("Resolution error"):
        ast = resolve_operator_refs(ast, tables_container)

    # Pass 2 of the normalization pipeline (epic #114): synthesize canonical
    # __giql_canon_* wrapper CTEs for non-canonical interval operands of operators
    # that opt in via GIQL_CANONICALIZE; those operators are rewritten here, and
    # operators that do not opt in are left untouched.
    with _reraise_as_value_error("Canonicalization error"):
        ast = canonicalize_coordinates(ast)

    # Pass 3 of the normalization pipeline (epic #137): replace each GIQL operator
    # node that opts in (GIQL_EXPAND) and resolves a registered expander with the
    # AST that expander produces for the active target. Operators that are
    # unflagged or resolve no expander are left untouched and the generator renders
    # them via their legacy ``*_sql`` emitter as before.
    expand_pass = ExpandOperators(target, tables_container)
    with _reraise_as_value_error("Expansion error"):
        ast = expand_pass.transform(ast)

    with _reraise_as_value_error("Transpilation error"):
        sql = generator.generate(ast)

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
    and generator propagate verbatim (so the dialect's targeted error messages
    survive the boundary) while wrapping unexpected exceptions in a uniform
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
