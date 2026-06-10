"""Transpile GIQL queries to SQL.

This module provides the main entry point for transpiling GIQL queries
to standard SQL.
"""

from contextlib import contextmanager
from typing import Iterator
from typing import Literal
from typing import overload

from sqlglot import parse_one

from giql.canonicalizer import canonicalize_coordinates
from giql.dialect import GIQLDialect
from giql.generators import BaseGIQLGenerator
from giql.resolver import resolve_operator_refs
from giql.table import Table
from giql.table import Tables
from giql.transformer import ClusterTransformer
from giql.transformer import IntersectsBinnedJoinTransformer
from giql.transformer import IntersectsDuckDBIEJoinTransformer
from giql.transformer import MergeTransformer


@overload
def transpile(
    giql: str,
    tables: list[str | Table] | None = None,
    *,
    dialect: None = None,
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
    dialect: Literal["duckdb"] | None = None,
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
    dialect : Literal["duckdb"] | None
        Optional target dialect. When set to ``"duckdb"``, column-to-column
        ``INTERSECTS`` joins (INNER, SEMI, or ANTI) are transpiled into a
        per-chromosome dynamic-SQL pattern (``SET VARIABLE`` +
        ``query(getvariable(...))``) that DuckDB plans through its
        range-join family (``IE_JOIN`` / ``PIECEWISE_MERGE_JOIN``).
        Mutually exclusive with ``intersects_bin_size``. Defaults to
        ``None`` (the generic binned equi-join path). Hard-error projection
        shapes raise ``ValueError`` at transpile time; see the performance
        guide for the full enumeration.
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
    if dialect is not None and dialect != "duckdb":
        raise ValueError(f"Unknown dialect: {dialect!r}. Supported: 'duckdb' or None.")
    if dialect == "duckdb" and intersects_bin_size is not None:
        raise ValueError(
            "intersects_bin_size has no effect with dialect='duckdb'; "
            "the DuckDB dialect uses an IEJoin per-partition plan instead "
            "of the binned equi-join. Pass one or the other, not both."
        )

    tables_container = _build_tables(tables)

    with _reraise_as_value_error("Parse error", query=giql):
        ast = parse_one(giql, dialect=GIQLDialect)

    # Falls back to the binned plan for unsupported shapes — see
    # IntersectsDuckDBIEJoinTransformer.transform_to_sql for the complete
    # fallback set.
    if dialect == "duckdb":
        duckdb_transformer = IntersectsDuckDBIEJoinTransformer(tables_container)
        with _reraise_as_value_error("Transformation error"):
            duckdb_sql = duckdb_transformer.transform_to_sql(ast)
        if duckdb_sql is not None:
            return duckdb_sql

    intersects_transformer = IntersectsBinnedJoinTransformer(
        tables_container,
        bin_size=intersects_bin_size,
    )
    merge_transformer = MergeTransformer(tables_container)
    cluster_transformer = ClusterTransformer(tables_container)
    generator = BaseGIQLGenerator(tables=tables_container)

    with _reraise_as_value_error("Transformation error"):
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
    # __giql_canon_* wrapper CTEs for non-canonical interval operands of
    # opted-in operators (GIQL_CANONICALIZE). No operator opts in yet, so this
    # is a strict no-op until the per-operator port issues (#122, #123) land.
    with _reraise_as_value_error("Canonicalization error"):
        ast = canonicalize_coordinates(ast)

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
