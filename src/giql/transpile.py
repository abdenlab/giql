"""Transpile GIQL queries to SQL.

This module provides the main entry point for transpiling GIQL queries
to standard SQL.
"""

from sqlglot import parse_one

from giql.dialect import GIQLDialect
from giql.generators import BaseGIQLGenerator
from giql.table import Table
from giql.table import Tables
from giql.transformer import ClusterTransformer
from giql.transformer import MergeTransformer


def _build_tables(tables: list[str | Table] | None) -> Tables:
    """Build a Tables container from table specifications.

    Parameters
    ----------
    tables : list[str | Table] | None
        Table specifications. Strings use default column mappings.
        Table objects provide custom column mappings.

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


def transpile(
    giql: str,
    tables: list[str | Table] | None = None,
) -> str:
    """Transpile a GIQL query to SQL.

    Parses the GIQL syntax and converts it to standard SQL-92 compatible
    output (uses LATERAL joins where needed for operations like NEAREST).

    Parameters
    ----------
    giql : str
        The GIQL query string containing genomic extensions like
        INTERSECTS, CONTAINS, WITHIN, CLUSTER, MERGE, or NEAREST.
    tables : list[str | Table] | None
        Table configurations. Strings use default column mappings
        (chrom, start, end, strand). Table objects provide custom
        column name mappings.

    Returns
    -------
    str
        The transpiled SQL query.

    Raises
    ------
    ValueError
        If the query cannot be parsed or transpiled.

    Examples
    --------
    Basic usage with default column mappings::

        sql = transpile(
            "SELECT * FROM peaks WHERE interval INTERSECTS 'chr1:1000-2000'",
            tables=["peaks"]
        )

    Custom table configuration::

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
            ]
        )
    """
    # Build tables container
    tables_container = _build_tables(tables)

    # Initialize transformers with table configurations
    merge_transformer = MergeTransformer(tables_container)
    cluster_transformer = ClusterTransformer(tables_container)

    # Initialize generator with table configurations
    generator = BaseGIQLGenerator(tables=tables_container)

    # Parse GIQL query
    try:
        ast = parse_one(giql, dialect=GIQLDialect)
    except Exception as e:
        raise ValueError(f"Parse error: {e}\nQuery: {giql}") from e

    # Apply transformations (MERGE first, then CLUSTER)
    try:
        # MERGE transformation (which may internally use CLUSTER)
        ast = merge_transformer.transform(ast)
        # CLUSTER transformation for any standalone CLUSTER expressions
        ast = cluster_transformer.transform(ast)
    except Exception as e:
        raise ValueError(f"Transformation error: {e}") from e

    # Generate SQL
    try:
        sql = generator.generate(ast)
    except Exception as e:
        raise ValueError(f"Transpilation error: {e}") from e

    return sql
