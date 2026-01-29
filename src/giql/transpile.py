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


def _build_tables(tables: list[str] | dict[str, Table] | None) -> Tables:
    """Build a Tables container from table specifications.

    Parameters
    ----------
    tables : list[str] | dict[str, Table] | None
        Table specifications. Strings use default column mappings.
        Dict maps table names to Table configurations.

    Returns
    -------
    Tables
        Container with all tables registered.
    """
    container = Tables()

    if tables is None:
        return container

    if isinstance(tables, dict):
        for name, table in tables.items():
            container.register(name, table)
    else:
        for name in tables:
            container.register(name, Table())

    return container


def transpile(
    giql: str,
    tables: list[str] | dict[str, Table] | None = None,
) -> str:
    """Transpile a GIQL query to SQL.

    Parses the GIQL syntax and converts it to standard SQL-92 compatible
    output (uses LATERAL joins where needed for operations like NEAREST).

    Parameters
    ----------
    giql : str
        The GIQL query string containing genomic extensions like
        INTERSECTS, CONTAINS, WITHIN, CLUSTER, MERGE, or NEAREST.
    tables : list[str] | dict[str, Table] | None
        Table configurations. A list of strings uses default column mappings
        (chromosome, start_pos, end_pos, strand). A dict maps table names
        to Table objects for custom column name mappings.

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
            tables={
                "peaks": Table(
                    genomic_col="interval",
                    chrom_col="chrom",
                    start_col="start",
                    end_col="end",
                )
            }
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
