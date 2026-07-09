"""Shared DuckDB per-chromosome dynamic-SQL scaffolding (#217).

Several DuckDB operator overrides reach the fast ``IE_JOIN`` range-join operator
by partitioning a range join per chromosome, which removes the low-cardinality
``chrom`` equality that would otherwise force a ``HASH_JOIN``. Because GIQL does
not know a table's chromosomes at transpile time, the per-chromosome
``UNION ALL`` is assembled at runtime: a ``string_agg`` over the distinct
chromosomes builds the dynamic SQL into a session variable, which the outer query
reads back through ``query(getvariable(...))``.

This module centralizes the target-agnostic scaffolding — the chromosome
string-literal interpolation, a collision-resistant variable name, the
``SET VARIABLE`` builder, and the dynamic-relation reference. The
operator-specific parts (the per-chromosome SQL template, the chromosome
partition source, the empty-schema fallback, and the outer projection) stay with
each expander. First extracted from
:mod:`giql.expanders.intersects_duckdb`; reused by DISJOIN (#216) and NEAREST.
"""

from __future__ import annotations

from uuid import uuid4

#: A DuckDB expression that renders the ``string_agg`` partition's ``chrom``
#: column as a single-quoted, injection-safe SQL string literal, for
#: interpolation into a per-chromosome ``WHERE chrom = '<c>'`` filter.
CHROM_LITERAL = "'''' || replace(chrom, '''', '''''') || ''''"


def sql_escape(text: str) -> str:
    """Return *text* with single quotes doubled for safe SQL string-literal use."""
    return text.replace("'", "''")


def new_variable_name(prefix: str) -> str:
    """Return a collision-resistant DuckDB session-variable name.

    The full uuid4 hex (128 bits) keeps the name unique even across many
    ``transpile()`` calls interleaved in one DuckDB session -- session variables
    are global session state, so a token collision would silently rebind a
    wrapper query to a different dynamic relation.
    """
    return f"{prefix}_{uuid4().hex}"


def set_variable_statement(
    var_name: str,
    per_chrom_template: str,
    chrom_partition: str,
    empty_schema: str,
) -> str:
    """Build the ``SET VARIABLE`` statement that assembles the per-chrom dynamic SQL.

    *per_chrom_template* is a SQL string-concat expression that ``string_agg``
    evaluates once per row of *chrom_partition* (a ``SELECT DISTINCT <chrom>``
    source), producing one per-chromosome ``SELECT`` joined by ``UNION ALL``.
    When no chromosome passes the partition the aggregate is ``NULL``, so it
    falls back to *empty_schema* -- a ``... WHERE FALSE`` query that lets DuckDB
    resolve every output column to its real declared type.
    """
    return (
        f"SET VARIABLE {var_name} = COALESCE((\n"
        f"  SELECT string_agg(\n"
        f"    {per_chrom_template},\n"
        f"    ' UNION ALL '\n"
        f"  )\n"
        f"  FROM ({chrom_partition})\n"
        f"), '{sql_escape(empty_schema)}')"
    )


def dynamic_relation(var_name: str, alias: str) -> str:
    """Return the ``query(getvariable('<var>')) AS <alias>`` relation reference."""
    return f"query(getvariable('{var_name}')) AS {alias}"
