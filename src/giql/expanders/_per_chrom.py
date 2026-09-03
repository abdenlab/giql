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

import hashlib

#: A DuckDB expression that renders the ``string_agg`` partition's ``chrom``
#: column as a single-quoted, injection-safe SQL string literal, for
#: interpolation into a per-chromosome ``WHERE chrom = '<c>'`` filter.
CHROM_LITERAL = "'''' || replace(chrom, '''', '''''') || ''''"

#: Partition-count ceiling above which the per-chromosome rewrite degrades to a
#: single non-partitioned query. The rewrite emits one ``UNION ALL`` branch per
#: distinct chromosome, so its cost scales with that cardinality while the plain
#: predicate's does not: the partitioned form wins by orders of magnitude on a
#: primary assembly and loses by orders of magnitude on a scaffold-level one.
#: The ceiling sits above every standard reference assembly's primary-chromosome
#: count and below the measured crossover (#95).
MAX_PER_CHROM_PARTITIONS = 256


def sql_escape(text: str) -> str:
    """Return *text* with single quotes doubled for safe SQL string-literal use."""
    return text.replace("'", "''")


def _content_address(value_sql: str) -> str:
    """Return the session-variable name that addresses *value_sql*.

    The name is a digest of the variable's own rendered value, so re-transpiling
    the same query shape rebinds the same name to a byte-identical value instead
    of minting a fresh one. That bounds a session's variable set by the number of
    distinct query shapes rather than the number of calls -- which matters
    because session variables live for the connection and the emitted script
    cannot release them: the final statement has to be the ``SELECT`` for DuckDB
    to return its result.

    Addressing the *rendered value* rather than the inputs that built it is what
    makes the property hold by construction. Hashing an argument list would
    silently lose the property the day an argument is added and not added to the
    digest, and nothing in the public surface can vary those arguments
    independently enough for a test to catch it.

    Uniqueness is what the previous ``uuid4`` naming bought (#169), and 128 bits
    of SHA-256 preserves it: two different values still get different names, so
    outputs from many ``transpile()`` calls stay safe to interleave in one
    session. Two *identical* values now share a name, where the rebind is a no-op.
    """
    return f"__giql_iejoin_{hashlib.sha256(value_sql.encode()).hexdigest()[:32]}"


def set_variable_statement(
    per_chrom_template: str,
    chrom_partition: str,
    empty_schema: str,
    unpartitioned_query: str,
    *,
    tail_branch: str | None = None,
    max_partitions: int = MAX_PER_CHROM_PARTITIONS,
) -> tuple[str, str]:
    """Return ``(variable name, SET VARIABLE statement)`` for the per-chrom dynamic SQL.

    The name is minted here rather than by the caller, and derived from the
    rendered value rather than from these parameters, so it cannot drift from the
    content it addresses.

    *per_chrom_template* is a SQL string-concat expression that ``string_agg``
    evaluates once per row of *chrom_partition* (a ``SELECT DISTINCT <chrom>``
    source), producing one per-chromosome ``SELECT`` joined by ``UNION ALL``.
    When no chromosome passes the partition the aggregate is ``NULL``, so it
    falls back to *empty_schema* -- a ``... WHERE FALSE`` query that lets DuckDB
    resolve every output column to its real declared type.

    *tail_branch*, when given, is a non-partitioned branch always unioned after
    the per-chromosome ones, and emitted alone when no chromosome passes the
    partition. It supersedes *empty_schema*, which a tail branch makes
    unreachable: the tail is itself a real relation over the source table, so it
    resolves every output column's declared type on its own.

    *unpartitioned_query* is the same result expressed as one query with the
    chromosome equality inlined rather than partitioned out. It is selected at
    execution time when the partition holds more than *max_partitions* rows,
    because branch count -- not row count -- is what makes the partitioned form
    expensive, and the partition's cardinality is a data property no
    transpile-time check can see (#95).
    """
    partition_count = f"(SELECT COUNT(*) FROM ({chrom_partition}))"
    aggregated = (
        f"(\n"
        f"    SELECT string_agg(\n"
        f"      {per_chrom_template},\n"
        f"      ' UNION ALL '\n"
        f"    )\n"
        f"    FROM ({chrom_partition})\n"
        f"  )"
    )
    if tail_branch is None:
        partitioned = f"COALESCE({aggregated}, '{sql_escape(empty_schema)}')"
    else:
        partitioned = (
            f"COALESCE({aggregated} || ' UNION ALL ', '') || '{sql_escape(tail_branch)}'"
        )
    # Built before the name, and referencing no part of it, so the digest below
    # provably covers everything that determines the variable's value.
    value_sql = (
        f"CASE\n"
        f"  WHEN {partition_count} > {max_partitions}\n"
        f"  THEN '{sql_escape(unpartitioned_query)}'\n"
        f"  ELSE {partitioned}\n"
        f"END"
    )
    var_name = _content_address(value_sql)
    return var_name, f"SET VARIABLE {var_name} = {value_sql}"


def dynamic_relation(var_name: str, alias: str) -> str:
    """Return the ``query(getvariable('<var>')) AS <alias>`` relation reference."""
    return f"query(getvariable('{var_name}')) AS {alias}"
