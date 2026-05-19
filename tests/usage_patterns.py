"""Canonical query-usage patterns for GIQL operators.

A shared catalogue of the ways an operator can appear *within* a SQL query --
standalone, filtered, joined, wrapped in a CTE, and so on. The functional suite
(``test_usage_patterns.py``) parametrises over this catalogue, executing every
pattern against a matrix of database engines and snapshotting the results.

Each operator class has a descriptor -- an :class:`OperatorUsage` subclass --
carrying its query templates, its ``str.format`` slot values, the fixture data
its queries run against, and its canonical pattern list. :func:`render`,
:func:`templated_patterns`, and :func:`rendered_cases` dispatch through the
descriptor, so a new operator class is a new subclass with no edits to those
helpers.

Only the table-function class (``DISJOIN``) is populated. ``AggregateUsage``,
``PredicateUsage``, and ``ScalarUsage`` are defined but unpopulated extension
points.

The DISJOIN catalogue is crossed against a second dimension: the
:data:`DISJOIN_PROFILES` collection of named semantic profiles, one per
DISJOIN edge case (point intervals, adjacent intervals, custom column names,
non-canonical encodings, and so on). The functional suite parametrises every
pattern over every profile, snapshotting one manifest per profile.

:func:`rendered_cases` tags each case with ``@pytest.mark.usage`` and, for a
pattern GIQL cannot transpile, a strict ``xfail``.
"""

from __future__ import annotations

import abc
from collections.abc import Mapping
from dataclasses import dataclass
from enum import Enum
from types import MappingProxyType
from typing import ClassVar

import pytest

from giql.table import Table


class OperatorType(Enum):
    """The structural class of a GIQL operator -- the query position it occupies."""

    TABLE_FUNCTION = "table_function"
    """A FROM-clause table function yielding rows (DISJOIN, NEAREST)."""

    AGGREGATE = "aggregate"
    """A SELECT-list aggregate or window operator (MERGE, CLUSTER)."""

    PREDICATE = "predicate"
    """A boolean predicate in WHERE / ON / HAVING (INTERSECTS, CONTAINS, WITHIN)."""

    SCALAR = "scalar"
    """A scalar function usable in any value expression (DISTANCE)."""


class UsagePattern(Enum):
    """A canonical way of placing or combining an operator within a query."""

    STANDALONE = "standalone"
    """The operator is the sole FROM source: ``SELECT * FROM OP(...)``."""

    ALIASED_PROJECTION = "aliased_projection"
    """The result is aliased and specific output columns are projected."""

    OUTER_WHERE = "outer_where"
    """An outer ``WHERE`` filters the result on a real predicate."""

    OUTER_ORDER_LIMIT = "outer_order_limit"
    """An outer ``ORDER BY`` / ``LIMIT`` is applied to the result."""

    OUTER_GROUP_BY = "outer_group_by"
    """An outer ``GROUP BY`` aggregates the result by an output column."""

    OUTER_HAVING = "outer_having"
    """An outer ``GROUP BY`` ... ``HAVING`` filters the grouped result."""

    OUTER_COUNT = "outer_count"
    """An outer ``SELECT COUNT(*)`` collapses the result to a single row."""

    AGGREGATED_INTERVAL = "aggregated_interval"
    """An aggregate over the operator's derived interval columns (coverage)."""

    OUTER_DISTINCT = "outer_distinct"
    """``SELECT DISTINCT`` is applied over the result."""

    OUTER_SPATIAL_WHERE = "outer_spatial_where"
    """An outer interval-overlap predicate filters the result to a region."""

    WINDOWED_PROJECTION = "windowed_projection"
    """A window function is computed over the operator's result."""

    WRAPPING_CTE = "wrapping_cte"
    """The operator is wrapped in a CTE and consumed downstream."""

    ORDERED_CTE = "ordered_cte"
    """The operator is wrapped in a CTE with an inner ``ORDER BY``, consumed outside."""

    DERIVED_SUBQUERY = "derived_subquery"
    """The operator's result is nested as a derived-table subquery."""

    NESTED_AGGREGATE = "nested_aggregate"
    """An aggregate consumes a derived-table subquery of the result."""

    JOINED = "joined"
    """The result is interval-overlap-joined to another relation."""

    OUTER_LEFT_JOIN = "outer_left_join"
    """The result is ``LEFT JOIN``-ed to another relation."""

    SELF_JOIN_OVERLAP = "self_join_overlap"
    """The result is interval-overlap-joined back to the operator's base table."""

    SET_UNION = "set_union"
    """The result is combined with another relation via ``UNION``."""

    SET_DIFFERENCE = "set_difference"
    """The result is combined with another relation via ``EXCEPT``."""

    SUBQUERY_INPUT = "subquery_input"
    """The operator's target is a derived-table subquery."""

    CTE_INPUT = "cte_input"
    """The operator's target is a CTE defined in the outer query."""

    SUBQUERY_REFERENCE = "subquery_reference"
    """The operator's ``reference`` argument is a subquery."""

    CHAINED_INPUT = "chained_input"
    """The operator's target is another operator's result."""

    CHAINED_REFERENCE = "chained_reference"
    """The operator's ``reference`` argument is another operator's result."""

    CHAINED_COOCCURRING = "chained_cooccurring"
    """Two GIQL operator calls co-occur in one query."""


# Reasons GIQL cannot transpile certain patterns -- attached as xfail reasons so
# the catalogue enumerates the expected use while documenting the limitation.
_TARGET_MUST_BE_TABLE = (
    "GIQL requires a table-function target to be a registered base table; "
    "a subquery or CTE in target position is rejected during target resolution"
)
_NESTED_OPERATOR_TARGET = (
    "GIQL does not recognise a nested operator call in a table-function's "
    "target position; the inner call is emitted as a mangled identifier"
)

# SQL templates for the table-function class. Brace slots are filled from a
# TableFunctionUsage by render(). Aliases use a __ prefix to avoid colliding
# with operator-generated identifiers.
_TABLE_FUNCTION_TEMPLATES: dict[UsagePattern, str] = {
    UsagePattern.STANDALONE: "SELECT * FROM {op}",
    UsagePattern.ALIASED_PROJECTION: (
        "SELECT __d.{out_start} AS s, __d.{out_end} AS e FROM {op} AS __d"
    ),
    UsagePattern.OUTER_WHERE: (
        "SELECT * FROM {op} WHERE {out_end} - {out_start} >= 50"
    ),
    UsagePattern.OUTER_ORDER_LIMIT: (
        "SELECT * FROM {op} ORDER BY {out_start} LIMIT 1"
    ),
    UsagePattern.OUTER_GROUP_BY: (
        "SELECT {out_chrom}, COUNT(*) AS n FROM {op} GROUP BY {out_chrom}"
    ),
    UsagePattern.OUTER_HAVING: (
        "SELECT {out_chrom}, COUNT(*) AS n FROM {op} "
        "GROUP BY {out_chrom} HAVING COUNT(*) > 1"
    ),
    # A scalar COUNT(*) over the whole result: collapses to one deterministic
    # row regardless of the result's own row order.
    UsagePattern.OUTER_COUNT: "SELECT COUNT(*) AS n FROM {op}",
    UsagePattern.AGGREGATED_INTERVAL: (
        "SELECT {out_chrom}, SUM({out_end} - {out_start}) AS total_bp "
        "FROM {op} GROUP BY {out_chrom}"
    ),
    UsagePattern.OUTER_DISTINCT: (
        "SELECT DISTINCT {out_chrom}, {out_start}, {out_end} FROM {op}"
    ),
    # A region filter over the operator's output. GIQL's INTERSECTS does not
    # work here: applied to a table-function alias it silently resolves to the
    # alias's default chrom/start/end (the passed-through parent interval), not
    # the disjoin_* segment -- so an explicit half-open overlap predicate is
    # the correct spatial filter on the operator's own output.
    UsagePattern.OUTER_SPATIAL_WHERE: (
        "SELECT * FROM {op} "
        "WHERE {out_chrom} = 'chr1' AND {out_start} < 50 AND {out_end} > 20"
    ),
    # COUNT(*) OVER, not ROW_NUMBER(): a window aggregate is deterministic,
    # whereas ROW_NUMBER over a non-total ORDER BY (duplicate segments share a
    # start) would assign ranks non-deterministically and flake the snapshot.
    UsagePattern.WINDOWED_PROJECTION: (
        "SELECT *, COUNT(*) OVER (PARTITION BY {out_chrom}) AS chrom_n "
        "FROM {op}"
    ),
    UsagePattern.WRAPPING_CTE: (
        "WITH __u AS (SELECT * FROM {op}) SELECT * FROM __u"
    ),
    # An inner ORDER BY inside the CTE. The outer query re-selects unordered;
    # the functional suite re-sorts every result by repr, so the snapshot is
    # stable even though the inner ORDER BY is not total.
    UsagePattern.ORDERED_CTE: (
        "WITH __o AS (SELECT * FROM {op} ORDER BY {out_start}) "
        "SELECT * FROM __o"
    ),
    UsagePattern.DERIVED_SUBQUERY: (
        "SELECT * FROM (SELECT * FROM {op}) AS __s"
    ),
    # An aggregate consuming a derived-table subquery of the result. COUNT(*)
    # and SUM are order-insensitive, so the single output row is deterministic.
    UsagePattern.NESTED_AGGREGATE: (
        "SELECT COUNT(*) AS n, SUM({out_end} - {out_start}) AS total_bp "
        "FROM (SELECT * FROM {op}) AS __s"
    ),
    UsagePattern.JOINED: (
        "SELECT __d.{out_chrom}, __d.{out_start}, __d.{out_end}, __o.name "
        "FROM {op} AS __d JOIN {join_table} AS __o "
        "ON __o.chrom = __d.{out_chrom} "
        "AND __o.start < __d.{out_end} AND __o.end > __d.{out_start}"
    ),
    # A LEFT JOIN: every result sub-interval survives, with NULLs where the
    # join_table has no overlapping row. The functional suite re-sorts rows.
    UsagePattern.OUTER_LEFT_JOIN: (
        "SELECT __d.{out_chrom}, __d.{out_start}, __d.{out_end}, __o.name "
        "FROM {op} AS __d LEFT JOIN {join_table} AS __o "
        "ON __o.chrom = __d.{out_chrom} "
        "AND __o.start < __d.{out_end} AND __o.end > __d.{out_start}"
    ),
    # The result interval-overlap-joined back to its own base target table --
    # a real relation join, not a nested operator call. The join addresses the
    # target's physical genomic columns (custom for a custom-column target);
    # the projected parent start proves a real target row was matched.
    UsagePattern.SELF_JOIN_OVERLAP: (
        "SELECT __d.{out_chrom}, __d.{out_start}, __d.{out_end}, "
        "__t.{tgt_start} AS parent_start "
        "FROM {op} AS __d JOIN {target} AS __t "
        "ON __t.{tgt_chrom} = __d.{out_chrom} "
        "AND __t.{tgt_start} < __d.{out_end} "
        "AND __t.{tgt_end} > __d.{out_start}"
    ),
    UsagePattern.SET_UNION: (
        "SELECT {out_chrom} AS c FROM {op} "
        "UNION ALL SELECT chrom AS c FROM {join_table}"
    ),
    UsagePattern.SET_DIFFERENCE: (
        "SELECT {out_chrom} AS c FROM {op} "
        "EXCEPT SELECT chrom AS c FROM {join_table}"
    ),
    UsagePattern.SUBQUERY_INPUT: (
        "SELECT * FROM {name}((SELECT * FROM {target}))"
    ),
    UsagePattern.CTE_INPUT: (
        "WITH __in AS (SELECT * FROM {target}) SELECT * FROM {name}(__in)"
    ),
    # The reference subquery projects the operand's genomic columns aliased to
    # the canonical chrom/start/end names. GIQL's subquery-reference resolution
    # assumes a canonical layout, so a custom-column operand must be aliased
    # explicitly here rather than passed through with SELECT *.
    UsagePattern.SUBQUERY_REFERENCE: (
        "SELECT * FROM {name}({target}, "
        "reference := (SELECT {ref_select} FROM {ref_operand}))"
    ),
    UsagePattern.CHAINED_INPUT: (
        "SELECT * FROM {name}((SELECT * FROM {name}({target})))"
    ),
    # The reference is another DISJOIN call; its output exposes the
    # convention-fixed disjoin_* columns, aliased here to canonical names so
    # the outer DISJOIN's subquery-reference resolution can resolve them.
    UsagePattern.CHAINED_REFERENCE: (
        "SELECT * FROM {name}({target}, "
        "reference := (SELECT {out_chrom} AS chrom, {out_start} AS \"start\", "
        "{out_end} AS \"end\" FROM {name}({ref_operand})))"
    ),
    UsagePattern.CHAINED_COOCCURRING: (
        "SELECT __a.{out_start} FROM {op} AS __a JOIN {op} AS __b "
        "ON __a.{out_chrom} = __b.{out_chrom} "
        "AND __a.{out_start} < __b.{out_end} AND __a.{out_end} > __b.{out_start}"
    ),
}


# The canonical physical schema for an interval table: the standard genomic
# layout the original DISJOIN fixtures and the bedtools-compatible defaults
# assume. A TableFixture with no explicit `columns` uses this.
CANONICAL_COLUMNS: tuple[tuple[str, str], ...] = (
    ("chrom", "VARCHAR"),
    ("start", "INTEGER"),
    ("end", "INTEGER"),
    ("name", "VARCHAR"),
)


@dataclass(frozen=True)
class TableFixture:
    """Fixture rows for one table the functional suite loads before executing.

    ``columns`` is the physical schema as ``(name, sql_type)`` pairs; it
    defaults to :data:`CANONICAL_COLUMNS` (``chrom`` VARCHAR, ``start``
    INTEGER, ``end`` INTEGER, ``name`` VARCHAR). Each row tuple must match the
    declared schema width and column order.

    ``table`` is an optional :class:`giql.table.Table` configuration. When
    present it carries custom column names and/or a non-canonical coordinate
    system, and the functional suite passes it to ``transpile`` so GIQL
    resolves the fixture's physical layout. When ``None`` the table uses
    GIQL's default canonical mapping.
    """

    name: str
    rows: tuple[tuple, ...]
    columns: tuple[tuple[str, str], ...] = CANONICAL_COLUMNS
    table: Table | None = None


class OperatorUsage(abc.ABC):
    """Descriptor base for one operator's usage-pattern slot values and fixtures.

    One subclass exists per :class:`OperatorType`. A subclass is a frozen
    dataclass providing, as attributes/fields, ``name`` (the operator name),
    ``operator_type``, ``canonical_patterns`` (the patterns its class is
    expected to cover), ``xfails`` (patterns GIQL cannot transpile, mapped to a
    reason), and ``fixtures`` (the table data its queries run against). It
    implements :meth:`tables`, :meth:`template_table`, and :meth:`slots`.

    :func:`render`, :func:`templated_patterns`, and :func:`rendered_cases`
    dispatch through those members, so adding an operator class needs no edits
    to the helpers themselves.
    """

    operator_type: ClassVar[OperatorType]
    canonical_patterns: ClassVar[tuple[UsagePattern, ...]] = ()
    xfails: ClassVar[Mapping[UsagePattern, str]] = MappingProxyType({})

    @property
    @abc.abstractmethod
    def tables(self) -> tuple[str | Table, ...]:
        """Return the table configs to register when transpiling this operator."""

    @abc.abstractmethod
    def template_table(self) -> dict[UsagePattern, str]:
        """Return the ``{pattern: SQL template}`` table for this operator's class."""

    @abc.abstractmethod
    def slots(self) -> dict[str, str]:
        """Return the ``str.format`` keyword values for this operator's templates."""


@dataclass(frozen=True)
class TableFunctionUsage(OperatorUsage):
    """Usage descriptor for a FROM-clause table-function operator (e.g. DISJOIN).

    ``target`` must be a registered base table -- GIQL rejects a subquery or CTE
    in target position. ``reference`` is the optional two-table-mode reference
    relation; ``join_table`` is a second interval table used by the ``JOINED``
    and set-operation patterns. ``fixtures`` supplies the rows of every table
    the patterns reference.

    ``profile_id`` names the semantic edge case the descriptor exercises; it is
    the manifest-key prefix that keeps every profile's snapshot rows distinct.
    """

    name: str
    target: str
    out_chrom: str
    out_start: str
    out_end: str
    fixtures: tuple[TableFixture, ...]
    reference: str | None = None
    join_table: str = "annotations"
    profile_id: str = "CANONICAL"

    operator_type: ClassVar[OperatorType] = OperatorType.TABLE_FUNCTION
    canonical_patterns: ClassVar[tuple[UsagePattern, ...]] = (
        UsagePattern.STANDALONE,
        UsagePattern.ALIASED_PROJECTION,
        UsagePattern.OUTER_WHERE,
        UsagePattern.OUTER_ORDER_LIMIT,
        UsagePattern.OUTER_GROUP_BY,
        UsagePattern.OUTER_HAVING,
        UsagePattern.OUTER_COUNT,
        UsagePattern.AGGREGATED_INTERVAL,
        UsagePattern.OUTER_DISTINCT,
        UsagePattern.OUTER_SPATIAL_WHERE,
        UsagePattern.WINDOWED_PROJECTION,
        UsagePattern.WRAPPING_CTE,
        UsagePattern.ORDERED_CTE,
        UsagePattern.DERIVED_SUBQUERY,
        UsagePattern.NESTED_AGGREGATE,
        UsagePattern.JOINED,
        UsagePattern.OUTER_LEFT_JOIN,
        UsagePattern.SELF_JOIN_OVERLAP,
        UsagePattern.SET_UNION,
        UsagePattern.SET_DIFFERENCE,
        UsagePattern.SUBQUERY_INPUT,
        UsagePattern.CTE_INPUT,
        UsagePattern.SUBQUERY_REFERENCE,
        UsagePattern.CHAINED_INPUT,
        UsagePattern.CHAINED_REFERENCE,
        UsagePattern.CHAINED_COOCCURRING,
    )
    xfails: ClassVar[Mapping[UsagePattern, str]] = MappingProxyType(
        {
            UsagePattern.SUBQUERY_INPUT: _TARGET_MUST_BE_TABLE,
            UsagePattern.CTE_INPUT: _TARGET_MUST_BE_TABLE,
            UsagePattern.CHAINED_INPUT: _NESTED_OPERATOR_TARGET,
        }
    )

    @property
    def op(self) -> str:
        """Return the operator call as it appears in a FROM clause."""
        if self.reference is None:
            return f"{self.name}({self.target})"
        return f"{self.name}({self.target}, reference := {self.reference})"

    @property
    def ref_operand(self) -> str:
        """Return the relation a reference-bearing template should reference."""
        return self.target if self.reference is None else self.reference

    def _genomic_columns(self, table_name: str) -> tuple[str, str, str]:
        """Return the ``(chrom, start, end)`` physical column names of a fixture.

        A fixture carrying a :class:`giql.table.Table` config reports that
        config's custom column names; otherwise the canonical ``chrom`` /
        ``start`` / ``end`` layout is assumed. Each name is double-quoted so a
        template can interpolate it safely even when the physical column is a
        SQL reserved word (the canonical layout's ``start`` and ``end`` are).
        """
        for fixture in self.fixtures:
            if fixture.name == table_name and fixture.table is not None:
                config = fixture.table
                names = (config.chrom_col, config.start_col, config.end_col)
                break
        else:
            names = ("chrom", "start", "end")
        return tuple(f'"{name}"' for name in names)

    @property
    def tables(self) -> tuple[str | Table, ...]:
        """Return the table configs this operator's patterns reference.

        A fixture carrying a custom :class:`giql.table.Table` config yields the
        ``Table`` object so GIQL resolves its physical columns and coordinate
        system; a fixture with the canonical layout yields its bare name.
        """
        return tuple(
            fixture.table if fixture.table is not None else fixture.name
            for fixture in self.fixtures
        )

    def template_table(self) -> dict[UsagePattern, str]:
        """Return the table-function template table."""
        return _TABLE_FUNCTION_TEMPLATES

    def slots(self) -> dict[str, str]:
        """Return the template slot values for this table-function usage.

        ``tgt_chrom`` / ``tgt_start`` / ``tgt_end`` are the target table's
        physical genomic column names (custom when the target fixture carries
        a :class:`giql.table.Table` config), so a template that joins the
        result back to the base target table addresses real columns.
        ``ref_select`` is a canonical-column projection of the reference
        operand: it aliases the operand's physical columns to ``chrom`` /
        ``start`` / ``end``, the column names GIQL's subquery / CTE reference
        resolution requires.
        """
        tgt_chrom, tgt_start, tgt_end = self._genomic_columns(self.target)
        ref_chrom, ref_start, ref_end = self._genomic_columns(self.ref_operand)
        ref_select = (
            f'{ref_chrom} AS chrom, {ref_start} AS "start", '
            f'{ref_end} AS "end"'
        )
        return {
            "name": self.name,
            "op": self.op,
            "target": self.target,
            "ref_operand": self.ref_operand,
            "ref_select": ref_select,
            "out_chrom": self.out_chrom,
            "out_start": self.out_start,
            "out_end": self.out_end,
            "tgt_chrom": tgt_chrom,
            "tgt_start": tgt_start,
            "tgt_end": tgt_end,
            "join_table": self.join_table,
        }


class AggregateUsage(OperatorUsage):
    """Extension point -- SELECT-list aggregate operators (MERGE, CLUSTER).

    Not yet populated. To populate: make this a frozen dataclass with slots for
    the aggregate call, the source table, an optional result alias (the window
    form), and operator arguments; set ``name``, ``canonical_patterns`` (an
    aggregate-specific pattern set), and ``fixtures``; implement
    :meth:`tables`, :meth:`template_table`, and :meth:`slots`. See
    :class:`TableFunctionUsage` -- keep ``operator_type`` and ``xfails``
    annotated ``ClassVar`` so the dataclass does not treat them as init fields.
    """

    operator_type: ClassVar[OperatorType] = OperatorType.AGGREGATE


class PredicateUsage(OperatorUsage):
    """Extension point -- WHERE/ON/HAVING boolean predicates (INTERSECTS, etc.).

    Not yet populated. A predicate's patterns (in a ``WHERE``, in a join ``ON``,
    negated, AND/OR-combined, ANY/ALL-quantified, literal-range vs. column
    operand) are distinct from a table function's; populate this with a
    predicate-specific ``canonical_patterns`` set and template table. See
    :class:`TableFunctionUsage` for the descriptor shape and the
    ``ClassVar`` caveat.
    """

    operator_type: ClassVar[OperatorType] = OperatorType.PREDICATE


class ScalarUsage(OperatorUsage):
    """Extension point -- scalar functions usable in value expressions (DISTANCE).

    Not yet populated. A scalar's patterns (in the SELECT list, in a ``WHERE``
    comparison, in ``ORDER BY``, wrapped in an aggregate) are distinct from a
    table function's; populate this with a scalar-specific ``canonical_patterns``
    set and template table. See :class:`TableFunctionUsage` for the descriptor
    shape and the ``ClassVar`` caveat.
    """

    operator_type: ClassVar[OperatorType] = OperatorType.SCALAR


def templated_patterns(usage: OperatorUsage) -> tuple[UsagePattern, ...]:
    """Return the canonical patterns :func:`render` can produce for ``usage``.

    Patterns are returned in canonical order, restricted to those the operator's
    class has a template for.
    """
    templates = usage.template_table()
    return tuple(
        pattern for pattern in usage.canonical_patterns if pattern in templates
    )


def render(pattern: UsagePattern, usage: OperatorUsage) -> str:
    """Render the GIQL query that exercises ``pattern`` for ``usage``.

    :param pattern:
        The usage pattern to render.
    :param usage:
        The operator's usage descriptor.
    :return:
        A complete GIQL query string.
    :raises ValueError:
        If the operator's class has no template for ``pattern`` -- for example
        an unpopulated operator class.
    """
    templates = usage.template_table()
    if pattern not in templates:
        raise ValueError(
            f"{usage.name}: usage pattern {pattern.name} has no template for "
            f"operator class {usage.operator_type.name}"
        )
    return templates[pattern].format(**usage.slots())


def rendered_cases(usage: OperatorUsage) -> list:
    """Return rendered usage cases for ``usage`` as ``pytest.param`` values.

    Each param carries ``(pattern, query)``, a stable id, and a
    ``pytest.mark.usage`` mark; patterns in ``usage.xfails`` additionally carry
    a strict ``pytest.mark.xfail`` with the documented reason. The result is
    ready to splat into ``@pytest.mark.parametrize(("pattern", "query"), ...)``.
    """
    cases: list = []
    for pattern in templated_patterns(usage):
        marks = [pytest.mark.usage(pattern)]
        if pattern in usage.xfails:
            marks.append(
                pytest.mark.xfail(reason=usage.xfails[pattern], strict=True)
            )
        cases.append(
            pytest.param(
                pattern,
                render(pattern, usage),
                id=f"{usage.name}-{pattern.name}",
                marks=marks,
            )
        )
    return cases


# Canonical fixture tables.
#
# The intervals are chosen so each pattern produces a non-trivial,
# discriminating result: the self-mode split, the reference-mode coverage
# filter, the overlap join, and the length/region filters all drop or
# transform some rows. The CANONICAL profile reuses these unchanged so its
# snapshot rows stay semantically identical to the originally committed
# manifests.
_FEATURES = TableFixture(
    "features",
    (("chr1", 0, 100, "a"), ("chr1", 60, 180, "b"), ("chr2", 10, 40, "c")),
)
_MASK = TableFixture(
    "mask",
    (("chr1", 0, 70, "m1"), ("chr1", 90, 200, "m2"), ("chr2", 0, 100, "m3")),
)
_ANNOTATIONS = TableFixture(
    "annotations",
    (("chr1", 0, 30, "g1"), ("chr2", 0, 1000, "g2")),
)


def _features_like(name: str, rows: tuple[tuple, ...]) -> TableFixture:
    """Build a canonical-layout interval fixture with the given name and rows."""
    return TableFixture(name, rows)


# A small canonical annotations table reused by profiles whose own target /
# reference fixtures already isolate the behaviour under test; it backs the
# JOINED, SELF_JOIN_OVERLAP, and set-operation patterns.
def _annotations(rows: tuple[tuple, ...] = _ANNOTATIONS.rows) -> TableFixture:
    """Build an ``annotations`` fixture (canonical layout) with the given rows."""
    return _features_like("annotations", rows)


def _table_function_profile(
    profile_id: str,
    *,
    features: tuple[tuple, ...],
    annotations: tuple[tuple, ...] = _ANNOTATIONS.rows,
    mask: tuple[tuple, ...] | None = None,
    reference: bool = False,
) -> TableFunctionUsage:
    """Build a canonical-layout DISJOIN profile from raw row tuples.

    When ``reference`` is true the profile is two-table-mode against a ``mask``
    table; otherwise it is self-mode. ``annotations`` backs the join and
    set-operation patterns.
    """
    feat_fixture = _features_like("features", features)
    ann_fixture = _annotations(annotations)
    if reference:
        if mask is None:
            raise ValueError(f"{profile_id}: reference-mode profile needs a mask")
        return TableFunctionUsage(
            name="DISJOIN",
            target="features",
            out_chrom="disjoin_chrom",
            out_start="disjoin_start",
            out_end="disjoin_end",
            fixtures=(feat_fixture, _features_like("mask", mask), ann_fixture),
            reference="mask",
            profile_id=profile_id,
        )
    return TableFunctionUsage(
        name="DISJOIN",
        target="features",
        out_chrom="disjoin_chrom",
        out_start="disjoin_start",
        out_end="disjoin_end",
        fixtures=(feat_fixture, ann_fixture),
        profile_id=profile_id,
    )


# DISJOIN semantic profiles.
#
# Each profile is a TableFunctionUsage exercising one DISJOIN edge case. The
# functional suite parametrises every usage pattern over every profile and
# snapshots one manifest per (profile, mode). CANONICAL reuses the original
# self / reference fixtures so its rows match the first committed manifests.

# CANONICAL -- the original benign data profile, self-mode and reference-mode.
DISJOIN_USAGE = TableFunctionUsage(
    name="DISJOIN",
    target="features",
    out_chrom="disjoin_chrom",
    out_start="disjoin_start",
    out_end="disjoin_end",
    fixtures=(_FEATURES, _ANNOTATIONS),
    profile_id="CANONICAL",
)
DISJOIN_REFERENCE_USAGE = TableFunctionUsage(
    name="DISJOIN",
    target="features",
    out_chrom="disjoin_chrom",
    out_start="disjoin_start",
    out_end="disjoin_end",
    fixtures=(_FEATURES, _MASK, _ANNOTATIONS),
    reference="mask",
    profile_id="CANONICAL",
)

# POINT_INTERVALS -- every target interval is zero-width (start == end). No
# cut produces a positive-length segment, so DISJOIN yields no rows.
_POINT_SELF = _table_function_profile(
    "POINT_INTERVALS",
    features=(("chr1", 10, 10, "p1"), ("chr1", 20, 20, "p2"), ("chr2", 5, 5, "p3")),
)
_POINT_REF = _table_function_profile(
    "POINT_INTERVALS",
    features=(("chr1", 10, 10, "p1"), ("chr1", 20, 20, "p2"), ("chr2", 5, 5, "p3")),
    mask=(("chr1", 0, 100, "m1"), ("chr2", 0, 100, "m2")),
    reference=True,
)

# ADJACENT_INTERVALS -- intervals that touch end-to-start but never overlap.
# Self-mode cuts only at shared boundaries; segments stay whole.
_ADJACENT_SELF = _table_function_profile(
    "ADJACENT_INTERVALS",
    features=(
        ("chr1", 0, 30, "a"),
        ("chr1", 30, 60, "b"),
        ("chr1", 60, 90, "c"),
    ),
)
_ADJACENT_REF = _table_function_profile(
    "ADJACENT_INTERVALS",
    features=(("chr1", 0, 90, "t"),),
    mask=(("chr1", 0, 30, "r1"), ("chr1", 30, 60, "r2"), ("chr1", 60, 90, "r3")),
    reference=True,
)

# NESTED_INTERVALS -- an interval fully contained in another. The inner
# interval's endpoints are interior cuts of the outer one.
_NESTED_SELF = _table_function_profile(
    "NESTED_INTERVALS",
    features=(
        ("chr1", 0, 100, "outer"),
        ("chr1", 30, 60, "inner"),
        ("chr1", 40, 50, "innermost"),
    ),
)
_NESTED_REF = _table_function_profile(
    "NESTED_INTERVALS",
    features=(("chr1", 0, 100, "outer"),),
    mask=(("chr1", 30, 60, "inner"), ("chr1", 40, 50, "innermost")),
    reference=True,
)

# IDENTICAL_INTERVALS -- several rows share the exact same interval. DISJOIN
# partitions cuts per (chrom, start, end), so each duplicate is split
# independently and contributes its own sub-interval rows.
_IDENTICAL_SELF = _table_function_profile(
    "IDENTICAL_INTERVALS",
    features=(
        ("chr1", 0, 50, "x1"),
        ("chr1", 0, 50, "x2"),
        ("chr1", 20, 80, "y"),
    ),
)
_IDENTICAL_REF = _table_function_profile(
    "IDENTICAL_INTERVALS",
    features=(("chr1", 0, 50, "x1"), ("chr1", 0, 50, "x2")),
    mask=(("chr1", 20, 80, "r1"), ("chr1", 20, 80, "r2")),
    reference=True,
)

# DUPLICATE_TARGETS -- duplicate target rows split against a single reference
# breakpoint set; each duplicate parent produces the same sub-intervals.
_DUPLICATE_SELF = _table_function_profile(
    "DUPLICATE_TARGETS",
    features=(
        ("chr1", 0, 100, "d1"),
        ("chr1", 0, 100, "d2"),
        ("chr1", 0, 100, "d3"),
        ("chr1", 40, 60, "cut"),
    ),
)
_DUPLICATE_REF = _table_function_profile(
    "DUPLICATE_TARGETS",
    features=(
        ("chr1", 0, 100, "d1"),
        ("chr1", 0, 100, "d2"),
        ("chr1", 0, 100, "d3"),
    ),
    mask=(("chr1", 40, 60, "r1"),),
    reference=True,
)

# MANY_CHROMOSOMES -- intervals spread across four chromosomes. No sub-interval
# may cross a chromosome boundary; cuts partition per chromosome.
_MANY_SELF = _table_function_profile(
    "MANY_CHROMOSOMES",
    features=(
        ("chr1", 0, 100, "a"),
        ("chr1", 50, 150, "b"),
        ("chr2", 0, 60, "c"),
        ("chr2", 30, 90, "d"),
        ("chr3", 10, 40, "e"),
        ("chrX", 5, 25, "f"),
    ),
    annotations=(
        ("chr1", 0, 1000, "g1"),
        ("chr2", 0, 1000, "g2"),
        ("chr3", 0, 1000, "g3"),
        ("chrX", 0, 1000, "g4"),
    ),
)
_MANY_REF = _table_function_profile(
    "MANY_CHROMOSOMES",
    features=(
        ("chr1", 0, 100, "a"),
        ("chr2", 0, 60, "c"),
        ("chr3", 10, 40, "e"),
        ("chrX", 5, 25, "f"),
    ),
    mask=(
        ("chr1", 50, 150, "r1"),
        ("chr2", 30, 90, "r2"),
        ("chr3", 0, 1000, "r3"),
        ("chrX", 0, 1000, "r4"),
    ),
    annotations=(
        ("chr1", 0, 1000, "g1"),
        ("chr2", 0, 1000, "g2"),
        ("chr3", 0, 1000, "g3"),
        ("chrX", 0, 1000, "g4"),
    ),
    reference=True,
)

# SINGLE_CHROMOSOME -- all intervals on one contig; a dense overlap stack.
_SINGLE_SELF = _table_function_profile(
    "SINGLE_CHROMOSOME",
    features=(
        ("chr1", 0, 80, "a"),
        ("chr1", 20, 100, "b"),
        ("chr1", 40, 120, "c"),
    ),
    annotations=(("chr1", 0, 1000, "g1"),),
)
_SINGLE_REF = _table_function_profile(
    "SINGLE_CHROMOSOME",
    features=(("chr1", 0, 120, "t"),),
    mask=(("chr1", 20, 60, "r1"), ("chr1", 60, 100, "r2")),
    annotations=(("chr1", 0, 1000, "g1"),),
    reference=True,
)

# REFERENCE_GAP -- the reference set leaves a hole inside a target interval.
# DISJOIN cuts the target at the gap's edges, then the coverage filter drops
# the sub-interval lying over the uncovered gap.
_GAP_SELF = _table_function_profile(
    "REFERENCE_GAP",
    features=(
        ("chr1", 0, 40, "a"),
        ("chr1", 60, 100, "b"),
        ("chr1", 0, 100, "spanning"),
    ),
)
_GAP_REF = _table_function_profile(
    "REFERENCE_GAP",
    features=(("chr1", 0, 100, "t"),),
    mask=(("chr1", 0, 40, "r1"), ("chr1", 60, 100, "r2")),
    reference=True,
)

# REFERENCE_SUPERSET -- the reference covers far more than the target. Every
# target sub-interval survives the coverage filter; cuts come from references
# interior to the target only.
_SUPERSET_SELF = _table_function_profile(
    "REFERENCE_SUPERSET",
    features=(
        ("chr1", 20, 80, "a"),
        ("chr1", 0, 200, "wide"),
    ),
    annotations=(("chr1", 0, 1000, "g1"),),
)
_SUPERSET_REF = _table_function_profile(
    "REFERENCE_SUPERSET",
    features=(("chr1", 40, 60, "t"),),
    mask=(("chr1", 0, 1000, "r1"), ("chr1", 0, 500, "r2")),
    annotations=(("chr1", 0, 1000, "g1"),),
    reference=True,
)

# EMPTY_TARGET -- the target table has no rows. DISJOIN yields nothing and
# must not error.
_EMPTY_TARGET_SELF = _table_function_profile(
    "EMPTY_TARGET",
    features=(),
)
_EMPTY_TARGET_REF = _table_function_profile(
    "EMPTY_TARGET",
    features=(),
    mask=(("chr1", 0, 100, "r1"),),
    reference=True,
)

# EMPTY_REFERENCE -- the reference table has no rows, so there are no
# breakpoints and no coverage; DISJOIN yields nothing.
_EMPTY_REFERENCE_SELF = _table_function_profile(
    "EMPTY_REFERENCE",
    features=(),
)
_EMPTY_REFERENCE_REF = _table_function_profile(
    "EMPTY_REFERENCE",
    features=(("chr1", 0, 100, "t"), ("chr1", 50, 150, "u")),
    mask=(),
    reference=True,
)

# CUSTOM_COLUMNS -- the target and reference store intervals under custom
# physical column names (seqid / lo / hi). The result must be identical to
# CANONICAL: the disjoin_* output columns are convention-fixed.
_CUSTOM_SCHEMA: tuple[tuple[str, str], ...] = (
    ("seqid", "VARCHAR"),
    ("lo", "INTEGER"),
    ("hi", "INTEGER"),
    ("label", "VARCHAR"),
)
_CUSTOM_FEATURES = TableFixture(
    "features",
    (("chr1", 0, 100, "a"), ("chr1", 60, 180, "b"), ("chr2", 10, 40, "c")),
    columns=_CUSTOM_SCHEMA,
    table=Table("features", chrom_col="seqid", start_col="lo", end_col="hi"),
)
_CUSTOM_MASK = TableFixture(
    "mask",
    (("chr1", 0, 70, "m1"), ("chr1", 90, 200, "m2"), ("chr2", 0, 100, "m3")),
    columns=_CUSTOM_SCHEMA,
    table=Table("mask", chrom_col="seqid", start_col="lo", end_col="hi"),
)
_CUSTOM_SELF = TableFunctionUsage(
    name="DISJOIN",
    target="features",
    out_chrom="disjoin_chrom",
    out_start="disjoin_start",
    out_end="disjoin_end",
    fixtures=(_CUSTOM_FEATURES, _ANNOTATIONS),
    profile_id="CUSTOM_COLUMNS",
)
_CUSTOM_REF = TableFunctionUsage(
    name="DISJOIN",
    target="features",
    out_chrom="disjoin_chrom",
    out_start="disjoin_start",
    out_end="disjoin_end",
    fixtures=(_CUSTOM_FEATURES, _CUSTOM_MASK, _ANNOTATIONS),
    reference="mask",
    profile_id="CUSTOM_COLUMNS",
)

# ONE_BASED_CLOSED_TARGET -- the target is stored 1-based closed. GIQL must
# canonicalize endpoints; the disjoin_* sub-intervals come back 0-based
# half-open. The 1-based-closed encoding of canonical [0,100) is [1,100],
# of [60,180) is [61,180], of [10,40) is [11,40].
_ONE_BASED_FEATURES = TableFixture(
    "features",
    (("chr1", 1, 100, "a"), ("chr1", 61, 180, "b"), ("chr2", 11, 40, "c")),
    table=Table("features", coordinate_system="1based", interval_type="closed"),
)
_ONE_BASED_SELF = TableFunctionUsage(
    name="DISJOIN",
    target="features",
    out_chrom="disjoin_chrom",
    out_start="disjoin_start",
    out_end="disjoin_end",
    fixtures=(_ONE_BASED_FEATURES, _ANNOTATIONS),
    profile_id="ONE_BASED_CLOSED_TARGET",
)
# 1-based-closed mask encoding of canonical [0,70)->[1,70], [90,200)->[91,200],
# [0,100)->[1,100].
_ONE_BASED_MASK = TableFixture(
    "mask",
    (("chr1", 1, 70, "m1"), ("chr1", 91, 200, "m2"), ("chr2", 1, 100, "m3")),
    table=Table("mask", coordinate_system="1based", interval_type="closed"),
)
_ONE_BASED_REF = TableFunctionUsage(
    name="DISJOIN",
    target="features",
    out_chrom="disjoin_chrom",
    out_start="disjoin_start",
    out_end="disjoin_end",
    fixtures=(_ONE_BASED_FEATURES, _ONE_BASED_MASK, _ANNOTATIONS),
    reference="mask",
    profile_id="ONE_BASED_CLOSED_TARGET",
)

# MIXED_ENCODING -- target and reference use different conventions. The target
# is 0-based half-open; the reference is 1-based closed. Both canonicalize to
# the same 0-based half-open space, so the result matches CANONICAL's.
_MIXED_FEATURES = _FEATURES  # canonical 0-based half-open
_MIXED_MASK = TableFixture(
    "mask",
    (("chr1", 1, 70, "m1"), ("chr1", 91, 200, "m2"), ("chr2", 1, 100, "m3")),
    table=Table("mask", coordinate_system="1based", interval_type="closed"),
)
_MIXED_SELF = TableFunctionUsage(
    name="DISJOIN",
    target="features",
    out_chrom="disjoin_chrom",
    out_start="disjoin_start",
    out_end="disjoin_end",
    fixtures=(_MIXED_FEATURES, _ANNOTATIONS),
    profile_id="MIXED_ENCODING",
)
_MIXED_REF = TableFunctionUsage(
    name="DISJOIN",
    target="features",
    out_chrom="disjoin_chrom",
    out_start="disjoin_start",
    out_end="disjoin_end",
    fixtures=(_MIXED_FEATURES, _MIXED_MASK, _ANNOTATIONS),
    reference="mask",
    profile_id="MIXED_ENCODING",
)


@dataclass(frozen=True)
class DisjoinProfile:
    """A named DISJOIN semantic profile: its self-mode and reference-mode usages.

    Each profile binds one edge-case fixture set to both a self-mode and a
    reference-mode :class:`TableFunctionUsage`, so the functional suite can
    parametrise every usage pattern over the same data in both modes.
    """

    profile_id: str
    self_usage: TableFunctionUsage
    reference_usage: TableFunctionUsage


# The full DISJOIN profile matrix. Profile ids are unique -- they prefix the
# manifest keys, so a collision would silently overwrite snapshot rows.
DISJOIN_PROFILES: tuple[DisjoinProfile, ...] = (
    DisjoinProfile("CANONICAL", DISJOIN_USAGE, DISJOIN_REFERENCE_USAGE),
    DisjoinProfile("POINT_INTERVALS", _POINT_SELF, _POINT_REF),
    DisjoinProfile("ADJACENT_INTERVALS", _ADJACENT_SELF, _ADJACENT_REF),
    DisjoinProfile("NESTED_INTERVALS", _NESTED_SELF, _NESTED_REF),
    DisjoinProfile("IDENTICAL_INTERVALS", _IDENTICAL_SELF, _IDENTICAL_REF),
    DisjoinProfile("DUPLICATE_TARGETS", _DUPLICATE_SELF, _DUPLICATE_REF),
    DisjoinProfile("MANY_CHROMOSOMES", _MANY_SELF, _MANY_REF),
    DisjoinProfile("SINGLE_CHROMOSOME", _SINGLE_SELF, _SINGLE_REF),
    DisjoinProfile("REFERENCE_GAP", _GAP_SELF, _GAP_REF),
    DisjoinProfile("REFERENCE_SUPERSET", _SUPERSET_SELF, _SUPERSET_REF),
    DisjoinProfile("EMPTY_TARGET", _EMPTY_TARGET_SELF, _EMPTY_TARGET_REF),
    DisjoinProfile("EMPTY_REFERENCE", _EMPTY_REFERENCE_SELF, _EMPTY_REFERENCE_REF),
    DisjoinProfile("CUSTOM_COLUMNS", _CUSTOM_SELF, _CUSTOM_REF),
    DisjoinProfile("ONE_BASED_CLOSED_TARGET", _ONE_BASED_SELF, _ONE_BASED_REF),
    DisjoinProfile("MIXED_ENCODING", _MIXED_SELF, _MIXED_REF),
)
