"""Target-engine model for GIQL transpilation.

A :class:`Target` is a first-class SQL target engine: it carries a
:class:`Capabilities` set describing what the engine supports and the
sqlglot output dialect used to serialize standard AST for that engine.

This is step 1 of epic #137. The targets and capability descriptors defined
here are the foundation that later steps build on — the operator-expander
registry is keyed by ``(target, operator)`` and emission choices become
capability lookups rather than scattered ``if dialect == ...`` branches.
At this step the model is wired into :func:`giql.transpile.transpile` only
to resolve the ``dialect`` parameter and drive the existing DuckDB IEJoin
gate; no emission behaviour changes.
"""

from dataclasses import dataclass
from typing import Literal

RangeJoinStrategy = Literal["binned", "iejoin"]


@dataclass(frozen=True)
class Capabilities:
    """Feature set of a SQL target engine.

    Each field is a portable choice that later steps of epic #137 turn into
    a capability lookup instead of a hardcoded dialect branch.

    Parameters
    ----------
    supports_lateral : bool
        Whether the engine supports ``LATERAL`` / correlated joins. Drives the
        NEAREST LATERAL-vs-window-function strategy (#142): a correlated NEAREST
        expands to a portable correlated ``LATERAL`` subquery where this holds
        and to a decorrelated window-function form where it does not. This
        capability is the single source of truth — the former
        ``BaseGIQLGenerator.SUPPORTS_LATERAL`` generator attribute has been
        removed.
    supports_star_replace : bool
        Whether the engine supports ``SELECT * REPLACE (...)``. Drives the
        coordinate-canonicalization output: ``* REPLACE`` where supported,
        an explicit portable projection otherwise (#143). Supported by
        DuckDB / BigQuery / Snowflake / ClickHouse; not by PostgreSQL,
        SQLite, or DataFusion.
    supports_qualify : bool
        Whether the engine supports the ``QUALIFY`` clause. Reserved: no
        emission path consumes it yet (a future window-function operator
        port would).
    range_join_strategy : RangeJoinStrategy
        The plan used for column-to-column INTERSECTS joins: ``"binned"``
        for the generic binned equi-join, ``"iejoin"`` for DuckDB's
        per-partition IEJoin plan. The IEJoin path covers INNER / SEMI /
        ANTI joins, with binned fallback for unsupported shapes.
    """

    supports_lateral: bool
    supports_star_replace: bool
    supports_qualify: bool
    range_join_strategy: RangeJoinStrategy


@dataclass(frozen=True)
class Target:
    """A SQL target engine.

    Subclasses declare the engine ``name``, the ``sqlglot_dialect`` used to
    serialize AST for that engine (``None`` selects sqlglot's default
    generic serialization), and the engine ``capabilities``.

    Targets are frozen, value-equal, and hashable: two ``DuckDBTarget()``
    instances compare equal and hash alike, so the operator-expander registry
    (#138) can key on a resolved target by value. Equality is class-scoped —
    ``GenericTarget() != DataFusionTarget()`` even where their fields overlap.
    """

    name: str
    sqlglot_dialect: str | None
    capabilities: Capabilities


@dataclass(frozen=True)
class GenericTarget(Target):
    """Portable SQL-92-ish target with no engine-specific features.

    This is the default target (``dialect=None``). Its capabilities are the
    conservative, maximally portable baseline that matches today's
    :class:`giql.generators.base.BaseGIQLGenerator` output.
    """

    name: str = "generic"
    sqlglot_dialect: str | None = None
    capabilities: Capabilities = Capabilities(
        supports_lateral=True,
        supports_star_replace=False,
        supports_qualify=False,
        range_join_strategy="binned",
    )


@dataclass(frozen=True)
class DuckDBTarget(Target):
    """DuckDB target.

    Serializes through sqlglot's ``duckdb`` dialect and uses the IEJoin
    per-partition plan for column-to-column INTERSECTS joins.
    """

    name: str = "duckdb"
    sqlglot_dialect: str | None = "duckdb"
    capabilities: Capabilities = Capabilities(
        supports_lateral=True,
        supports_star_replace=True,
        supports_qualify=True,
        range_join_strategy="iejoin",
    )


@dataclass(frozen=True)
class DataFusionTarget(Target):
    """Apache DataFusion target.

    sqlglot has no DataFusion dialect, so serialization falls back to the
    generic form (``sqlglot_dialect = None``) for now; #145 finalizes
    DataFusion serialization. The capability values below are conservative
    and provisional — they are validated against a real DataFusion engine
    when the operator migrations exercise them (#142, #145). DataFusion
    supports ``* EXCEPT`` / ``* EXCLUDE`` but not ``* REPLACE``.
    """

    name: str = "datafusion"
    sqlglot_dialect: str | None = None
    capabilities: Capabilities = Capabilities(
        supports_lateral=False,
        supports_star_replace=False,
        supports_qualify=False,
        range_join_strategy="binned",
    )


# Public dialect names only. ``generic`` is intentionally absent: ``None`` is
# the sole public way to select it (see :func:`resolve_target`).
_TARGETS_BY_NAME: dict[str, type[Target]] = {
    DuckDBTarget.name: DuckDBTarget,
    DataFusionTarget.name: DataFusionTarget,
}


def resolve_target(dialect: str | None) -> Target:
    """Resolve a ``dialect`` parameter to a :class:`Target` instance.

    Parameters
    ----------
    dialect : str | None
        The target dialect name. ``None`` resolves to :class:`GenericTarget`;
        ``"duckdb"`` and ``"datafusion"`` resolve to their respective targets.

    Returns
    -------
    Target
        The resolved target instance.

    Raises
    ------
    ValueError
        If *dialect* is not a recognized target name.
    """
    if dialect is None:
        return GenericTarget()

    target_cls = _TARGETS_BY_NAME.get(dialect)
    if target_cls is None:
        raise ValueError(
            f"Unknown dialect: {dialect!r}. Supported: 'duckdb', 'datafusion', or None."
        )
    return target_cls()
