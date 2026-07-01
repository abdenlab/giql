"""Target-engine model for GIQL transpilation.

A :class:`Target` is a first-class SQL target engine: it carries a
:class:`Capabilities` set describing what the engine supports and the
sqlglot output dialect used to serialize standard AST for that engine.

This is step 1 of epic #137. The targets and capability descriptors defined
here are the foundation that later steps build on — the operator-expander
registry is keyed by ``(target, operator)`` and emission choices are
capability lookups rather than scattered ``if dialect == ...`` branches.
As of #143/#145 the model actively drives emission: the INTERSECTS join
strategy (``range_join_strategy``), the NEAREST LATERAL-vs-window form
(``supports_lateral``), and the coordinate-canonicalization wrapper and
NEAREST/DISJOIN passthroughs (``supports_star_replace``) all read from the
active target's capabilities.
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
        capability is the single source of truth — the former generator-level
        ``SUPPORTS_LATERAL`` attribute has been removed.
    supports_star_replace : bool
        Whether the engine supports ``SELECT * REPLACE (...)``. Drives every
        coordinate-canonicalization site — the canonicalizer wrapper CTE and the
        DISJOIN / NEAREST passthroughs (#143/#145): ``* REPLACE`` where supported,
        the portable ``* EXCEPT`` projection otherwise. Supported by
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
    conservative, maximally portable baseline, serialized through sqlglot's
    dialect-less default (``sqlglot_dialect = None``).

    "SQL-92-ish", not strict SQL-92: because ``supports_star_replace=False``,
    every coordinate-canonicalization site over a **non-canonical** target falls
    back to a ``SELECT * EXCEPT (...)`` projection (re-appending the recomputed
    interval columns) — the canonicalizer wrapper CTE and the DISJOIN / NEAREST
    passthroughs alike. ``* EXCEPT`` is **not** SQL-92 and is **not
    DuckDB-runnable** — it is a DataFusion-family extension — so the generic
    target's non-canonical output runs only on an ``* EXCEPT``-capable engine.
    A canonical (0-based half-open) target passes the row through as a plain,
    fully portable ``SELECT *``.
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

    Serializing through the ``duckdb`` dialect makes null ordering **explicit** on
    ``ORDER BY`` / window terms, emitting ``NULLS FIRST`` where the generic
    (dialect-less) serialization leaves it implicit. This is result-preserving for
    every migrated operator because each sorts on a non-null coordinate key
    (chromosome / position, or a same-chromosome-filtered distance), so the
    explicit ``NULLS FIRST`` never reorders a NULL into or out of the top rows. A
    future operator ordering on a nullable key would need to re-verify this before
    relying on the DuckDB serialization being a semantic no-op.
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

    sqlglot has no DataFusion dialect, so serialization uses the generic form
    (``sqlglot_dialect = None``); this is the finalized strategy (#145) — the
    portable SQL the generic generator emits runs on DataFusion, verified
    end-to-end by the cross-target oracle (``tests/integration/datafusion/``):
    every operator at the default encoding, and the canonicalizing operators
    (DISJOIN / NEAREST) across all four coordinate encodings plus custom-column
    and strand schemas.

    The capability values below are validated against a real DataFusion engine by
    that oracle: ``supports_lateral=False`` (no correlated-LATERAL physical plan,
    so NEAREST takes the decorrelated window fallback), ``supports_star_replace=
    False`` (DataFusion supports ``* EXCEPT`` / ``* EXCLUDE`` but not
    ``* REPLACE``, so the canonicalizer and the NEAREST/DISJOIN passthroughs emit
    the portable ``* EXCEPT`` form), ``supports_qualify=False``, and the binned
    equi-join range strategy.

    One documented gap remains (#160, dependent on #146): a ``SELECT *`` /
    ``SELECT b.*`` over a correlated NEAREST exposes the fallback's reserved
    ``__giql_x_*`` columns, so the cross-target identity claim is narrowed to
    explicitly-projected queries until a query-level projection seam lands.
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
        ``"duckdb"`` and ``"datafusion"`` resolve to their respective built-in
        targets. Any other name is resolved against the plugin registry — a
        custom :class:`Target` declared through
        :meth:`giql.expander.ExpanderRegistry.register_target` (or as a side
        effect of :func:`giql.expander.register`) is selectable by its ``name``.

    Returns
    -------
    Target
        The resolved target instance.

    Raises
    ------
    ValueError
        If *dialect* is neither a built-in name nor a registered custom target.
    """
    if dialect is None:
        return GenericTarget()

    target_cls = _TARGETS_BY_NAME.get(dialect)
    if target_cls is not None:
        return target_cls()

    # A non-built-in name may be a custom target registered on the plugin hub.
    # Imported lazily so this module stays import-cycle-free: ``giql.expander``
    # imports ``Target`` / ``GenericTarget`` from here at module load.
    from giql.expander import REGISTRY

    registered = REGISTRY.target(dialect)
    if registered is not None:
        return registered

    raise ValueError(
        f"Unknown dialect: {dialect!r}. Supported: 'duckdb', 'datafusion', None, "
        "or a custom target registered via "
        "giql.expander.REGISTRY.register_target()."
    )
