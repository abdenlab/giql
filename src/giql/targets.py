"""Target-engine model for GIQL transpilation.

A :class:`Target` is a first-class SQL target engine: it carries a
:class:`Capabilities` set describing what the engine supports and the
sqlglot output dialect used to serialize standard AST for that engine.

This is step 1 of epic #137. The targets and capability descriptors defined
here are the foundation that later steps build on — the operator-expander
registry is keyed by ``(target, operator)`` and emission choices are
capability lookups rather than scattered ``if dialect == ...`` branches.
As of #143/#145 the model actively drives emission: the NEAREST LATERAL-vs-window
form (``supports_lateral``) and the coordinate-canonicalization wrapper and
NEAREST/DISJOIN passthroughs (``supports_star_replace``) read from the active
target's capabilities. The INTERSECTS join strategy is instead selected by the
operator-expander registry (#169): DuckDB registers a ``(DuckDBTarget, Intersects)``
IEJoin override, and every other target falls back to the generic naive-predicate
expander — no capability flag is consulted.
"""

from dataclasses import dataclass


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
    """

    supports_lateral: bool
    supports_star_replace: bool
    supports_qualify: bool


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
    passthroughs alike. A star-projected CLUSTER is a further ``* EXCEPT`` site
    (independent of canonicalization): it hides its synthesized
    ``__giql_is_new_cluster`` flag with ``* EXCEPT`` (#184). ``* EXCEPT`` is **not**
    SQL-92 and is **not DuckDB-runnable** — it is a DataFusion-family extension — so
    the generic target's output from any of these sites runs only on an
    ``* EXCEPT``-capable engine (transpile with ``dialect="duckdb"`` to execute on
    DuckDB, where it spells the exclusion ``EXCLUDE``). A canonical (0-based
    half-open) target with no star-projected CLUSTER passes the row through as a
    plain, fully portable ``SELECT *``.
    """

    name: str = "generic"
    sqlglot_dialect: str | None = None
    capabilities: Capabilities = Capabilities(
        supports_lateral=True,
        supports_star_replace=False,
        supports_qualify=False,
    )


@dataclass(frozen=True)
class DuckDBTarget(Target):
    """DuckDB target.

    Serializes through sqlglot's ``duckdb`` dialect and uses the IEJoin
    per-partition plan for column-to-column INTERSECTS joins, registered as the
    ``(DuckDBTarget, Intersects)`` override in the operator-expander registry
    (:mod:`giql.expanders.intersects_duckdb`, #169).

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
    the portable ``* EXCEPT`` form), and ``supports_qualify=False``. DataFusion
    registers no INTERSECTS override, so a column-to-column INTERSECTS join falls
    back to the generic naive overlap predicate — a plain ``ON`` condition
    DataFusion plans as a hash join keyed on ``chrom`` with the position
    inequalities as a residual join filter (#167).

    A ``SELECT *`` / ``SELECT b.*`` over a correlated NEAREST would otherwise expose
    the decorrelated fallback's reserved ``__giql_x_*`` rank/key columns; a statement
    finalizer wraps the enclosing ``SELECT`` in ``SELECT * EXCEPT (...)`` on the
    fallback path (#160), so the cross-target identity claim holds for star
    projections over a correlated NEAREST unconditionally. The finalizer re-locates
    its target join by a reserved ``meta`` tag and mints its reserved column names
    per fallback (#172), so the claim now holds through the two former residual
    compositions too: a correlated NEAREST re-surfaced by an enclosing ``SELECT *``
    *outside* its own SELECT (e.g. a wrapping ``CLUSTER``, whose ``copy()`` +
    ``transplant`` the tag survives) and two correlated NEAREST fallbacks in one
    query (whose now-distinct per-run reserved names each ``* EXCEPT`` independently).
    """

    name: str = "datafusion"
    sqlglot_dialect: str | None = None
    capabilities: Capabilities = Capabilities(
        supports_lateral=False,
        supports_star_replace=False,
        supports_qualify=False,
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
