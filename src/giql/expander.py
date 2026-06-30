"""The operator-expander registry and the ``ExpandOperators`` pass (epic #137).

This module is step 2 of epic #137. It lands the dispatch infrastructure that the
remaining steps migrate each operator onto, one at a time:

* an :class:`OperatorExpander` protocol — ``expand(node, ctx) -> exp.Expression``,
  the unit that turns one GIQL operator node into standard sqlglot AST for the
  active target;
* an :class:`ExpansionContext` carrying everything an expander needs — the pass-1
  :class:`~giql.resolver.OperatorResolution` metadata, the active
  :class:`~giql.targets.Target` and its :class:`~giql.targets.Capabilities`,
  collision-safe alias minting, and the registered :class:`~giql.table.Tables`;
* a :class:`ExpanderRegistry` keyed by ``(target, operator_type)`` resolving an
  expander through the fallback chain ``(target, op)`` → ``(generic, op)`` →
  the legacy ``*_sql`` emitter;
* a :func:`register` decorator — the **public extension hook** by which a user
  adds a target or overrides an operator for one, e.g.
  ``@register(DuckDBTarget, GIQLDisjoin)``;
* the registry's teardown seam — :meth:`ExpanderRegistry.unregister`,
  :meth:`ExpanderRegistry.clear`, and ``in`` membership (``__contains__``) — the
  **supported public API** for undoing a custom registration. A user who
  registers a custom expander (e.g. in a test fixture or a plugin's teardown)
  resets the process-wide :data:`REGISTRY` through these rather than reaching
  into private state: ``REGISTRY.unregister(DuckDBTarget(), GIQLDisjoin)`` drops
  one entry, ``REGISTRY.clear()`` drops them all, and
  ``(DuckDBTarget(), GIQLDisjoin) in REGISTRY`` tests for one;
* the :class:`ExpandOperators` pass, which runs after
  :func:`giql.canonicalizer.canonicalize_coordinates`, walks the AST, and
  replaces every opted-in GIQL operator node with the registry's expansion.

Gating (mirrors ``GIQL_CANONICALIZE``)
--------------------------------------
The pass is gated per operator by a ``GIQL_EXPAND`` class attribute on the
operator's expression class, exactly as pass 2 is gated by ``GIQL_CANONICALIZE``.
An operator takes the new AST-expansion path **only when both** hold:

* its class sets ``GIQL_EXPAND = True``, and
* the registry resolves an expander for ``(active target, operator type)`` —
  i.e. a ``(target, op)`` or ``(generic, op)`` expander is registered.

Otherwise it falls through to the legacy ``*_sql`` emitter on
:class:`giql.generators.base.BaseGIQLGenerator`. The built-in expanders register
at import time via :mod:`giql.expanders`; the pass rewrites a node only when
``GIQL_EXPAND=True`` **and** an expander resolves for ``(active target, operator
type)``, and is a no-op for any operator that is unflagged or has no registered
expander. A migration PR registers an expander, flips one operator's
``GIQL_EXPAND`` flag, and deletes that operator's ``*_sql`` method.
"""

from __future__ import annotations

from typing import Callable
from typing import Protocol
from typing import runtime_checkable

from sqlglot import exp
from sqlglot.helper import name_sequence

from giql.resolver import META_KEY
from giql.resolver import OperatorResolution
from giql.resolver import ResolutionError
from giql.table import Tables
from giql.targets import GenericTarget
from giql.targets import Target

__all__ = [
    "EXPAND_ALIAS_PREFIX",
    "ExpansionContext",
    "OperatorExpander",
    "ExpanderRegistry",
    "REGISTRY",
    "register",
    "expand_operators",
    "ExpandOperators",
]

#: The prefix every alias minted by :meth:`ExpansionContext.alias` carries. The
#: leading double underscore keeps synthesized names clear of user identifiers,
#: mirroring the ``__giql_canon_`` / ``__giql_dj_`` reserved prefixes.
EXPAND_ALIAS_PREFIX = "__giql_x_"


class ExpansionContext:
    """Everything an :class:`OperatorExpander` needs to expand one operator.

    A fresh context is built per operator node by :class:`ExpandOperators`. It
    bundles the node's pass-1 resolution metadata together with the active
    target, collision-safe alias minting, and the registered tables, so an
    expander never has to reach back into the pass or re-resolve anything.

    Attributes
    ----------
    node : exp.Expression
        The GIQL operator node being expanded.
    resolution : OperatorResolution | None
        The metadata :func:`giql.resolver.resolve_operator_refs` attached to
        *node* (its resolved reference slots, column operands, deferrals, and
        de-canonicalization output tables). ``None`` only if the node was never
        annotated — an internal invariant violation an expander may assert on.
    target : Target
        The active :class:`~giql.targets.Target`, exposing its
        ``capabilities`` and ``sqlglot_dialect``.
    tables : Tables
        The registered :class:`~giql.table.Tables` container.
    registry : ExpanderRegistry | None
        The registry the pass is resolving against. Carried so a whole-query
        rewrite (CLUSTER / MERGE) can re-enter :func:`expand_operators` over the
        SELECT it just restructured and expand sibling operators it copied into
        it, honoring a custom-registry pass run. ``None`` for a standalone
        context built outside the pass.
    """

    __slots__ = ("node", "resolution", "target", "tables", "registry", "_alias_seq")

    def __init__(
        self,
        node: exp.Expression,
        resolution: OperatorResolution | None,
        target: Target,
        tables: Tables,
        alias_seq: Callable[[], str] | None = None,
        registry: ExpanderRegistry | None = None,
    ) -> None:
        self.node = node
        self.resolution = resolution
        self.target = target
        self.tables = tables
        self.registry = registry
        # A single sequence is threaded across every context built for one
        # ``ExpandOperators`` run so aliases minted for sibling operators never
        # collide; a standalone context falls back to its own sequence.
        self._alias_seq = alias_seq or name_sequence(EXPAND_ALIAS_PREFIX)

    @property
    def capabilities(self):
        """The active target's :class:`~giql.targets.Capabilities`."""
        return self.target.capabilities

    def alias(self) -> str:
        """Mint a fresh, query-unique alias with the reserved expander prefix.

        Draws sequential names from a per-run :func:`sqlglot.helper.name_sequence`
        so two expanders (or two slots of one expander) never collide.
        """
        return self._alias_seq()


@runtime_checkable
class OperatorExpander(Protocol):
    """Turns one GIQL operator node into standard sqlglot AST for a target.

    An expander is any callable-bearing object whose :meth:`expand` takes the
    operator node and its :class:`ExpansionContext` and returns the sqlglot AST
    the node expands to for the active target. The returned expression replaces
    the operator node in the tree; the stock per-target serializer renders it.

    The protocol is :func:`~typing.runtime_checkable`, so the registry can assert
    a registered object satisfies it. A plain function is *not* an
    ``OperatorExpander`` (it has no ``expand`` method); register one by wrapping
    it (see :func:`register`, which accepts either form).

    An expander is **node-local**: ``expand(node, ctx) -> exp.Expression`` sees
    one operator node and returns the expression that replaces it in place. It
    cannot express a whole-query rewrite such as the INTERSECTS IEJoin fold,
    which restructures the surrounding query (joins, CTEs) rather than a single
    node. That fold is therefore deferred — it would need a separate
    query-level mechanism — and is handled by the pre-pass join transformers, not
    by an expander.
    """

    def expand(self, node: exp.Expression, ctx: ExpansionContext) -> exp.Expression: ...


#: A bare expander function — the form most expanders are written in. The
#: registry stores either an :class:`OperatorExpander` object or one of these.
ExpanderFn = Callable[[exp.Expression, ExpansionContext], exp.Expression]


def _as_callable(expander: OperatorExpander | ExpanderFn) -> ExpanderFn:
    """Normalize an expander to a plain ``(node, ctx) -> Expression`` callable."""
    # ``@runtime_checkable`` only checks an ``expand`` attribute *exists*, so an
    # object with a non-callable ``expand`` would pass the isinstance check; guard
    # that the attribute is actually callable before trusting the protocol branch.
    if isinstance(expander, OperatorExpander) and callable(expander.expand):
        return expander.expand
    if callable(expander):
        return expander
    raise TypeError(
        "An operator expander must be an OperatorExpander (an object with an "
        f"'expand' method) or a callable; got {type(expander).__name__}."
    )


class ExpanderRegistry:
    """Maps ``(target, operator_type)`` to the expander that handles it.

    Resolution follows the epic's fallback chain so a generic expander written
    once covers every target, and a per-target entry overrides it only where the
    engine genuinely differs:

    1. ``(target, op)`` — an expander registered for the *exact* active target;
    2. ``(generic, op)`` — the portable expander registered for
       :class:`~giql.targets.GenericTarget`;
    3. ``None`` — no expander, so the caller keeps the legacy ``*_sql`` emitter.

    The registry keys on the frozen, value-equal :class:`~giql.targets.Target`
    *instance* (two ``DuckDBTarget()`` compare and hash alike), so registering
    against ``DuckDBTarget`` and resolving against a freshly built
    ``DuckDBTarget()`` hit the same entry. Registering against the *class*
    ``DuckDBTarget`` is a convenience the decorator instantiates.
    """

    def __init__(self) -> None:
        self._expanders: dict[tuple[Target, type], ExpanderFn] = {}

    def register(
        self,
        target: Target,
        operator: type,
        expander: OperatorExpander | ExpanderFn,
    ) -> None:
        """Register *expander* for ``(target, operator)``, overriding any prior.

        Parameters
        ----------
        target : Target
            The target instance the expander handles. Use
            :class:`~giql.targets.GenericTarget` for the portable fallback.
        operator : type
            The GIQL operator expression class (e.g.
            :class:`~giql.expressions.GIQLDisjoin`).
        expander : OperatorExpander | ExpanderFn
            The expander object or function. A later registration for the same
            key replaces an earlier one (last-write-wins override).

        Notes
        -----
        A *non-generic* ``(target, operator)`` entry is intended to also act as a
        join-rewrite override for operators with a built-in whole-query join
        rewrite (notably :class:`~giql.expressions.Intersects`, whose binned
        equi-join / DuckDB IEJoin transformers run before expansion), letting a
        per-target expander assume responsibility for that rewrite. That bypass
        is intended for a future INTERSECTS consumer and is **not wired by any
        caller yet** — no transformer consults :meth:`has_override` here (see
        #141).
        """
        self._expanders[(target, operator)] = _as_callable(expander)

    def resolve(self, target: Target, operator: type) -> ExpanderFn | None:
        """Return the expander for ``(target, operator)`` via the fallback chain.

        Tries the exact ``(target, op)`` entry, then the
        ``(GenericTarget(), op)`` fallback, then ``None`` (legacy emitter).

        A non-generic exact ``(target, op)`` entry is also intended to act as a
        *join-rewrite override* for operators with a built-in whole-query join
        rewrite (notably :class:`~giql.expressions.Intersects`). That override is
        intended for a future INTERSECTS consumer and is **not wired by any
        caller yet** — resolution does not itself bypass the built-in
        binned / IEJoin transformers (see :meth:`register`, :meth:`has_override`,
        and #141).
        """
        fn = self._expanders.get((target, operator))
        if fn is not None:
            return fn
        if target != GenericTarget():
            fn = self._expanders.get((GenericTarget(), operator))
            if fn is not None:
                return fn
        return None

    def has_override(self, target: Target, operator: type) -> bool:
        """Whether an exact non-generic ``(target, operator)`` entry is registered.

        Returns ``True`` only when *target* is not :class:`~giql.targets.GenericTarget`
        and an exact ``(target, operator)`` entry is registered; the portable
        ``(GenericTarget(), operator)`` fallback is *not* an override and does not
        count here.

        Such an entry is intended to mark a target-specific override that
        supersedes built-in handling (e.g. taking responsibility for the
        whole-query join rewrite the built-in transformers would otherwise
        perform). That mechanism is intended for a future INTERSECTS consumer and
        is **not wired by any caller yet** — no transformer consults this method
        in the current pipeline (see #141).
        """
        return target != GenericTarget() and (target, operator) in self._expanders

    def unregister(self, target: Target, operator: type) -> None:
        """Drop the ``(target, operator)`` entry if present.

        Part of the registry's **public teardown seam**: a caller that
        registered a custom expander (a plugin, or a test fixture) undoes one
        registration with this rather than mutating private state. Resolving the
        same key afterward falls back through the chain to ``None``. Dropping an
        absent key is a no-op.
        """
        self._expanders.pop((target, operator), None)

    def clear(self) -> None:
        """Drop every registration.

        The bulk form of :meth:`unregister`, and the **public reset** for the
        process-wide :data:`REGISTRY` — e.g. a test fixture that saves and
        restores the registry around a body that registers custom expanders.
        """
        self._expanders.clear()

    def snapshot(self) -> dict[tuple[Target, type], ExpanderFn]:
        """Return a shallow copy of the current registrations.

        The save half of a save/restore seam that supports test
        baseline-isolation: capture the baseline with this and hand it back to
        :meth:`restore` afterward, so the built-in expanders registered at import
        survive an isolating fixture that would otherwise :meth:`clear` them
        permanently.

        The returned dict is a fresh mapping (mutating it does not affect the
        registry), keyed by the same ``(target, operator)`` tuples.
        """
        return dict(self._expanders)

    def restore(self, snapshot: dict[tuple[Target, type], ExpanderFn]) -> None:
        """Replace all registrations with those captured by :meth:`snapshot`.

        The restore half of the save/restore seam that supports test
        baseline-isolation. Drops every current entry and re-installs exactly the
        *snapshot* contents, so a fixture can return the registry to a previously
        captured baseline regardless of what its body registered or cleared.
        """
        self._expanders.clear()
        self._expanders.update(snapshot)

    def __contains__(self, key: tuple[Target, type]) -> bool:
        """Whether an *exact* ``(target, operator)`` entry is registered.

        Tests the exact key only — it does **not** walk the generic fallback
        chain :meth:`resolve` follows (``(GenericTarget(), op) in registry`` is
        the only way ``(target, op)`` reports present via the generic entry).
        Part of the public teardown seam, for asserting a registration landed or
        was torn down.
        """
        return key in self._expanders

    def __len__(self) -> int:
        """The number of registered ``(target, operator)`` entries.

        Part of the public introspection surface (alongside ``in``): emptiness is
        observable without reaching into private state, so a teardown fixture or
        leak guard can assert the registry is clear via ``len(registry)`` or
        ``bool(registry)``.
        """
        return len(self._expanders)

    def __bool__(self) -> bool:
        """Whether any expander is registered (``False`` when empty)."""
        return bool(self._expanders)


#: The process-wide registry the :func:`register` decorator writes to and the
#: :class:`ExpandOperators` pass reads from. The built-in expanders register into
#: it at import time via :mod:`giql.expanders`; the pass rewrites a node only when
#: an expander resolves here (and the operator is flagged ``GIQL_EXPAND``).
REGISTRY = ExpanderRegistry()


def register(
    target: type[Target] | Target, operator: type
) -> Callable[[OperatorExpander | ExpanderFn], OperatorExpander | ExpanderFn]:
    """Register an operator expander for ``(target, operator)`` — the extension hook.

    The public decorator by which a built-in or user-supplied expander is added
    to the process-wide :data:`REGISTRY`::

        @register(DuckDBTarget, GIQLDisjoin)
        def expand_disjoin(node, ctx): ...

    Parameters
    ----------
    target : type[Target] | Target
        The target the expander handles, given as the target *class* (the usual
        form — it is instantiated with its default capabilities) or an explicit
        instance. Pass :class:`~giql.targets.GenericTarget` for the portable
        expander every target falls back to.
    operator : type
        The GIQL operator expression class the expander handles.

    Returns
    -------
    Callable
        A decorator that registers its argument and returns it unchanged, so the
        decorated object stays usable on its own.
    """
    if isinstance(target, type):
        try:
            target_instance: Target = target()
        except TypeError as exc:
            raise TypeError(
                f"register() could not instantiate target class {target.__name__!r} "
                "with no arguments: a custom Target subclass with required fields "
                "cannot be defaulted. Pass an instance (e.g. "
                f"@register({target.__name__}(...), ...)) or give its fields "
                "defaults so the class form can construct it."
            ) from exc
    else:
        target_instance = target

    def decorator(
        expander: OperatorExpander | ExpanderFn,
    ) -> OperatorExpander | ExpanderFn:
        REGISTRY.register(target_instance, operator, expander)
        return expander

    return decorator


# The GIQL operator expression classes the pass inspects. Membership alone does
# not opt an operator in: the per-class ``GIQL_EXPAND`` flag plus a registered
# expander do (see module docstring). Imported lazily as a precaution; there is no
# real cycle (``canonicalizer.py`` imports the same classes eagerly at module
# level), so a later step may hoist these to module scope to match that sibling.
def _giql_operators() -> tuple[type, ...]:
    """Return the GIQL operator classes, imported lazily as a precaution."""
    from giql.expressions import Contains
    from giql.expressions import GIQLCluster
    from giql.expressions import GIQLDisjoin
    from giql.expressions import GIQLDistance
    from giql.expressions import GIQLMerge
    from giql.expressions import GIQLNearest
    from giql.expressions import Intersects
    from giql.expressions import SpatialSetPredicate
    from giql.expressions import Within

    return (
        GIQLDisjoin,
        GIQLNearest,
        GIQLDistance,
        GIQLCluster,
        GIQLMerge,
        Intersects,
        Contains,
        Within,
        SpatialSetPredicate,
    )


def expand_operators(
    expression: exp.Expression,
    target: Target,
    tables: Tables,
    registry: ExpanderRegistry | None = None,
) -> exp.Expression:
    """Replace each opted-in GIQL operator with its registry expansion.

    Pass 3 of the normalization pipeline (epic #137). Runs after
    :func:`giql.canonicalizer.canonicalize_coordinates`. For every GIQL operator
    node it dispatches to the new AST-expansion path **only when** the operator's
    class sets ``GIQL_EXPAND = True`` *and* the registry resolves an expander for
    ``(target, operator type)`` through its fallback chain; otherwise the node is
    left untouched and the legacy ``*_sql`` emitter handles it.

    The pass mutates and returns *expression* in place. It touches only nodes
    whose operator is flagged ``GIQL_EXPAND`` and resolves an expander; for every
    other operator it is a no-op, leaving the emitted SQL byte-identical, so the
    existing suite is the migration oracle.

    Parameters
    ----------
    expression : exp.Expression
        The pass-1/pass-2-annotated AST.
    target : Target
        The active target whose expanders to dispatch to.
    tables : Tables
        The registered table configurations, threaded into each context.
    registry : ExpanderRegistry | None
        The registry to resolve against; defaults to the process-wide
        :data:`REGISTRY`.

    Returns
    -------
    exp.Expression
        The same *expression*, with each opted-in operator node that resolves an
        expander replaced by its target-specific expansion; nodes that are
        unflagged or resolve no expander are left untouched.
    """
    reg = registry if registry is not None else REGISTRY
    operators = _giql_operators()
    alias_seq = name_sequence(EXPAND_ALIAS_PREFIX)

    # Collect first, then mutate: replacing nodes mid-walk is unsafe.
    pending: list[tuple[exp.Expression, ExpanderFn]] = []
    for node in expression.walk():
        if not isinstance(node, operators):
            continue
        if not getattr(node, "GIQL_EXPAND", False):
            continue
        fn = reg.resolve(target, type(node))
        if fn is None:
            continue
        pending.append((node, fn))

    # Replace deepest-first: replacing an ancestor node detaches any collected
    # descendant, whose later ``node.replace`` would be a silent no-op (the
    # detached node's parent is gone), so the inner operator would never expand.
    # Mutating leaves first keeps every still-collected ancestor attached. A
    # detached node is skipped defensively as a second line of defence (#154).
    pending.sort(key=lambda item: _node_depth(item[0]), reverse=True)

    for node, fn in pending:
        if node is not expression and node.parent is None:
            # An ancestor expansion already detached this node from the tree;
            # replacing it would be a no-op, so skip it.
            continue
        resolution = node.meta.get(META_KEY)
        if not isinstance(resolution, OperatorResolution):
            # Pass 1 guarantees every operator node carries its resolution; a
            # missing or malformed one is an internal invariant violation, not a
            # user error, so fail loudly rather than minting a None-resolution
            # context the expander would then dereference blindly.
            raise ResolutionError(
                f"{type(node).__name__} reached the ExpandOperators pass without "
                "valid resolution metadata; pass 1 (resolve_operator_refs) must "
                "run first and annotate every operator node."
            )
        ctx = ExpansionContext(node, resolution, target, tables, alias_seq, registry=reg)
        replacement = fn(node, ctx)
        if not isinstance(replacement, exp.Expression):
            raise TypeError(
                f"expander {fn!r} for {type(node).__name__} returned "
                f"{type(replacement).__name__}, not exp.Expression"
            )
        if replacement is not node:
            node.replace(replacement)

    return expression


def _node_depth(node: exp.Expression) -> int:
    """The number of ancestors above *node* (root is depth 0).

    Used to order the collect-then-replace loop deepest-first so replacing an
    ancestor never detaches a still-pending descendant before its own replace.
    """
    depth = 0
    parent = node.parent
    while parent is not None:
        depth += 1
        parent = parent.parent
    return depth


class ExpandOperators:
    """Callable wrapper for :func:`expand_operators`, parallel to the transformers.

    The transpile pipeline composes transformer *objects* bound to their
    :class:`~giql.table.Tables`; this wrapper gives the expansion pass the same
    shape — construct it with the active target and tables, then call it on the
    AST — so it slots into the pipeline uniformly after
    :class:`giql.canonicalizer.canonicalize_coordinates`.
    """

    def __init__(
        self,
        target: Target,
        tables: Tables,
        registry: ExpanderRegistry | None = None,
    ) -> None:
        self.target = target
        self.tables = tables
        self.registry = registry

    def transform(self, expression: exp.Expression) -> exp.Expression:
        """Run the expansion pass over *expression* and return it."""
        return expand_operators(expression, self.target, self.tables, self.registry)
