"""The operator-expander registry and the ``ExpandOperators`` pass (epic #137).

This module provides the dispatch infrastructure every GIQL operator expands
through (epic #137, complete):

* an :class:`OperatorExpander` protocol — ``expand(node, ctx) -> exp.Expression``,
  the unit that turns one GIQL operator node into standard sqlglot AST for the
  active target;
* an :class:`ExpansionContext` carrying everything an expander needs — the pass-1
  :class:`~giql.resolver.OperatorResolution` metadata, the active
  :class:`~giql.targets.Target` and its :class:`~giql.targets.Capabilities`,
  collision-safe alias minting, and the registered :class:`~giql.table.Tables`;
* a :class:`ExpanderRegistry` keyed by ``(target, operator_type)`` resolving an
  expander through the fallback chain ``(target, op)`` → ``(generic, op)``;
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
  replaces every GIQL operator node with the registry's expansion.

Dispatch
--------
With every operator migrated (epic #137 complete), the pass rewrites *every* GIQL
operator node: it resolves an expander for ``(active target, operator type)``
through the registry's fallback chain and replaces the node with the AST that
expander returns. The built-in expanders register at import time via
:mod:`giql.expanders`, so a built-in operator always resolves at least a
``(generic, op)`` expander; a resolution miss is therefore an internal invariant
violation (an operator with no registered expander) and the pass raises rather
than leaving an un-serializable GIQL node in the tree. There is no per-operator
opt-in flag and no legacy ``*_sql`` fallback — both were migration scaffolding
removed once every operator was migrated (#146).
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Callable
from typing import Protocol
from typing import runtime_checkable

from sqlglot import exp
from sqlglot.helper import name_sequence

from giql.expressions import Contains
from giql.expressions import GIQLCluster
from giql.expressions import GIQLDisjoin
from giql.expressions import GIQLDistance
from giql.expressions import GIQLMerge
from giql.expressions import GIQLNearest
from giql.expressions import Intersects
from giql.expressions import SpatialSetPredicate
from giql.expressions import Within
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
    "StatementFinalizer",
    "ExpanderRegistry",
    "RegistrySnapshot",
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

    A node-local expander that needs to rewrite the *enclosing* statement (rather
    than just replace its own node) registers a :data:`StatementFinalizer` via
    :meth:`add_statement_finalizer`; the pass applies every finalizer to the
    statement after all node-local replacements. This is the query-level seam for
    a target whose expansion must reshape the enclosing statement — for example to
    project internal helper columns out of a surfacing ``SELECT *``.
    """

    __slots__ = (
        "node",
        "resolution",
        "target",
        "tables",
        "registry",
        "_alias_seq",
        "_finalizers",
    )

    def __init__(
        self,
        node: exp.Expression,
        resolution: OperatorResolution | None,
        target: Target,
        tables: Tables,
        alias_seq: Callable[[], str] | None = None,
        registry: ExpanderRegistry | None = None,
        finalizers: list[StatementFinalizer] | None = None,
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
        # A single finalizer list is likewise shared across one run's contexts so
        # a finalizer registered while expanding one node is applied once, after
        # every node-local replacement; a standalone context gets its own (inert
        # unless someone drives it manually).
        self._finalizers = finalizers if finalizers is not None else []

    @property
    def capabilities(self):
        """The active target's :class:`~giql.targets.Capabilities`."""
        return self.target.capabilities

    def add_statement_finalizer(self, finalizer: StatementFinalizer) -> None:
        """Register a query-level :data:`StatementFinalizer` for this run.

        The **query-level seam**: an expander is node-local — its return value
        replaces only its own node — so a target that must rewrite the *enclosing*
        statement (for example to project internal helper columns out of a
        surfacing ``SELECT *``) registers a finalizer here instead.
        :func:`expand_operators` applies every registered finalizer to the
        statement, in registration order, *after* all node-local replacements
        complete; each receives the current statement root and returns the
        (possibly new) root.

        A finalizer's returned root is emitted verbatim — beyond a type check that
        it is an :class:`~sqlglot.expressions.Expression`, it is **not**
        semantically re-validated — so it must not reference columns or relations
        absent from the projection it rewrites: a wrapper over an absent column
        builds without error but fails at engine runtime. Finalizers registered on a standalone context (one built
        outside the pass) are collected but never applied.
        """
        self._finalizers.append(finalizer)

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

    An expander's **return value** is node-local: ``expand(node, ctx) ->
    exp.Expression`` returns the one expression that replaces the operator node in
    place. When a target additionally needs to rewrite the *enclosing* statement —
    for example to project internal helper columns away from a surfacing
    ``SELECT *`` — the expander registers a query-level :data:`StatementFinalizer`
    via :meth:`ExpansionContext.add_statement_finalizer`, applied to the statement
    after every node-local replacement. The INTERSECTS IEJoin whole-query fold is a
    separate concern still handled by the pre-pass join transformers, not by an
    expander.
    """

    def expand(self, node: exp.Expression, ctx: ExpansionContext) -> exp.Expression: ...


#: A bare expander function — the form most expanders are written in. The
#: registry stores either an :class:`OperatorExpander` object or one of these.
ExpanderFn = Callable[[exp.Expression, ExpansionContext], exp.Expression]

#: A query-level statement finalizer: ``finalize(root) -> root``. An expander
#: registers one via :meth:`ExpansionContext.add_statement_finalizer` to wrap or
#: rewrite the enclosing statement after every node-local replacement; the pass
#: applies each in registration order and threads the (possibly new) root through.
#: Used, for example, to project internal helper columns out of a surfacing
#: ``SELECT *`` / ``b.*``.
StatementFinalizer = Callable[[exp.Expression], exp.Expression]


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


@dataclass(frozen=True)
class RegistrySnapshot:
    """An opaque save/restore token for :class:`ExpanderRegistry` state.

    Bundles the ``(target, operator)`` expander map and the declared-target map
    captured by :meth:`ExpanderRegistry.snapshot`. Hand it back to
    :meth:`ExpanderRegistry.restore` to return the registry to that state, or
    compare two snapshots for equality to assert the registry returned to a
    baseline (the leak-guard use).

    A snapshot is compared by equality and handed to :meth:`restore` only; it is
    **not** hashable (its two ``dict`` fields are mutable), so ``hash(snapshot)``
    raises — never use one as a set member or dict key.
    """

    expanders: dict[tuple[Target, type], ExpanderFn]
    targets: dict[str, Target]


class ExpanderRegistry:
    """Maps ``(target, operator_type)`` to the expander that handles it.

    Resolution follows the epic's fallback chain so a generic expander written
    once covers every target, and a per-target entry overrides it only where the
    engine genuinely differs:

    1. ``(target, op)`` — an expander registered for the *exact* active target;
    2. ``(generic, op)`` — the portable expander registered for
       :class:`~giql.targets.GenericTarget`;
    3. ``None`` — no expander resolves. Every built-in operator registers a
       ``(generic, op)`` expander, so this only arises for a cleared registry or
       an operator with no generic fallback; the pass raises rather than emitting
       (there is no legacy ``*_sql`` fallback).

    The registry keys on the frozen, value-equal :class:`~giql.targets.Target`
    *instance* (two ``DuckDBTarget()`` compare and hash alike), so registering
    against ``DuckDBTarget`` and resolving against a freshly built
    ``DuckDBTarget()`` hit the same entry. Registering against the *class*
    ``DuckDBTarget`` is a convenience the decorator instantiates.
    """

    def __init__(self) -> None:
        self._expanders: dict[tuple[Target, type], ExpanderFn] = {}
        self._targets_by_name: dict[str, Target] = {}

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
        Registering an expander also declares *target* by name (see
        :meth:`register_target`), so a custom target becomes selectable via
        ``transpile(dialect=target.name)`` as soon as any expander is registered
        against it — a user overriding one operator need not also declare the
        target separately.

        A *non-generic* ``(target, operator)`` entry additionally acts as a
        *join-rewrite override* (:meth:`has_override`) for operators with a
        built-in whole-query join rewrite (notably
        :class:`~giql.expressions.Intersects`, whose binned equi-join / DuckDB
        IEJoin transformers run before expansion), letting a per-target expander
        assume responsibility for that rewrite — :func:`giql.transpile.transpile`
        consults :meth:`has_override` and bypasses the built-in transformers when
        it holds.
        """
        self._expanders[(target, operator)] = _as_callable(expander)
        self.register_target(target)

    def register_target(self, target: Target) -> None:
        """Declare *target* so ``transpile(dialect=target.name)`` resolves it.

        The registry doubles as the target plugin hub: a custom
        :class:`~giql.targets.Target` registered here (directly, or as a side
        effect of :meth:`register`) is resolvable by name through
        :func:`giql.targets.resolve_target`. This is the path for a
        *capability-only* target — one that overrides no operators (so it never
        calls :meth:`register`) but declares a distinct
        :class:`~giql.targets.Capabilities` set / ``sqlglot_dialect``.

        The built-in :class:`~giql.targets.GenericTarget` is intentionally **not**
        declared: it is the portable fallback, selected only via ``dialect=None``
        (see :func:`giql.targets.resolve_target`), never by name. The skip is
        *by value* (``target == GenericTarget()``), so registering an expander
        against ``GenericTarget`` is a no-op here — but a *custom*
        :class:`~giql.targets.Target` subclass is declared under its own ``name``
        even if that name happens to be ``"generic"`` (it does not compare equal to
        the built-in). A later registration for the same name replaces an earlier
        one.
        """
        if target == GenericTarget():
            return
        self._targets_by_name[target.name] = target

    def target(self, name: str) -> Target | None:
        """Return the registered custom target named *name*, or ``None``.

        Consulted by :func:`giql.targets.resolve_target` for any dialect name
        that is not a built-in, so a registered custom target is selectable via
        ``transpile(dialect=name)``.
        """
        return self._targets_by_name.get(name)

    def resolve(self, target: Target, operator: type) -> ExpanderFn | None:
        """Return the expander for ``(target, operator)`` via the fallback chain.

        Tries the exact ``(target, op)`` entry, then the
        ``(GenericTarget(), op)`` fallback, then ``None``. ``None`` means no
        expander is registered; the :class:`ExpandOperators` pass treats that as
        an internal invariant violation and raises (there is no legacy emitter).

        A non-generic exact ``(target, op)`` entry additionally acts as a
        *join-rewrite override* for operators with a built-in whole-query join
        rewrite (notably :class:`~giql.expressions.Intersects`); resolution itself
        does not bypass the built-in binned / IEJoin transformers —
        :func:`giql.transpile.transpile` does, gated on :meth:`has_override` (see
        :meth:`register` and #141).
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

        Such an entry marks a target-specific override that supersedes built-in
        handling (e.g. taking responsibility for the whole-query join rewrite the
        built-in transformers would otherwise perform).
        :func:`giql.transpile.transpile` consults this method for
        ``(target, Intersects)`` and, when it holds, bypasses the built-in
        binned / IEJoin join transformers so the ``Intersects`` node flows to the
        registered expander instead (see #141).
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
        """Drop every registration — both expanders and declared targets.

        The bulk form of :meth:`unregister`, and the **public reset** for the
        process-wide :data:`REGISTRY` — e.g. a test fixture that saves and
        restores the registry around a body that registers custom expanders or
        targets.
        """
        self._expanders.clear()
        self._targets_by_name.clear()

    def snapshot(self) -> RegistrySnapshot:
        """Return a shallow copy of the current registrations and targets.

        The save half of a save/restore seam that supports test
        baseline-isolation: capture the baseline with this and hand it back to
        :meth:`restore` afterward, so the built-in registrations survive an
        isolating fixture that would otherwise :meth:`clear` them permanently.

        The returned :class:`RegistrySnapshot` is a fresh, opaque value bundling
        both the ``(target, operator)`` expander map and the declared-target map;
        mutating the registry afterward does not affect it, and two snapshots
        compare equal when the registry is in the same state.
        """
        return RegistrySnapshot(dict(self._expanders), dict(self._targets_by_name))

    def restore(self, snapshot: RegistrySnapshot) -> None:
        """Replace all registrations with those captured by :meth:`snapshot`.

        The restore half of the save/restore seam that supports test
        baseline-isolation. Drops every current entry and re-installs exactly the
        *snapshot* contents (expanders and declared targets), so a fixture can
        return the registry to a previously captured baseline regardless of what
        its body registered or cleared.
        """
        self._expanders.clear()
        self._expanders.update(snapshot.expanders)
        self._targets_by_name.clear()
        self._targets_by_name.update(snapshot.targets)

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
#: it at import time via :mod:`giql.expanders`; the pass rewrites every GIQL
#: operator node by resolving its expander here.
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


# The GIQL operator expression classes the pass inspects. Every one is migrated,
# so each resolves an expander through the registry (see module docstring).
# Imported eagerly at module scope, as ``canonicalizer.py`` imports the same
# classes: ``giql.expressions`` does not import this module, so there is no cycle.
_GIQL_OPERATORS: tuple[type, ...] = (
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
    """Replace every GIQL operator node with its registry expansion.

    Pass 3 of the normalization pipeline (epic #137). Runs after
    :func:`giql.canonicalizer.canonicalize_coordinates`. Every GIQL operator is
    migrated, so the pass dispatches *every* operator node: it resolves an expander
    for ``(target, operator type)`` through the registry's fallback chain and
    replaces the node with the AST that expander returns. A resolution miss is an
    internal invariant violation (a built-in operator always has at least a
    ``(generic, op)`` expander) and raises — there is no legacy ``*_sql`` fallback.

    The pass mutates *expression* in place for the node-local replacements, then
    applies any :data:`StatementFinalizer` an expander registered (via
    :meth:`ExpansionContext.add_statement_finalizer`) to the statement in
    registration order. A finalizer may return a *new* root, so callers must use
    the return value rather than assume in-place mutation.

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
        *expression* with each operator node replaced by its target-specific
        expansion, after any registered statement finalizers have run — the same
        object when no finalizer replaced the root, otherwise the finalized root.
    """
    reg = registry if registry is not None else REGISTRY
    operators = _GIQL_OPERATORS
    alias_seq = name_sequence(EXPAND_ALIAS_PREFIX)
    # Shared across every context this run builds: an expander that must rewrite
    # the enclosing statement appends a finalizer here, applied after the walk.
    finalizers: list[StatementFinalizer] = []

    # Collect first, then mutate: replacing nodes mid-walk is unsafe.
    pending: list[tuple[exp.Expression, ExpanderFn]] = []
    for node in expression.walk():
        if not isinstance(node, operators):
            continue
        fn = reg.resolve(target, type(node))
        if fn is None:
            # Every built-in operator registers at least a ``(generic, op)``
            # expander at import time, so a miss means the registry was cleared or
            # a custom target shadowed an operator without providing one — an
            # internal/config error, not user input. Fail loudly rather than
            # leaving an un-serializable GIQL node for the stock serializer (there
            # is no legacy ``*_sql`` fallback anymore).
            raise ValueError(
                f"No expander registered for {type(node).__name__} on target "
                f"{target.name!r}; the ExpandOperators pass cannot leave a GIQL "
                "operator node un-expanded."
            )
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
        ctx = ExpansionContext(
            node,
            resolution,
            target,
            tables,
            alias_seq,
            registry=reg,
            finalizers=finalizers,
        )
        replacement = fn(node, ctx)
        if not isinstance(replacement, exp.Expression):
            raise TypeError(
                f"expander {fn!r} for {type(node).__name__} returned "
                f"{type(replacement).__name__}, not exp.Expression"
            )
        if replacement is not node:
            node.replace(replacement)

    # Apply any query-level finalizers an expander registered, in registration
    # order, once every node-local replacement is in place. A finalizer receives
    # the current root and returns the (possibly new) root to thread forward.
    for finalize in finalizers:
        expression = finalize(expression)
        if not isinstance(expression, exp.Expression):
            # Mirror the node-local return guard: a finalizer that forgets to
            # return the root would otherwise make the pass silently return the
            # non-Expression far from the cause.
            raise TypeError(
                f"statement finalizer {finalize!r} returned "
                f"{type(expression).__name__}, not exp.Expression"
            )

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
