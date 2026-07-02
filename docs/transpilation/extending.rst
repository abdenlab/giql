Extending GIQL: custom targets and operators
=============================================

Every GIQL operator is expanded to standard SQL by a small **expander** function
selected from a process-wide registry, keyed by ``(target, operator)``. That
registry is a supported public extension point. Through it you can:

- **register a custom target** — a new SQL engine with its own capabilities and
  sqlglot output dialect — and select it with ``transpile(dialect="<name>")``; and
- **override an operator's expander** for a target, so that operator emits SQL
  tailored to that engine.

Both hooks live in the ``giql.expander`` and ``giql.targets`` modules, and the key
symbols are re-exported from the top-level package:

.. code-block:: python

   from giql import (
       transpile,
       register,        # decorator: register/override an operator expander
       REGISTRY,        # the process-wide plugin registry
       Target,          # base class for a custom target
       Capabilities,    # a target's feature set
       GenericTarget,   # the portable fallback target
       ExpansionContext,
   )


Targets and capabilities
------------------------

A :class:`~giql.targets.Target` is a first-class SQL engine. It carries a
``name``, the ``sqlglot_dialect`` used to serialize the final AST (``None``
selects sqlglot's portable default), and a :class:`~giql.targets.Capabilities`
set. Capabilities — not scattered ``if dialect == ...`` branches — drive the
portable emission choices:

- ``supports_lateral`` — correlated ``LATERAL`` joins (NEAREST uses a decorrelated
  window-function form when this is ``False``);
- ``supports_star_replace`` — ``SELECT * REPLACE (...)`` (used for coordinate
  canonicalization and the DISJOIN / NEAREST passthroughs; the portable
  ``SELECT * EXCEPT (...)`` form is emitted otherwise);
- ``supports_qualify`` — the ``QUALIFY`` clause;
- ``range_join_strategy`` — ``"binned"`` or ``"iejoin"`` for column-to-column
  INTERSECTS joins.

Define a custom target by subclassing :class:`~giql.targets.Target` as a frozen
dataclass. Give every field a default so the class is constructible with no
arguments (the ``@register`` decorator instantiates a target class for you):

.. code-block:: python

   from dataclasses import dataclass
   from giql import Target, Capabilities

   @dataclass(frozen=True)
   class PostgresTarget(Target):
       name: str = "postgres"
       sqlglot_dialect: str | None = "postgres"
       capabilities: Capabilities = Capabilities(
           supports_lateral=True,
           supports_star_replace=False,   # -> portable "* EXCEPT" canonicalization
           supports_qualify=True,
           range_join_strategy="binned",
       )


Registering and selecting a custom target
-----------------------------------------

Declare the target on the registry, then select it by name. A **capability-only**
target — one that reuses every built-in operator expander and differs only in its
capabilities and dialect — needs nothing more than
:meth:`~giql.expander.ExpanderRegistry.register_target`:

.. code-block:: python

   from giql import REGISTRY, transpile

   REGISTRY.register_target(PostgresTarget())

   sql = transpile(
       "SELECT * FROM peaks WHERE interval WITHIN 'chr1:1000-5000'",
       tables=["peaks"],
       dialect="postgres",
   )

``dialect="postgres"`` resolves the registered target; its
``supports_star_replace=False`` capability selects the portable projection form,
and its ``sqlglot_dialect="postgres"`` serializes the final SQL. (The built-in
names ``"duckdb"`` and ``"datafusion"`` and ``None`` for the generic target still
resolve as before; ``"generic"`` is *not* a selectable name — ``None`` is the
sole way to select the generic target.)


Overriding an operator's expander
---------------------------------

To tailor one operator to an engine, register an expander for
``(target, operator)`` with the :func:`~giql.expander.register` decorator. An
expander takes the operator node and its :class:`~giql.expander.ExpansionContext`
and returns the sqlglot AST that replaces the node:

.. code-block:: python

   from sqlglot import parse_one
   from giql import register, ExpansionContext
   from giql.expressions import Within

   @register(PostgresTarget, Within)
   def within_between(node, ctx: ExpansionContext):
       # ctx.resolution exposes the operator's resolved column operands, and
       # ctx.capabilities / ctx.target the active engine's capabilities.
       col = ctx.resolution.column("this")
       lo, hi = 1000, 5000  # (a real expander parses the node's range)
       return parse_one(f"{col.start} BETWEEN {lo} AND {hi}")

Registering an expander also **declares its target by name** as a side effect, so
``@register(PostgresTarget, Within)`` alone makes ``dialect="postgres"``
resolvable — a separate ``register_target`` call is only needed for a target that
overrides no operators.

Resolution follows a fallback chain: an exact ``(target, operator)`` expander wins,
otherwise the built-in ``(GenericTarget(), operator)`` expander runs. So you only
register the operators that genuinely differ for your engine; everything else
reuses the portable built-ins.

.. note::

   ``ctx.resolution.column("this").start`` (and ``.end`` / ``.chrom``) return
   **quoted** physical column names (e.g. ``"start"``, or ``a."start"`` when the
   operand is aliased), ready to splice into a SQL string that you ``parse_one``.
   Do not pass them to ``sqlglot.exp.column()``,
   which would quote them a second time. Building fragments as strings and parsing
   once — as the built-in expanders do — is the simplest correct approach.


The node-local boundary
-----------------------

An expander's **return value** is node-local: ``expand(node, ctx)`` returns the
one expression that replaces the operator node in place. It cannot *return* a
reshaped enclosing query. An expander may still restructure the query it sits in
as a side effect and then return the node unchanged — the built-in CLUSTER and
MERGE expanders do exactly this, rewriting their single-table ``SELECT`` in place.
What no expander can express is a rewrite that **adds or reshapes joins** across
relations: the DuckDB IEJoin plan for column-to-column INTERSECTS joins is handled
by a capability-gated pre-pass transformer, not an expander, because it restructures
the surrounding join. A general query-level expander seam for such join rewrites is
planned future work.


Undoing a registration
-----------------------

The registry exposes a teardown seam so a plugin or a test fixture can undo its
registrations rather than reaching into private state:

.. code-block:: python

   REGISTRY.unregister(PostgresTarget(), Within)   # drop one expander entry
   REGISTRY.clear()                                 # drop all entries and targets

   saved = REGISTRY.snapshot()                      # save/restore around a body
   try:
       ...                                          # register custom entries
   finally:
       REGISTRY.restore(saved)

See :func:`giql.transpile.transpile`, :func:`giql.expander.register`, and
:class:`giql.targets.Capabilities` in the :doc:`api-reference`.
