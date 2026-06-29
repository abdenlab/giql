"""Custom AST expression nodes for genomic operations.

This module defines custom SQLGlot expression nodes for GIQL spatial operators.

It also declares the per-operator :class:`SlotSpec` descriptors that the
``ResolveOperatorRefs`` normalization pass (:mod:`giql.resolver`) consumes. Each
operator class carries a ``GIQL_SLOTS`` tuple naming the slots the operator
resolves (``this`` / ``reference`` / ``expression`` / ``ranges``) together with
the AST shapes that slot accepts. Keeping the acceptance rules declarative on the
expression class — rather than duplicated across resolver code — is the structural
change epic #114 is built around.
"""

from dataclasses import dataclass
from typing import Literal

from sqlglot import exp
from sqlglot.errors import ParseError

# The shapes an operator slot may accept. The first three form the
# "registered table / CTE / subquery" trichotomy that ``sqlglot``'s scope
# machinery distinguishes natively and that :class:`giql.resolver.ResolvedRef`
# models; the remaining three describe the column / literal / implicit-outer
# forms used by NEAREST, DISTANCE, and the spatial predicates.
SlotShape = Literal[
    "registered_table",
    "cte",
    "subquery",
    "literal_range",
    "column",
    "implicit_outer",
]

# The subset of slot shapes that resolve to a table-shaped relation and are
# therefore representable as a :class:`giql.resolver.ResolvedRef`. A slot whose
# accepted shapes are a subset of these is a "reference slot" that the
# ``ResolveOperatorRefs`` pass resolves and attaches metadata for.
TABLE_SHAPES: frozenset[SlotShape] = frozenset({"registered_table", "cte", "subquery"})


@dataclass(frozen=True, slots=True)
class SlotSpec:
    """Declarative description of one resolvable operator slot.

    Parameters
    ----------
    arg : str
        The :attr:`arg_types` key naming the slot on the operator expression
        (e.g. ``"this"``, ``"reference"``, ``"expression"``, ``"ranges"``).
    accepts : frozenset[SlotShape]
        The AST shapes the slot accepts.
    required : bool
        Whether the slot must be present on a well-formed operator.
    is_sequence : bool
        Whether the slot holds a list of expressions (e.g. the ``ranges`` slot
        of :class:`SpatialSetPredicate`) rather than a single expression.
    """

    arg: str
    accepts: frozenset[SlotShape]
    required: bool = False
    is_sequence: bool = False

    @property
    def is_ref_slot(self) -> bool:
        """Whether this slot resolves to a table-shaped :class:`ResolvedRef`.

        A reference slot is one whose accepted shapes are all drawn from the
        registered-table / CTE / subquery trichotomy. These are the slots for
        which the ``ResolveOperatorRefs`` pass resolves and attaches a
        :class:`giql.resolver.ResolvedRef`; column / literal / implicit
        slots are declared but their resolution is deferred to later epic #114
        migration steps.
        """
        return bool(self.accepts) and self.accepts <= TABLE_SHAPES


class GenomicRange(exp.Expression):
    """Represents a parsed genomic range.

    Examples:
        'chr1:1000-2000'
        'chr1:[1000,2000)'
        'chr1:[1001,2000]'
    """

    arg_types = {
        "chromosome": True,
        "start": True,
        "end": True,
        "strand": False,
        "coord_system": False,
    }


class SpatialPredicate(exp.Binary):
    """Base class for spatial predicates."""

    pass


#: Opt the spatial predicates, DISTANCE, and NEAREST into the
#: CanonicalizeCoordinates pass (epic #114 step 8, issue #123). With this flag
#: set, pass 2 canonicalizes each operator's interval operands: a non-table column
#: operand (DISTANCE / predicate operand, NEAREST's column / implicit-outer
#: reference) has its resolution metadata rewritten in place to canonical 0-based
#: half-open arithmetic, and a registered-table reference slot (NEAREST's target)
#: is wrapped in a ``__giql_canon_*`` CTE. The emitter then consumes already-
#: canonical fragments with no in-emitter canonicalization. Identity (0-based
#: half-open) operands are left untouched and the emitted SQL stays byte-identical.
_CANONICALIZE = True

#: Default opt-out from the ``ExpandOperators`` pass (epic #137, step 2). The
#: per-operator ``GIQL_EXPAND`` flag mirrors ``GIQL_CANONICALIZE``: an operator
#: takes the new AST-expansion path only when it sets ``GIQL_EXPAND = True`` *and*
#: an expander is registered for it; otherwise the legacy ``*_sql`` emitter runs.
#: This is the opt-out default: an operator inherits it and stays on the legacy
#: emitter until its migration step flips the flag to ``True`` alongside its
#: registered expander. Operators already migrated override it on their own class.
_EXPAND = False


class Intersects(SpatialPredicate):
    """INTERSECTS spatial predicate.

    Example: column INTERSECTS 'chr1:1000-2000'
    """

    GIQL_CANONICALIZE = _CANONICALIZE
    #: Migrated to the ExpandOperators registry (#141). A literal-range or
    #: residual column-to-column INTERSECTS *predicate* expands through
    #: ``giql.expanders.intersects``; a column-to-column INTERSECTS *join* is
    #: consumed by the capability-gated binned / IEJoin pre-pass transformers
    #: before this pass runs, so the predicate expander never sees it.
    GIQL_EXPAND = True

    GIQL_SLOTS = (
        SlotSpec("this", frozenset({"column"}), required=True),
        SlotSpec("expression", frozenset({"literal_range", "column"}), required=True),
    )


class Contains(SpatialPredicate):
    """CONTAINS spatial predicate.

    Example: column CONTAINS 'chr1:1500'
    """

    GIQL_CANONICALIZE = _CANONICALIZE
    #: Migrated to the ExpandOperators registry (#141); expands through
    #: ``giql.expanders.intersects``.
    GIQL_EXPAND = True

    GIQL_SLOTS = (
        SlotSpec("this", frozenset({"column"}), required=True),
        SlotSpec("expression", frozenset({"literal_range", "column"}), required=True),
    )


class Within(SpatialPredicate):
    """WITHIN spatial predicate.

    Example: column WITHIN 'chr1:1000-5000'
    """

    GIQL_CANONICALIZE = _CANONICALIZE
    #: Migrated to the ExpandOperators registry (#141); expands through
    #: ``giql.expanders.intersects``.
    GIQL_EXPAND = True

    GIQL_SLOTS = (
        SlotSpec("this", frozenset({"column"}), required=True),
        SlotSpec("expression", frozenset({"literal_range", "column"}), required=True),
    )


class SpatialSetPredicate(exp.Expression):
    """Spatial predicates with set quantifiers.

    Examples:
        column INTERSECTS ANY('chr1:1000-2000', 'chr1:5000-6000')
        column CONTAINS ALL('chr1:1500', 'chr1:1600')
    """

    arg_types = {
        "this": True,
        "operator": True,
        "quantifier": True,
        "ranges": True,
    }

    GIQL_CANONICALIZE = _CANONICALIZE
    #: Migrated to the ExpandOperators registry (#141); expands through
    #: ``giql.expanders.intersects``.
    GIQL_EXPAND = True

    GIQL_SLOTS = (
        SlotSpec("this", frozenset({"column"}), required=True),
        SlotSpec(
            "ranges",
            frozenset({"literal_range"}),
            required=True,
            is_sequence=True,
        ),
    )


def _split_named_and_positional(args):
    """Separate named parameters (:= and =>) from positional arguments."""
    kwargs = {}
    positional_args = []
    for arg in args:
        if isinstance(arg, (exp.PropertyEQ, exp.Kwarg)):
            param_name = arg.this.name if hasattr(arg.this, "name") else str(arg.this)
            kwargs[param_name.lower()] = arg.expression
        else:
            positional_args.append(arg)
    return kwargs, positional_args


class GIQLCluster(exp.Func):
    """CLUSTER window function for assigning cluster IDs to overlapping intervals.

    Implicitly partitions by chromosome and orders by start position.

    The optional ``predicate`` argument is a boolean expression evaluated
    between each interval and its sorted predecessor; intervals are only kept
    in the same cluster when they are adjacent *and* the predicate holds. Bare
    columns resolve to the current interval; the predecessor's value of a
    column is referenced with ``PREV(column)``.

    Examples:
        CLUSTER(interval)
        CLUSTER(interval, 1000)
        CLUSTER(interval, stranded := true)
        CLUSTER(interval, 1000, stranded := true)
        CLUSTER(interval, predicate := depth = PREV(depth))
    """

    arg_types = {
        "this": True,  # genomic column
        "distance": False,  # maximum distance between features
        "stranded": False,  # strand-specific clustering
        "predicate": False,  # pairwise boolean gate (current row vs PREV(col))
    }

    # Inert today: the CLUSTER/MERGE transformers rewrite these nodes before the
    # ExpandOperators pass runs, so the pass never sees a GIQLCluster to dispatch
    # and this flag is not a live opt-in. It is forward-looking for #144, which
    # migrates these operators onto the expander registry.
    GIQL_EXPAND = _EXPAND

    @classmethod
    def from_arg_list(cls, args):
        kwargs, positional_args = _split_named_and_positional(args)
        if len(positional_args) > 0:
            kwargs["this"] = positional_args[0]
        if len(positional_args) > 1:
            kwargs["distance"] = positional_args[1]
        if "this" not in kwargs:
            raise ParseError(
                "CLUSTER requires a genomic interval column as its first argument."
            )
        return cls(**kwargs)


class GIQLMerge(exp.Func):
    """MERGE aggregate function for combining overlapping intervals.

    Merges overlapping or bookended intervals into single intervals.
    Built on top of CLUSTER operation.

    The optional ``predicate`` argument gates merging on a pairwise boolean
    expression between each interval and its sorted predecessor (see
    :class:`GIQLCluster`); ``PREV(column)`` references the predecessor's value
    of a column. When the predicate tests equality of a value this yields a
    run-length encoding of the input interval sequence.

    Examples:
        MERGE(interval)
        MERGE(interval, 1000)
        MERGE(interval, stranded := true)
        MERGE(interval, predicate := depth = PREV(depth))
        MERGE(interval, predicate := strand = PREV(strand) AND name = PREV(name))
    """

    arg_types = {
        "this": True,  # genomic column
        "distance": False,  # maximum distance between features
        "stranded": False,  # strand-specific merging
        "predicate": False,  # pairwise boolean gate (current row vs PREV(col))
    }

    # Inert today: the CLUSTER/MERGE transformers rewrite these nodes before the
    # ExpandOperators pass runs, so the pass never sees a GIQLMerge to dispatch
    # and this flag is not a live opt-in. It is forward-looking for #144, which
    # migrates these operators onto the expander registry.
    GIQL_EXPAND = _EXPAND

    @classmethod
    def from_arg_list(cls, args):
        kwargs, positional_args = _split_named_and_positional(args)
        if len(positional_args) > 0:
            kwargs["this"] = positional_args[0]
        if len(positional_args) > 1:
            kwargs["distance"] = positional_args[1]
        if "this" not in kwargs:
            raise ParseError(
                "MERGE requires a genomic interval column as its first argument."
            )
        return cls(**kwargs)


class GIQLDistance(exp.Func):
    """DISTANCE function for calculating genomic distances between intervals.

    Generates SQL CASE expression that computes distance between two genomic
    intervals, with optional strand-specific and signed (directional) modes.

    Examples:
        DISTANCE(a.interval, b.interval)
        DISTANCE(a.interval, 'chr1:1000-2000')
        DISTANCE(a.interval, b.interval, stranded := true)
        DISTANCE(a.interval, b.interval, signed := true)
        DISTANCE(a.interval, b.interval, stranded := true, signed := true)
    """

    arg_types = {
        "this": True,  # Required: interval_a (column ref or literal range)
        "expression": True,  # Required: interval_b (column ref or literal range)
        "stranded": False,  # Optional: boolean for strand-specific distance
        "signed": False,  # Optional: boolean for directional distance
    }

    GIQL_CANONICALIZE = _CANONICALIZE
    # Migrated to the registry's AST-expansion path (epic #137, issue #140): the
    # generic expander in giql.expanders.distance builds the CASE; the legacy
    # giqldistance_sql emitter is gone.
    GIQL_EXPAND = True

    GIQL_SLOTS = (
        SlotSpec("this", frozenset({"column"}), required=True),
        SlotSpec("expression", frozenset({"column"}), required=True),
    )

    @classmethod
    def from_arg_list(cls, args):
        kwargs, positional_args = _split_named_and_positional(args)
        if len(positional_args) >= 1:
            kwargs["this"] = positional_args[0]
        if len(positional_args) >= 2:
            kwargs["expression"] = positional_args[1]
        return cls(**kwargs)


class GIQLNearest(exp.Func):
    """NEAREST function for finding k-nearest genomic features.

    Generates SQL for k-nearest neighbor queries using LATERAL joins
    (PostgreSQL/DuckDB) or window functions (SQLite).

    Examples:
        NEAREST(genes, k := 3)
        NEAREST(genes, reference := peaks.interval, k := 5)
        NEAREST(genes, reference := 'chr1:1000-2000', k := 3)
        NEAREST(genes, k := 5, max_distance := 100000, stranded := true)
    """

    arg_types = {
        "this": True,  # Required: target table name
        "reference": False,  # Optional: position reference (column or literal)
        "k": False,  # Optional: number of neighbors (default=1)
        "max_distance": False,  # Optional: distance threshold
        "stranded": False,  # Optional: strand-specific search
        "signed": False,  # Optional: directional distance
    }

    #: Opt NEAREST into the CanonicalizeCoordinates pass (epic #114 step 8, issue
    #: #123). Pass 2 wraps a non-canonical *target* (the ``this`` registered-table
    #: ref slot) in a ``__giql_canon_*`` CTE — so the emitter reads canonical
    #: target columns and de-canonicalizes its ``*`` row passthrough back to the
    #: declared encoding — and canonicalizes a non-table ``column`` /
    #: ``implicit_outer`` reference's metadata in place. Identity (0-based
    #: half-open) operands are left untouched and the emitted SQL stays
    #: byte-identical.
    GIQL_CANONICALIZE = _CANONICALIZE
    #: Migrated to the ExpandOperators pass (epic #137, issue #142): NEAREST is
    #: expanded by ``giql.expanders.nearest`` — the portable correlated LATERAL
    #: subquery where ``supports_lateral`` holds, a decorrelated window-function
    #: form otherwise. The legacy ``giqlnearest_sql`` emitter has been removed.
    GIQL_EXPAND = True

    GIQL_SLOTS = (
        SlotSpec("this", frozenset({"registered_table"}), required=True),
        SlotSpec(
            "reference",
            frozenset({"literal_range", "column", "implicit_outer"}),
        ),
    )

    @classmethod
    def from_arg_list(cls, args):
        kwargs, positional_args = _split_named_and_positional(args)
        if len(positional_args) >= 1:
            kwargs["this"] = positional_args[0]
        if "this" not in kwargs:
            raise ParseError("NEAREST requires a target table as its first argument.")
        return cls(**kwargs)


class GIQLDisjoin(exp.Func):
    """DISJOIN table function for splitting intervals at reference breakpoints.

    Generates SQL that cuts each target interval at every reference breakpoint
    strictly interior to it, so each resulting sub-interval is fully contained
    by every reference interval it overlaps. The target row passes through
    intact and the sub-interval is appended as ``disjoin_chrom`` /
    ``disjoin_start`` / ``disjoin_end``. When ``reference`` is omitted it
    defaults to the target set.

    Examples:
        DISJOIN(features)
        DISJOIN(features, reference := other)
        DISJOIN(features, reference => other)
        DISJOIN(features, reference := (SELECT ...))
    """

    arg_types = {
        "this": True,  # Required: target table name
        "reference": False,  # Optional: reference table/CTE name or subquery
    }

    #: Opt DISJOIN into the CanonicalizeCoordinates pass (epic #114 step 7,
    #: issue #122). With this flag set, pass 2 wraps every non-canonical
    #: interval-bearing operand in a canonical ``__giql_canon_*`` CTE and
    #: rewrites the slot to point at it, so the emitter consumes already-canonical
    #: 0-based half-open columns instead of canonicalizing inline. Identity
    #: (0-based half-open) operands are left unwrapped and the emitted SQL stays
    #: byte-identical.
    GIQL_CANONICALIZE = True
    #: Opt DISJOIN into the ExpandOperators pass (epic #137, issue #143). DISJOIN
    #: is migrated to a registered expander (``giql.expanders.disjoin``), so the
    #: pass replaces the node with the expander's AST and the legacy
    #: ``giqldisjoin_sql`` emitter is removed.
    GIQL_EXPAND = True

    GIQL_SLOTS = (
        SlotSpec("this", frozenset({"registered_table"}), required=True),
        SlotSpec(
            "reference",
            frozenset({"registered_table", "cte", "subquery"}),
        ),
    )

    @classmethod
    def from_arg_list(cls, args):
        kwargs, positional_args = _split_named_and_positional(args)
        unknown = set(kwargs) - set(cls.arg_types)
        if unknown:
            raise ParseError(
                f"DISJOIN got unexpected named argument(s): "
                f"{', '.join(sorted(unknown))}. Valid arguments: reference."
            )
        if len(positional_args) > 1:
            raise ParseError(
                "DISJOIN accepts at most one positional argument (the target "
                f"table); got {len(positional_args)}. Pass the reference set "
                "as 'reference := ...'."
            )
        if len(positional_args) >= 1:
            kwargs["this"] = positional_args[0]
        if "this" not in kwargs:
            raise ParseError("DISJOIN requires a target table as its first argument.")
        return cls(**kwargs)
