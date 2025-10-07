"""Custom AST expression nodes for genomic operations.

This module defines custom SQLGlot expression nodes for GIQL spatial operators.
"""

from sqlglot import exp


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


class Intersects(SpatialPredicate):
    """INTERSECTS spatial predicate.

    Example: column INTERSECTS 'chr1:1000-2000'
    """

    pass


class Contains(SpatialPredicate):
    """CONTAINS spatial predicate.

    Example: column CONTAINS 'chr1:1500'
    """

    pass


class Within(SpatialPredicate):
    """WITHIN spatial predicate.

    Example: column WITHIN 'chr1:1000-5000'
    """

    pass


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
