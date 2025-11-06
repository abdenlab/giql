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


class GIQLCluster(exp.Func):
    """CLUSTER window function for assigning cluster IDs to overlapping intervals.

    Implicitly partitions by chromosome and orders by start position.

    Examples:
        CLUSTER(position)
        CLUSTER(position, 1000)
        CLUSTER(position, stranded=true)
        CLUSTER(position, 1000, stranded=true)
    """

    arg_types = {
        "this": True,  # genomic column
        "distance": False,  # maximum distance between features
        "stranded": False,  # strand-specific clustering
    }

    @classmethod
    def from_arg_list(cls, args):
        """Parse argument list, handling named parameters.

        :param args:
            List of arguments from parser
        :return:
            GIQLCluster instance with properly mapped arguments
        """
        kwargs = {}
        positional_args = []

        # Separate named (EQ) and positional arguments
        for arg in args:
            if isinstance(arg, exp.EQ):
                # Named parameter: extract name and value
                param_name = (
                    arg.this.name if isinstance(arg.this, exp.Column) else str(arg.this)
                )
                kwargs[param_name.lower()] = arg.expression
            else:
                positional_args.append(arg)

        # Map positional arguments
        if len(positional_args) > 0:
            kwargs["this"] = positional_args[0]
        if len(positional_args) > 1:
            kwargs["distance"] = positional_args[1]

        return cls(**kwargs)


class GIQLMerge(exp.Func):
    """MERGE aggregate function for combining overlapping intervals.

    Merges overlapping or bookended intervals into single intervals.
    Built on top of CLUSTER operation.

    Examples:
        MERGE(position)
        MERGE(position, 1000)
        MERGE(position, stranded=true)
    """

    arg_types = {
        "this": True,  # genomic column
        "distance": False,  # maximum distance between features
        "stranded": False,  # strand-specific merging
    }

    @classmethod
    def from_arg_list(cls, args):
        """Parse argument list, handling named parameters.

        :param args: List of arguments from parser
        :return: GIQLMerge instance with properly mapped arguments
        """
        kwargs = {}
        positional_args = []

        # Separate named (EQ) and positional arguments
        for arg in args:
            if isinstance(arg, exp.EQ):
                # Named parameter: extract name and value
                param_name = (
                    arg.this.name if isinstance(arg.this, exp.Column) else str(arg.this)
                )
                kwargs[param_name.lower()] = arg.expression
            else:
                positional_args.append(arg)

        # Map positional arguments
        if len(positional_args) > 0:
            kwargs["this"] = positional_args[0]
        if len(positional_args) > 1:
            kwargs["distance"] = positional_args[1]

        return cls(**kwargs)


class GIQLDistance(exp.Func):
    """DISTANCE function for calculating genomic distances between intervals.

    Generates SQL CASE expression that computes distance between two genomic
    intervals, with optional strand-specific and signed (directional) modes.

    Examples:
        DISTANCE(a.position, b.position)
        DISTANCE(a.position, 'chr1:1000-2000')
        DISTANCE(a.position, b.position, stranded=true)
        DISTANCE(a.position, b.position, signed=true)
        DISTANCE(a.position, b.position, stranded=true, signed=true)
    """

    arg_types = {
        "this": True,  # Required: interval_a (column ref or literal range)
        "expression": True,  # Required: interval_b (column ref or literal range)
        "stranded": False,  # Optional: boolean for strand-specific distance
        "signed": False,  # Optional: boolean for directional distance
    }

    @classmethod
    def from_arg_list(cls, args):
        """Parse argument list, handling named parameters.

        :param args:
            List of arguments from parser
        :return:
            GIQLDistance instance with properly mapped arguments
        """
        kwargs = {}
        positional_args = []

        # Separate named (EQ) and positional arguments
        for arg in args:
            if isinstance(arg, exp.EQ):
                # Named parameter: extract name and value
                param_name = (
                    arg.this.name if isinstance(arg.this, exp.Column) else str(arg.this)
                )
                kwargs[param_name.lower()] = arg.expression
            else:
                positional_args.append(arg)

        # Map positional arguments
        if len(positional_args) >= 1:
            kwargs["this"] = positional_args[0]
        if len(positional_args) >= 2:
            kwargs["expression"] = positional_args[1]

        return cls(**kwargs)
