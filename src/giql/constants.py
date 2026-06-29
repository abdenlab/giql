"""Default constants for GIQL.

This module defines default column names and other constants used throughout GIQL.
"""

# Default genomic column names
DEFAULT_CHROM_COL = "chrom"
DEFAULT_START_COL = "start"
DEFAULT_END_COL = "end"
DEFAULT_STRAND_COL = "strand"
DEFAULT_GENOMIC_COL = "interval"

# Default bin size for INTERSECTS binned equi-join optimization
DEFAULT_BIN_SIZE = 10_000

#: The reserved alias prefix naming DISJOIN's internal CTEs (``__giql_dj_ref`` /
#: ``__giql_dj_tgt`` / ``__giql_dj_cuts`` / ...). Follows the same
#: reserved-prefix scheme as :data:`giql.canonicalizer.CANON_PREFIX`; the leading
#: double underscore keeps
#: the namespace clear of user identifiers. A single source of truth so the
#: resolver's reserved-prefix guard and the DISJOIN expander's CTE names and
#: guard cannot drift apart.
DJ_PREFIX = "__giql_dj_"
