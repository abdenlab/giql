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

#: The reserved column name for CLUSTER's synthesized per-row "new cluster" flag
#: (``CASE ... END AS __giql_is_new_cluster``), consumed by the outer
#: ``SUM(...) OVER (...)`` window that assigns cluster ids. The reserved ``__giql_``
#: prefix keeps it clear of user columns: a user column named ``is_new_cluster``
#: coexisting with ``SELECT *`` otherwise mis-binds the SUM (#161). A single source
#: of truth so the alias and its ``SUM`` consumer cannot drift apart.
CLUSTER_FLAG_COL = "__giql_is_new_cluster"

#: The reserved column name for the cluster id MERGE synthesizes
#: (``CLUSTER(...) AS __giql_cluster_id``) and then groups by. Reserved-prefixed for
#: the same reason as :data:`CLUSTER_FLAG_COL`; a single source of truth so the alias
#: and its ``GROUP BY`` consumer cannot drift apart.
CLUSTER_ID_COL = "__giql_cluster_id"

#: The reserved alias prefix CLUSTER materializes a star's aliased sibling items
#: under, inside the ``__giql_lag_calc`` subquery, so the outer star can EXCEPT them
#: (and the outer projection reference them) without colliding with a same-named base
#: column or another sibling. Reserved-prefixed for the same reason as
#: :data:`CLUSTER_FLAG_COL`; ``SELECT *, expr AS score`` would otherwise strip the base
#: ``score`` from the star, and ``SELECT *, 1 AS x, 2 AS x`` would duplicate the EXCEPT
#: entry (#190).
CLUSTER_SIBLING_PREFIX = "__giql_sibling_"
