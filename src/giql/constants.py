"""Default constants for GIQL.

This module defines default column names and other constants used throughout GIQL.
"""

# Default genomic column names
DEFAULT_CHROM_COL = "chrom"
DEFAULT_START_COL = "start"
DEFAULT_END_COL = "end"
DEFAULT_STRAND_COL = "strand"
DEFAULT_GENOMIC_COL = "interval"

# Canonical column tuple a CTE or subquery input to DISJOIN is contractually
# assumed to expose (0-based half-open).
CANONICAL_DEFAULT_COLS: tuple[str, str, str] = (
    DEFAULT_CHROM_COL,
    DEFAULT_START_COL,
    DEFAULT_END_COL,
)

# Default bin size for INTERSECTS binned equi-join optimization
DEFAULT_BIN_SIZE = 10_000
