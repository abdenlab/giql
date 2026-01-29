"""Table configuration for GIQL transpilation.

This module defines the Table dataclass for configuring genomic table schemas.
"""

from dataclasses import dataclass
from typing import Literal

from giql.constants import DEFAULT_CHROM_COL
from giql.constants import DEFAULT_END_COL
from giql.constants import DEFAULT_GENOMIC_COL
from giql.constants import DEFAULT_START_COL
from giql.constants import DEFAULT_STRAND_COL


@dataclass
class Table:
    """Genomic table configuration for transpilation.

    This class defines how genomic intervals are stored in a database table,
    mapping a pseudo-column name (genomic_col) to the physical columns that
    store chromosome, start, end, and optionally strand information.

    Parameters
    ----------
    genomic_col : str
        The pseudo-column name used in GIQL queries to reference the genomic
        interval (default: "interval").
    chrom_col : str
        The physical column name storing chromosome/contig (default: "chromosome").
    start_col : str
        The physical column name storing interval start position
        (default: "start_pos").
    end_col : str
        The physical column name storing interval end position
        (default: "end_pos").
    strand_col : str | None
        The physical column name storing strand information, or None if the
        table has no strand column (default: "strand").
    coordinate_system : Literal["0based", "1based"]
        The coordinate system used for positions (default: "0based").
    interval_type : Literal["half_open", "closed"]
        The interval endpoint convention (default: "half_open").

    Examples
    --------
    Using default column names (via transpile)::

        sql = transpile(query, tables=["peaks"])

    Using custom column names::

        sql = transpile(
            query,
            tables={
                "variants": Table(
                    genomic_col="position",
                    chrom_col="chr",
                    start_col="pos_start",
                    end_col="pos_end",
                    strand_col=None,  # No strand column
                    coordinate_system="1based",
                    interval_type="closed",
                )
            }
        )
    """

    genomic_col: str = DEFAULT_GENOMIC_COL
    chrom_col: str = DEFAULT_CHROM_COL
    start_col: str = DEFAULT_START_COL
    end_col: str = DEFAULT_END_COL
    strand_col: str | None = DEFAULT_STRAND_COL
    coordinate_system: Literal["0based", "1based"] = "0based"
    interval_type: Literal["half_open", "closed"] = "half_open"

    def __post_init__(self) -> None:
        """Validate field values after initialization."""
        if self.coordinate_system not in ("0based", "1based"):
            raise ValueError(
                f"coordinate_system must be '0based' or '1based', "
                f"got {self.coordinate_system!r}"
            )
        if self.interval_type not in ("half_open", "closed"):
            raise ValueError(
                f"interval_type must be 'half_open' or 'closed', "
                f"got {self.interval_type!r}"
            )


class Tables:
    """Container for Table configurations.

    Provides lookup of Table objects by name for use during transpilation.
    """

    def __init__(self) -> None:
        self._tables: dict[str, Table] = {}

    def register(self, name: str, table: Table) -> None:
        """Register a table configuration.

        Parameters
        ----------
        name : str
            The table name to register.
        table : Table
            Table configuration to register.
        """
        self._tables[name] = table

    def get(self, name: str) -> Table | None:
        """Get a table configuration by name.

        Parameters
        ----------
        name : str
            Table name to look up.

        Returns
        -------
        Table | None
            Table configuration if found, None otherwise.
        """
        return self._tables.get(name)

    def __contains__(self, name: str) -> bool:
        return name in self._tables

    def __iter__(self):
        return iter(self._tables.values())
