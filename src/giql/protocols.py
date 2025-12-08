"""Protocol definitions for GIQL.

This module defines protocols for type checking and interface compatibility.
"""

from typing import Any
from typing import Protocol
from typing import Sequence


class CursorLike(Protocol):
    """Protocol for DB-API 2.0 compatible cursors.

    Based on PEP 249: https://peps.python.org/pep-0249/

    This protocol defines the minimal interface required for database cursors
    that can be used with GIQL. All DB-API 2.0 compliant drivers (SQLite,
    PostgreSQL, MySQL, DuckDB) implement this interface.
    """

    @property
    def description(
        self,
    ) -> (
        Sequence[
            tuple[str, Any, Any | None, Any | None, Any | None, Any | None, Any | None]
        ]
        | None
    ):
        """Column descriptions.

        A sequence of 7-tuples describing each column:
        (name, type_code, display_size, internal_size, precision, scale, null_ok)

        Only 'name' is required; other values may be None.
        Returns None if no operation has been performed yet.
        """
        ...

    @property
    def rowcount(self) -> int:
        """Number of rows affected by last operation.

        Returns -1 if no operation has been performed or if the count
        cannot be determined.
        """
        ...

    def fetchone(self) -> tuple[Any, ...] | None:
        """Fetch the next row of a query result set.

        Returns a tuple representing the next row, or None when no more
        rows are available.
        """
        ...

    def fetchmany(self, size: int = 1) -> list[tuple[Any, ...]]:
        """Fetch the next set of rows of a query result set.

        Returns a list of tuples. An empty list is returned when no more
        rows are available.

        :param size:
            Number of rows to fetch (default: 1)
        """
        ...

    def fetchall(self) -> list[tuple[Any, ...]]:
        """Fetch all remaining rows of a query result set.

        Returns a list of tuples. An empty list is returned when no rows
        are available.
        """
        ...

    def close(self) -> None:
        """Close the cursor.

        Makes the cursor unusable for further operations.
        """
        ...
