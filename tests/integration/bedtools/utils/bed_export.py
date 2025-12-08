"""DuckDB loading utilities for genomic intervals.

This module provides functions for loading genomic intervals into DuckDB tables.
"""

from typing import List
from typing import Tuple


def load_intervals(
    conn,
    table_name: str,
    intervals: List[Tuple[str, int, int, str | None, int | None, str | None]],
):
    """Load intervals into DuckDB table.

    Args:
        conn: DuckDB connection
        table_name: Name of table to create
        intervals: List of (chrom, start, end, name, score, strand) tuples
                   where name, score, and strand can be None

    Note:
        Creates a new table with GIQL's default column names for genomic data:
        chromosome, start_pos, end_pos, name, score, strand
    """
    # Create table with GIQL's default column names
    conn.execute(f"""
        CREATE TABLE {table_name} (
            chromosome VARCHAR,
            start_pos INTEGER,
            end_pos INTEGER,
            name VARCHAR,
            score INTEGER,
            strand VARCHAR
        )
    """)

    # Insert intervals
    conn.executemany(f"INSERT INTO {table_name} VALUES (?,?,?,?,?,?)", intervals)
