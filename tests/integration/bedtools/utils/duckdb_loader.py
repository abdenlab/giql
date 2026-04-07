"""Utilities for loading interval data into DuckDB tables."""


def load_intervals(conn, table_name: str, intervals: list[tuple]) -> None:
    """Load interval tuples into a DuckDB table.

    Creates a table with GIQL default column names (chrom, start, end,
    name, score, strand) and inserts the provided intervals.

    Args:
        conn: DuckDB connection
        table_name: Name of the table to create.  Must be a simple
            identifier -- this is test-only code with controlled inputs.
        intervals: List of 6-element tuples
            (chrom, start, end, name, score, strand)
    """
    conn.execute(f"""
        CREATE TABLE {table_name} (
            chrom VARCHAR,
            "start" INTEGER,
            "end" INTEGER,
            name VARCHAR,
            score INTEGER,
            strand VARCHAR
        )
    """)
    conn.executemany(f"INSERT INTO {table_name} VALUES (?,?,?,?,?,?)", intervals)
