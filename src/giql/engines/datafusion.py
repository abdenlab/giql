"""DataFusion execution engine for GIQL queries.

Thin wrapper composing GIQL's transpile() with Apache DataFusion's
SessionContext to provide an integrated query engine over Parquet files
and Arrow tables.
"""

from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING

from giql.table import Table
from giql.transpile import transpile

if TYPE_CHECKING:
    import pyarrow


class DataFusionEngine:
    """Execute GIQL queries via Apache DataFusion.

    Uses GIQL's transpiler to convert queries to SQL, then executes
    them against DataFusion's SessionContext. Data can be registered
    from Parquet files or in-memory Arrow tables.

    Parameters
    ----------
    None

    Examples
    --------
    ::

        engine = DataFusionEngine()
        engine.register_parquet("peaks", "peaks.parquet")
        result = engine.query(
            "SELECT * FROM peaks WHERE interval INTERSECTS 'chr1:1000-2000'"
        )
    """

    def __init__(self) -> None:
        try:
            import datafusion
        except ImportError:
            raise ImportError(
                "datafusion is required for DataFusionEngine. "
                "Install it with: pip install 'giql[datafusion]'"
            ) from None

        self._ctx = datafusion.SessionContext()
        self._tables: dict[str, Table] = {}

    def register_parquet(
        self,
        name: str,
        path: str | Path,
        table: Table | None = None,
    ) -> None:
        """Register a Parquet file as a named table.

        Parameters
        ----------
        name : str
            Table name for use in queries.
        path : str or Path
            Path to the Parquet file.
        table : Table or None
            Optional GIQL Table config with custom column mappings.
            If None, default column mappings are used.
        """
        self._ctx.register_parquet(name, str(path))
        self._tables[name] = table if table is not None else Table(name)

    def register_arrow(
        self,
        name: str,
        data: pyarrow.RecordBatch | pyarrow.Table,
        table: Table | None = None,
    ) -> None:
        """Register an Arrow table or batch as a named table.

        Parameters
        ----------
        name : str
            Table name for use in queries.
        data : pyarrow.RecordBatch or pyarrow.Table
            Arrow data to register.
        table : Table or None
            Optional GIQL Table config with custom column mappings.
            If None, default column mappings are used.
        """
        import pyarrow

        if isinstance(data, pyarrow.Table):
            batches = data.to_batches()
        else:
            batches = [data]

        self._ctx.register_record_batches(name, [batches])
        self._tables[name] = table if table is not None else Table(name)

    def query(self, giql_query: str) -> pyarrow.Table:
        """Execute a GIQL query and return results as an Arrow table.

        Parameters
        ----------
        giql_query : str
            GIQL query string.

        Returns
        -------
        pyarrow.Table
            Query results.
        """
        sql = transpile(giql_query, tables=list(self._tables.values()))
        return self._ctx.sql(sql).to_arrow_table()

    def sql(self, raw_sql: str) -> pyarrow.Table:
        """Execute raw SQL directly against DataFusion.

        Useful for setup statements or debugging transpiled output.

        Parameters
        ----------
        raw_sql : str
            SQL query string.

        Returns
        -------
        pyarrow.Table
            Query results.
        """
        return self._ctx.sql(raw_sql).to_arrow_table()

    @property
    def context(self):
        """The underlying DataFusion SessionContext.

        Exposed for advanced use cases and benchmarking.
        """
        return self._ctx
