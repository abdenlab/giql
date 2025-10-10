"""Multi-backend query engine for GIQL.

This module provides the main query engine that supports multiple SQL databases
through transpilation of GIQL syntax to standard SQL.
"""

from typing import Literal

import pandas as pd
from sqlglot import parse_one

from giql.dialect import GIQLDialect
from giql.generators import BaseGIQLGenerator
from giql.generators import GIQLDuckDBGenerator
from giql.schema import ColumnInfo
from giql.schema import SchemaInfo
from giql.schema import TableSchema
from giql.transformer import ClusterTransformer
from giql.transformer import MergeTransformer

DialectType = Literal["duckdb", "sqlite"]


class GIQLEngine:
    """Multi-backend GIQL query engine.

    Supports multiple SQL databases through transpilation of GIQL syntax
    to standard SQL. Can work with DuckDB, SQLite, and other backends.

    Examples
    --------
    Query a pandas DataFrame with DuckDB::

        import pandas as pd

        df = pd.DataFrame(
            {
                "id": [1, 2, 3],
                "chromosome": ["chr1", "chr1", "chr2"],
                "start_pos": [1500, 10500, 500],
                "end_pos": [1600, 10600, 600],
            }
        )
        with GIQLEngine(target_dialect="duckdb") as engine:
            engine.conn.register("variants", df)
            result = engine.query(
                "SELECT * FROM variants WHERE position INTERSECTS 'chr1:1000-2000'"
            )

    Load from CSV::

        with GIQLEngine(target_dialect="duckdb") as engine:
            engine.load_csv("variants", "variants.csv")
            result = engine.query(
                "SELECT * FROM variants WHERE position INTERSECTS 'chr1:1000-2000'"
            )

    Using SQLite backend::

        with GIQLEngine(target_dialect="sqlite", db_path="data.db") as engine:
            result = engine.query(
                "SELECT * FROM variants WHERE position INTERSECTS 'chr1:1000-2000'"
            )
    """

    def __init__(
        self,
        target_dialect: DialectType | str = "duckdb",
        connection=None,
        db_path: str = ":memory:",
        verbose: bool = False,
        **dialect_options,
    ):
        """Initialize engine.

        :param target_dialect:
            Target SQL dialect ('duckdb', 'sqlite', 'standard')
        :param connection:
            Existing database connection (optional)
        :param db_path:
            Database path or connection string
        :param verbose:
            Print transpiled SQL
        :param dialect_options:
            Additional options for specific dialects
        """
        self.target_dialect = target_dialect
        self.verbose = verbose
        self.schema_info = SchemaInfo()
        self.dialect_options = dialect_options

        # Initialize connection
        if connection:
            self.conn = connection
            self.owns_connection = False
        else:
            self.conn = self._create_connection(db_path)
            self.owns_connection = True

        # Get appropriate generator
        self.generator = self._get_generator()

        # Initialize query transformers
        self.cluster_transformer = ClusterTransformer(self.schema_info)
        self.merge_transformer = MergeTransformer(self.schema_info)

    def _create_connection(self, db_path: str):
        """Create database connection based on target dialect.

        :param db_path:
            Path to database file or connection string
        :return:
            Connection object for the specified database backend
        :raises ImportError:
            If the required database driver is not installed
        :raises ValueError:
            If the dialect is not supported
        """
        if self.target_dialect == "duckdb":
            try:
                import duckdb

                return duckdb.connect(db_path)
            except ImportError:
                raise ImportError("DuckDB not installed.")

        elif self.target_dialect == "sqlite":
            import sqlite3

            return sqlite3.connect(db_path)

        else:
            raise ValueError(
                f"Unsupported dialect: {self.target_dialect}. Supported: duckdb, sqlite"
            )

    def _get_generator(self):
        """Get generator for target dialect.

        :return:
            SQL generator instance configured for the target dialect
        """
        generators = {
            "duckdb": GIQLDuckDBGenerator,
            "sqlite": BaseGIQLGenerator,
            "standard": BaseGIQLGenerator,
        }

        generator_class = generators.get(self.target_dialect, BaseGIQLGenerator)
        return generator_class(schema_info=self.schema_info, **self.dialect_options)

    def register_table_schema(
        self,
        table_name: str,
        columns: dict[str, str],
        genomic_column: str = "position",
        chrom_col: str = "chromosome",
        start_col: str = "start_pos",
        end_col: str = "end_pos",
        strand_col: str | None = None,
    ):
        """Register schema for a table.

        This method tells the engine how genomic ranges are stored in the table,
        mapping logical genomic column names to physical column names.

        :param table_name:
            Table name
        :param columns:
            Dict of column_name -> type
        :param genomic_column:
            Logical name for genomic position
        :param chrom_col:
            Physical chromosome column
        :param start_col:
            Physical start position column
        :param end_col:
            Physical end position column
        :param strand_col:
            Physical strand column (optional)
        """
        column_infos = {}

        for col_name, col_type in columns.items():
            if col_name == genomic_column:
                column_infos[col_name] = ColumnInfo(
                    name=col_name,
                    type=col_type,
                    is_genomic=True,
                    chrom_col=chrom_col,
                    start_col=start_col,
                    end_col=end_col,
                    strand_col=strand_col,
                )
            else:
                column_infos[col_name] = ColumnInfo(
                    name=col_name, type=col_type, is_genomic=False
                )

        table_schema = TableSchema(table_name, column_infos)
        self.schema_info.register_table(table_name, table_schema)

    def load_csv(self, table_name: str, file_path: str):
        """Load CSV file into database.

        :param table_name:
            Name to assign to the table
        :param file_path:
            Path to the CSV file
        """
        if self.target_dialect == "duckdb":
            self.conn.execute(
                f"CREATE TABLE {table_name} "
                f"AS SELECT * FROM read_csv_auto('{file_path}')"
            )
        elif self.target_dialect == "sqlite":
            # Use pandas for SQLite
            df = pd.read_csv(file_path)
            df.to_sql(table_name, self.conn, if_exists="replace", index=False)

        if self.verbose:
            print(f"Loaded {table_name} from {file_path}")

    def load_parquet(self, table_name: str, file_path: str):
        """Load Parquet file into database.

        :param table_name:
            Name to assign to the table
        :param file_path:
            Path to the Parquet file
        """
        if self.target_dialect == "duckdb":
            self.conn.execute(
                f"CREATE TABLE {table_name} AS SELECT * FROM read_parquet('{file_path}')"
            )
        else:
            df = pd.read_parquet(file_path)
            df.to_sql(table_name, self.conn, if_exists="replace", index=False)

        if self.verbose:
            print(f"Loaded {table_name} from {file_path}")

    def query(self, giql: str) -> pd.DataFrame:
        """Execute a GIQL query.

        Parses the GIQL syntax, transpiles to target SQL dialect,
        and executes the query.

        :param giql:
            Query string with GIQL genomic extensions
        :return:
            Query results as a pandas DataFrame
        :raises ValueError:
            If the query cannot be parsed, transpiled, or executed
        """
        # Parse with GIQL dialect
        try:
            ast = parse_one(giql, dialect=GIQLDialect)
        except Exception as e:
            raise ValueError(f"Parse error: {e}\nQuery: {giql}")

        # Transform query (MERGE first, then CLUSTER)
        try:
            # Apply MERGE transformation (which may internally use CLUSTER)
            ast = self.merge_transformer.transform(ast)
            # Apply CLUSTER transformation for any standalone CLUSTER expressions
            ast = self.cluster_transformer.transform(ast)
        except Exception as e:
            raise ValueError(f"Transformation error: {e}")

        # Transpile to target dialect
        try:
            target_sql = self.generator.generate(ast)
        except Exception as e:
            raise ValueError(f"Transpilation error: {e}")

        if self.verbose:
            print(f"\n{'=' * 60}")
            print(f"Target Dialect: {self.target_dialect}")
            print("\nOriginal GIQL:")
            print(giql)
            print("\nTranspiled SQL:")
            print(target_sql)
            print(f"{'=' * 60}\n")

        # Execute
        try:
            return pd.read_sql(target_sql, self.conn)
        except Exception as e:
            raise ValueError(f"Execution error: {e}\nSQL: {target_sql}")

    def execute_raw(self, sql: str) -> pd.DataFrame:
        """Execute raw SQL directly, bypassing GIQL parsing.

        :param sql:
            Raw SQL query string
        :return:
            Query results as a pandas DataFrame
        """
        return pd.read_sql(sql, self.conn)

    def close(self):
        """Close database connection.

        Only closes connections created by the engine. If an external
        connection was provided during initialization, it is not closed.
        """
        if self.owns_connection and self.conn:
            self.conn.close()

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()
