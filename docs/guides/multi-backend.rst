Multi-Backend Guide
===================

GIQL supports multiple database backends, allowing you to run the same genomic
queries against different database systems. This guide covers backend selection,
configuration, and backend-specific considerations.

.. contents::
   :local:
   :depth: 2

Supported Backends
------------------

GIQL currently supports the following database backends:

.. list-table::
   :header-rows: 1
   :widths: 20 20 60

   * - Backend
     - Status
     - Best For
   * - DuckDB
     - Full Support
     - Analytics, large datasets, in-memory processing
   * - SQLite
     - Full Support
     - Lightweight, embedded, portable databases
   * - PostgreSQL
     - Planned
     - Production deployments, shared databases

Selecting a Backend
-------------------

DuckDB (Recommended)
~~~~~~~~~~~~~~~~~~~~

DuckDB is the recommended backend for most use cases. It provides excellent
performance for analytical queries and handles large genomic datasets efficiently.

.. code-block:: python

   from giql import GIQLEngine

   # In-memory DuckDB (default)
   with GIQLEngine(target_dialect="duckdb") as engine:
       engine.load_csv("features", "features.bed")
       # ... register schemas and query

   # Persistent DuckDB database
   with GIQLEngine(target_dialect="duckdb", db_path="my_data.duckdb") as engine:
       # Data persists between sessions
       pass

**Advantages:**

- Fast analytical query performance
- Efficient columnar storage
- Good support for large datasets
- Rich SQL feature set
- In-memory and persistent options

**Best for:**

- Interactive analysis
- Large BED/VCF files
- Complex aggregations
- One-time analysis pipelines

SQLite
~~~~~~

SQLite is a lightweight, embedded database suitable for smaller datasets or
when portability is important.

.. code-block:: python

   # In-memory SQLite
   with GIQLEngine(target_dialect="sqlite") as engine:
       pass

   # Persistent SQLite database
   with GIQLEngine(target_dialect="sqlite", db_path="my_data.db") as engine:
       pass

**Advantages:**

- Zero configuration
- Single-file database
- Widely compatible
- Small memory footprint

**Best for:**

- Small to medium datasets
- Portable analysis
- Embedded applications
- Simple workflows

Backend Configuration
---------------------

In-Memory vs Persistent
~~~~~~~~~~~~~~~~~~~~~~~

Both DuckDB and SQLite support in-memory and persistent modes:

.. code-block:: python

   # In-memory (data lost when engine closes)
   with GIQLEngine(target_dialect="duckdb") as engine:
       engine.load_csv("features", "features.bed")
       # Data exists only during this session

   # Persistent (data saved to disk)
   with GIQLEngine(target_dialect="duckdb", db_path="analysis.duckdb") as engine:
       engine.load_csv("features", "features.bed")
       # Data persists after engine closes

   # Reopen persistent database
   with GIQLEngine(target_dialect="duckdb", db_path="analysis.duckdb") as engine:
       # Previous data is available
       cursor = engine.execute("SELECT * FROM features LIMIT 5")

Connection Options
~~~~~~~~~~~~~~~~~~

Pass additional connection options to the underlying database:

.. code-block:: python

   # DuckDB with custom settings
   with GIQLEngine(
       target_dialect="duckdb",
       db_path="analysis.duckdb",
       read_only=False,
   ) as engine:
       pass

Writing Portable Queries
------------------------

Query Compatibility
~~~~~~~~~~~~~~~~~~~

GIQL queries are portable across backends. The same query works on any
supported database:

.. code-block:: python

   query = """
       SELECT a.*, b.name AS gene
       FROM variants a
       JOIN genes b ON a.interval INTERSECTS b.interval
       WHERE a.quality >= 30
   """

   # Works on DuckDB
   with GIQLEngine(target_dialect="duckdb") as engine:
       # ... setup ...
       cursor = engine.execute(query)

   # Same query works on SQLite
   with GIQLEngine(target_dialect="sqlite") as engine:
       # ... setup ...
       cursor = engine.execute(query)

SQL Dialect Differences
~~~~~~~~~~~~~~~~~~~~~~~

While GIQL queries are portable, the generated SQL differs between backends.
Use ``transpile()`` to see the backend-specific SQL:

.. code-block:: python

   query = "SELECT * FROM features WHERE interval INTERSECTS 'chr1:1000-2000'"

   # DuckDB SQL
   with GIQLEngine(target_dialect="duckdb") as engine:
       engine.register_table_schema("features", {...}, genomic_column="interval")
       print(engine.transpile(query))

   # SQLite SQL (may differ slightly)
   with GIQLEngine(target_dialect="sqlite") as engine:
       engine.register_table_schema("features", {...}, genomic_column="interval")
       print(engine.transpile(query))

Backend-Specific Features
~~~~~~~~~~~~~~~~~~~~~~~~~

Some SQL features may only be available on certain backends:

.. list-table::
   :header-rows: 1
   :widths: 40 20 20 20

   * - Feature
     - DuckDB
     - SQLite
     - Notes
   * - Window functions
     - Yes
     - Yes
     - Full support
   * - CTEs (WITH clause)
     - Yes
     - Yes
     - Full support
   * - LATERAL joins
     - Yes
     - Limited
     - Used by NEAREST
   * - STRING_AGG
     - Yes
     - GROUP_CONCAT
     - Different function names

Migrating Between Backends
--------------------------

Exporting Data
~~~~~~~~~~~~~~

Export data from one backend for import into another:

.. code-block:: python

   # Export from DuckDB
   with GIQLEngine(target_dialect="duckdb", db_path="source.duckdb") as engine:
       cursor = engine.execute("SELECT * FROM features")
       import pandas as pd
       df = pd.DataFrame(cursor.fetchall(),
                        columns=[desc[0] for desc in cursor.description])
       df.to_csv("features_export.csv", index=False)

   # Import to SQLite
   with GIQLEngine(target_dialect="sqlite", db_path="target.db") as engine:
       engine.load_csv("features", "features_export.csv")
       engine.register_table_schema("features", {...}, genomic_column="interval")

Schema Compatibility
~~~~~~~~~~~~~~~~~~~~

Ensure schema definitions work across backends:

.. code-block:: python

   # Use portable type names
   schema = {
       "chromosome": "VARCHAR",   # Works on all backends
       "start_pos": "BIGINT",     # Maps to appropriate integer type
       "end_pos": "BIGINT",
       "name": "VARCHAR",
       "score": "FLOAT",          # Maps to appropriate float type
   }

   # Same schema works on both backends
   for dialect in ["duckdb", "sqlite"]:
       with GIQLEngine(target_dialect=dialect) as engine:
           engine.load_csv("features", "features.csv")
           engine.register_table_schema("features", schema, genomic_column="interval")

Performance Comparison
----------------------

Backend Performance Characteristics
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. list-table::
   :header-rows: 1
   :widths: 30 35 35

   * - Operation
     - DuckDB
     - SQLite
   * - Large table scans
     - Excellent (columnar)
     - Good
   * - Complex joins
     - Excellent
     - Good
   * - Aggregations
     - Excellent
     - Good
   * - Small queries
     - Good
     - Excellent
   * - Memory usage
     - Higher
     - Lower
   * - Startup time
     - Faster
     - Fast

Choosing the Right Backend
~~~~~~~~~~~~~~~~~~~~~~~~~~

**Choose DuckDB when:**

- Working with large datasets (millions of features)
- Running complex analytical queries
- Performing heavy aggregations
- Memory is not constrained

**Choose SQLite when:**

- Working with smaller datasets
- Need maximum portability
- Memory is constrained
- Simple query patterns

Using External Connections
--------------------------

Connecting to Existing Databases
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Connect to databases created outside of GIQL:

.. code-block:: python

   # Connect to existing DuckDB database
   with GIQLEngine(target_dialect="duckdb", db_path="existing.duckdb") as engine:
       # Register schemas for existing tables
       engine.register_table_schema(
           "my_existing_table",
           {
               "chromosome": "VARCHAR",
               "start_pos": "BIGINT",
               "end_pos": "BIGINT",
               "name": "VARCHAR",
           },
           genomic_column="interval",
       )

       # Query existing data with GIQL operators
       cursor = engine.execute("""
           SELECT * FROM my_existing_table
           WHERE interval INTERSECTS 'chr1:1000-2000'
       """)

Using Transpiled SQL Externally
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Generate SQL for use with external database connections:

.. code-block:: python

   import duckdb

   # Get transpiled SQL from GIQL
   with GIQLEngine(target_dialect="duckdb") as engine:
       engine.register_table_schema("features", {...}, genomic_column="interval")
       sql = engine.transpile("""
           SELECT * FROM features
           WHERE interval INTERSECTS 'chr1:1000-2000'
       """)

   # Execute with external connection
   conn = duckdb.connect("my_database.duckdb")
   result = conn.execute(sql).fetchall()
   conn.close()

This is useful when integrating GIQL with existing database workflows or
when you need more control over the database connection.
