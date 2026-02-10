Execution engines
=================

GIQL transpiles genomic queries to SQL that can be executed on any database
backend. This guide covers backend-specific considerations and tips.

.. contents::
   :local:
   :depth: 1

Supported Backends
------------------

GIQL generates SQL that works across database systems:

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

Using with DuckDB
-----------------

DuckDB is recommended for most genomic analysis use cases. It provides excellent
performance for analytical queries and handles large genomic datasets efficiently.

.. code-block:: python

   import duckdb
   from giql import transpile

   sql = transpile(
       """
       SELECT * FROM features
       WHERE interval INTERSECTS 'chr1:1000-2000'
       """,
       tables=["features"],
   )

   conn = duckdb.connect()
   conn.execute("CREATE TABLE features AS SELECT * FROM read_csv('features.bed', delim='\t')")
   result = conn.execute(sql).fetchdf()

**Advantages:**

- Fast analytical query performance
- Efficient columnar storage
- Good support for large datasets
- Rich SQL feature set
- In-memory and persistent options

Using with SQLite
-----------------

SQLite is a lightweight, embedded database suitable for smaller datasets or
when portability is important.

.. code-block:: python

   import sqlite3
   from giql import transpile

   sql = transpile(
       """
       SELECT * FROM features
       WHERE interval INTERSECTS 'chr1:1000-2000'
       """,
       tables=["features"],
   )

   conn = sqlite3.connect("data.db")
   cursor = conn.execute(sql)
   for row in cursor:
       print(row)

**Advantages:**

- Zero configuration
- Single-file database
- Widely compatible
- Small memory footprint

Writing Portable Queries
------------------------

Query Compatibility
~~~~~~~~~~~~~~~~~~~

GIQL queries are portable across backends. The same GIQL query produces SQL
that works on any supported database:

.. code-block:: python

   from giql import transpile

   query = """
       SELECT a.*, b.name AS gene
       FROM variants a
       JOIN genes b ON a.interval INTERSECTS b.interval
       WHERE a.quality >= 30
   """

   # Same GIQL query works for any backend
   sql = transpile(query, tables=["variants", "genes"])

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
