Quick Start
===========

Installation
------------

Install GIQL using pip:

.. code-block:: bash

   pip install giql

Or with optional dependencies:

.. code-block:: bash

   pip install giql[duckdb]  # For DuckDB support

Basic Usage
-----------

Expected Schema
~~~~~~~~~~~~~~~

GIQL works with genomic data stored in tables with separate columns for chromosome,
start position, and end position. The typical schema includes:

* **chromosome**: Chromosome identifier (e.g., 'chr1', 'chr2', 'chrX')
* **start_pos**: Start position of the genomic interval (0-based, inclusive)
* **end_pos**: End position of the genomic interval (0-based, exclusive, half-open)
* **strand** (optional): Strand orientation ('+', '-', or '.')

You must register the table schema with GIQL, mapping the logical genomic column
(used in queries) to the physical columns in your table:

.. code-block:: python

   engine.register_table_schema(
       "table_name",
       {
           "chromosome": "VARCHAR",
           "start_pos": "BIGINT",
           "end_pos": "BIGINT",
           "strand": "VARCHAR",        # Optional
           # ... other columns ...
       },
       genomic_column="interval",      # Logical name used in queries
   )

After registration, you can use ``interval`` in your GIQL queries, and the engine
will automatically map it to the ``chromosome``, ``start_pos``, and ``end_pos``
columns.

Query with DuckDB
~~~~~~~~~~~~~~~~~

.. code-block:: python

   from giql import GIQLEngine

   with GIQLEngine(target_dialect="duckdb") as engine:
       # Load CSV file into database
       engine.load_csv("variants", "variants.csv")

       # Register schema mapping
       engine.register_table_schema(
           "variants",
           {
               "id": "INTEGER",
               "chromosome": "VARCHAR",
               "start_pos": "BIGINT",
               "end_pos": "BIGINT",
           },
           genomic_column="interval",
       )

       # Query using the logical 'interval' column (returns cursor for streaming)
       cursor = engine.execute("""
           SELECT * FROM variants
           WHERE interval INTERSECTS 'chr1:1000-2000'
       """)

       # Process results lazily
       for row in cursor:
           print(row)

       # Or materialize to pandas DataFrame
       import pandas as pd
       cursor = engine.execute("SELECT ...")
       df = pd.DataFrame(cursor.fetchall(), columns=[desc[0] for desc in cursor.description])

Query with SQLite
~~~~~~~~~~~~~~~~~

.. code-block:: python

   from giql import GIQLEngine

   with GIQLEngine(target_dialect="sqlite", db_path="data.db") as engine:
       cursor = engine.execute("""
           SELECT * FROM variants
           WHERE interval INTERSECTS 'chr1:1000-2000'
       """)

       # Iterate results
       for row in cursor:
           print(row)

Spatial Operators
-----------------

INTERSECTS
~~~~~~~~~~

Check if genomic ranges overlap:

.. code-block:: sql

   SELECT * FROM variants
   WHERE interval INTERSECTS 'chr1:1000-2000'

CONTAINS
~~~~~~~~

Check if a range contains a point or another range:

.. code-block:: sql

   SELECT * FROM variants
   WHERE interval CONTAINS 'chr1:1500'

WITHIN
~~~~~~

Check if a range is within another range:

.. code-block:: sql

   SELECT * FROM variants
   WHERE interval WITHIN 'chr1:1000-5000'

Set Quantifiers
---------------

ANY
~~~

Match any of the specified ranges:

.. code-block:: sql

   SELECT * FROM variants
   WHERE interval INTERSECTS ANY('chr1:1000-2000', 'chr1:5000-6000')

ALL
~~~

Match all of the specified ranges:

.. code-block:: sql

   SELECT * FROM variants
   WHERE interval CONTAINS ALL('chr1:1500', 'chr1:1600')

Column-to-Column Joins
----------------------

Join tables on genomic position:

.. code-block:: sql

   SELECT v.*, g.name
   FROM variants v
   INNER JOIN genes g ON v.interval INTERSECTS g.interval

Transpiling to SQL
------------------

The ``transpile()`` method converts GIQL queries to standard SQL without executing them.
This is useful for debugging, understanding the generated SQL, or integrating with external tools:

.. code-block:: python

   from giql import GIQLEngine

   with GIQLEngine(target_dialect="duckdb") as engine:
       # Register table schema
       engine.register_table_schema(
           "variants",
           {
               "chromosome": "VARCHAR",
               "start_pos": "BIGINT",
               "end_pos": "BIGINT",
           },
           genomic_column="interval",
       )

       # Transpile GIQL to SQL
       sql = engine.transpile("""
           SELECT * FROM variants
           WHERE interval INTERSECTS 'chr1:1000-2000'
       """)

       print(sql)
       # Output: SELECT * FROM variants WHERE chromosome = 'chr1' AND start_pos < 2000 AND end_pos > 1000

Different target dialects generate different SQL:

.. code-block:: python

   # DuckDB dialect
   with GIQLEngine(target_dialect="duckdb") as engine:
       sql = engine.transpile("SELECT * FROM variants WHERE interval INTERSECTS 'chr1:1000-2000'")
       # Generates DuckDB-optimized SQL

   # SQLite dialect
   with GIQLEngine(target_dialect="sqlite") as engine:
       sql = engine.transpile("SELECT * FROM variants WHERE interval INTERSECTS 'chr1:1000-2000'")
       # Generates SQLite-compatible SQL

The transpiled SQL can be executed directly on your database or used with other tools.
Use ``verbose=True`` when creating the engine to see detailed transpilation information:

.. code-block:: python

   with GIQLEngine(target_dialect="duckdb", verbose=True) as engine:
       sql = engine.transpile("SELECT * FROM variants WHERE interval INTERSECTS 'chr1:1000-2000'")
       # Prints detailed information about the transpilation process
