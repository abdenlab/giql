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
       genomic_column="position",      # Logical name used in queries
   )

After registration, you can use ``position`` in your GIQL queries, and the engine
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
           genomic_column="position",
       )

       # Query using the logical 'position' column (returns cursor for streaming)
       cursor = engine.execute("""
           SELECT * FROM variants
           WHERE position INTERSECTS 'chr1:1000-2000'
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
           WHERE position INTERSECTS 'chr1:1000-2000'
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
   WHERE position INTERSECTS 'chr1:1000-2000'

CONTAINS
~~~~~~~~

Check if a range contains a point or another range:

.. code-block:: sql

   SELECT * FROM variants
   WHERE position CONTAINS 'chr1:1500'

WITHIN
~~~~~~

Check if a range is within another range:

.. code-block:: sql

   SELECT * FROM variants
   WHERE position WITHIN 'chr1:1000-5000'

Set Quantifiers
---------------

ANY
~~~

Match any of the specified ranges:

.. code-block:: sql

   SELECT * FROM variants
   WHERE position INTERSECTS ANY('chr1:1000-2000', 'chr1:5000-6000')

ALL
~~~

Match all of the specified ranges:

.. code-block:: sql

   SELECT * FROM variants
   WHERE position CONTAINS ALL('chr1:1500', 'chr1:1600')

Column-to-Column Joins
----------------------

Join tables on genomic position:

.. code-block:: sql

   SELECT v.*, g.name
   FROM variants v
   INNER JOIN genes g ON v.position INTERSECTS g.position
