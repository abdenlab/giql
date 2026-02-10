Quick Start
===========

GIQL provides a familiar SQL syntax for bioinformatics workflows, allowing
you to express complex genomic range operations without writing intricate
SQL expressions. GIQL queries read naturally, making your analysis code 
easier to review and share. GIQL operators follow established conventions 
around genomic spatial relationships, so the semantics are familiar and 
predictable.

- **Spatial operators**: INTERSECTS, CONTAINS, WITHIN for range relationships
- **Distance operators**: DISTANCE, NEAREST for proximity queries
- **Aggregation operators**: CLUSTER, MERGE for combining intervals
- **Set quantifiers**: ANY, ALL for multi-range queries
- **Range parsing**: Understands genomic range strings and coordinate systems
- **Transpilation**: Converts GIQL to standard SQL-92 compatible output for execution on any backend

Installation
------------

Install GIQL using pip:

.. code-block:: bash

   pip install giql

Basic Usage
-----------

Table Configuration
~~~~~~~~~~~~~~~~~~~

GIQL works with genomic data stored in tables with separate columns for chromosome,
start position, and end position. The default column names are:

* **chrom**: Chromosome identifier (e.g., 'chr1', 'chr2', 'chrX')
* **start**: Start position of the genomic interval (0-based, inclusive)
* **end**: End position of the genomic interval (0-based, exclusive, half-open)
* **strand** (optional): Strand orientation ('+', '-', or '.')

If your table uses the default column names, you can pass just the table name
as a string. For custom column names, use a ``Table`` object:

.. code-block:: python

   from giql import Table, transpile

   # Default column names (chrom, start, end, strand)
   sql = transpile(query, tables=["peaks"])

   # Custom column names
   sql = transpile(
       query,
       tables=[
           Table(
               "variants",
               genomic_col="interval",
               chrom_col="chromosome",
               start_col="start_pos",
               end_col="end_pos",
           )
       ],
   )

After configuration, you can use the genomic pseudo-column (default: ``interval``)
in your GIQL queries, and the transpiler will automatically expand it to the
physical column comparisons.

Query with DuckDB
~~~~~~~~~~~~~~~~~

.. code-block:: python

   import duckdb
   from giql import transpile

   sql = transpile(
       """
       SELECT * FROM variants
       WHERE interval INTERSECTS 'chr1:1000-2000'
       """,
       tables=["variants"],
   )

   conn = duckdb.connect()
   conn.execute("CREATE TABLE variants AS SELECT * FROM read_csv('variants.csv')")
   df = conn.execute(sql).fetchdf()

Query with SQLite
~~~~~~~~~~~~~~~~~

.. code-block:: python

   import sqlite3
   from giql import transpile

   sql = transpile(
       """
       SELECT * FROM variants
       WHERE interval INTERSECTS 'chr1:1000-2000'
       """,
       tables=["variants"],
   )

   conn = sqlite3.connect("data.db")
   cursor = conn.execute(sql)
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
