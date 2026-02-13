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

Test if ranges overlap.

.. code-block:: sql

   -- Against literal
   interval INTERSECTS 'chr1:1000-2000'

   -- Column to column
   a.interval INTERSECTS b.interval

   -- In JOIN
   JOIN table ON a.interval INTERSECTS b.interval

CONTAINS
~~~~~~~~

Test if one range fully contains another.

.. code-block:: sql

   -- Range contains point
   interval CONTAINS 'chr1:1500'

   -- Range contains range
   interval CONTAINS 'chr1:1200-1800'

   -- Column to column
   gene.interval CONTAINS exon.interval

WITHIN
~~~~~~

Test if one range is fully within another.

.. code-block:: sql

   -- Range within literal
   interval WITHIN 'chr1:1000-5000'

   -- Column to column
   exon.interval WITHIN gene.interval

Distance Operators
------------------

DISTANCE
~~~~~~~~

Calculate distance between two positions.

.. code-block:: sql

   DISTANCE(a.interval, b.interval)

Returns:

- ``0`` for overlapping ranges
- Positive integer (gap in bp) for non-overlapping
- ``NULL`` for different chromosomes

NEAREST
~~~~~~~

Find k-nearest neighbors.

.. code-block:: sql

   -- Basic syntax
   CROSS JOIN LATERAL NEAREST(
       target_table,
       reference=source.interval,
       k=N
   ) AS alias

   -- With parameters
   NEAREST(
       target_table,
       reference=interval,
       k=5,
       max_distance=100000,
       stranded=true,
       signed=true
   )

   -- Standalone
   SELECT * FROM NEAREST(table, reference='chr1:1000-2000', k=5)

Parameters:

- ``k``: Number of neighbors (default: 1)
- ``max_distance``: Maximum distance threshold
- ``stranded``: Same-strand only (default: false)
- ``signed``: Signed distances (default: false)

Aggregation Operators
---------------------

CLUSTER
~~~~~~~

Assign cluster IDs to overlapping intervals.

.. code-block:: sql

   -- Basic
   CLUSTER(interval) AS cluster_id

   -- With distance
   CLUSTER(interval, 1000) AS cluster_id

   -- Strand-specific
   CLUSTER(interval, stranded=true) AS cluster_id

   -- Combined
   CLUSTER(interval, 1000, stranded=true) AS cluster_id

MERGE
~~~~~

Combine overlapping intervals.

.. code-block:: sql

   -- Basic
   SELECT MERGE(interval) FROM table

   -- With distance
   SELECT MERGE(interval, 1000) FROM table

   -- Strand-specific
   SELECT MERGE(interval, stranded=true) FROM table

   -- With aggregations
   SELECT MERGE(interval), COUNT(*), AVG(score) FROM table

Set Quantifiers
---------------

ANY
~~~

Match any of multiple ranges.

.. code-block:: sql

   interval INTERSECTS ANY('chr1:1000-2000', 'chr2:5000-6000')
   interval CONTAINS ANY('chr1:1500', 'chr1:2500')
   interval WITHIN ANY('chr1:0-10000', 'chr2:0-10000')

ALL
~~~

Match all of multiple ranges.

.. code-block:: sql

   interval CONTAINS ALL('chr1:1500', 'chr1:1600', 'chr1:1700')
   interval INTERSECTS ALL('chr1:1000-1100', 'chr1:1050-1150')

Query Patterns
--------------

Basic Filter
~~~~~~~~~~~~

.. code-block:: sql

   SELECT * FROM table
   WHERE interval INTERSECTS 'chr1:1000-2000'

Join
~~~~

.. code-block:: sql

   SELECT a.*, b.name
   FROM table_a a
   JOIN table_b b ON a.interval INTERSECTS b.interval

Left Outer Join
~~~~~~~~~~~~~~~

.. code-block:: sql

   SELECT a.*, b.name
   FROM table_a a
   LEFT JOIN table_b b ON a.interval INTERSECTS b.interval

Exclusion (NOT IN)
~~~~~~~~~~~~~~~~~~

.. code-block:: sql

   SELECT a.*
   FROM table_a a
   LEFT JOIN table_b b ON a.interval INTERSECTS b.interval
   WHERE b.chrom IS NULL

Count Overlaps
~~~~~~~~~~~~~~

.. code-block:: sql

   SELECT a.*, COUNT(b.name) AS overlap_count
   FROM table_a a
   LEFT JOIN table_b b ON a.interval INTERSECTS b.interval
   GROUP BY a.chrom, a.start, a."end", ...

K-Nearest Neighbors
~~~~~~~~~~~~~~~~~~~

.. code-block:: sql

   SELECT source.*, nearest.name, nearest.distance
   FROM source
   CROSS JOIN LATERAL NEAREST(target, reference=source.interval, k=5) AS nearest

Clustering
~~~~~~~~~~

.. code-block:: sql

   SELECT *, CLUSTER(interval) AS cluster_id
   FROM table
   ORDER BY chrom, start

Merging
~~~~~~~

.. code-block:: sql

   SELECT MERGE(interval), COUNT(*) AS count
   FROM table
