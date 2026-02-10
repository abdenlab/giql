Syntax Reference
================

Quick reference for GIQL syntax and operators.

.. contents::
   :local:
   :depth: 1

Genomic Range Literals
----------------------

Format
~~~~~~

Genomic ranges are specified as string literals:

.. code-block:: text

   'chromosome:start-end'

Examples
~~~~~~~~

.. code-block:: sql

   'chr1:1000-2000'        -- Range on chr1 from 1000 to 2000
   'chr1:1000'             -- Point at position 1000
   'chrX:50000-100000'     -- Range on chrX
   'chr1:0-1000000'        -- First megabase of chr1

Coordinate System
~~~~~~~~~~~~~~~~~

- **0-based start**: First base is position 0
- **Half-open interval**: [start, end) - start inclusive, end exclusive
- Range ``chr1:100-200`` covers bases 100 through 199

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
