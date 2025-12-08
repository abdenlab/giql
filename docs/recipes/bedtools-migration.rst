Bedtools Migration Guide
========================

This guide maps bedtools commands to their GIQL equivalents. If you're familiar
with bedtools and want to replicate specific operations in GIQL, use this
reference to find the corresponding query patterns.

.. contents::
   :local:
   :depth: 2

Quick Reference Table
---------------------

.. list-table::
   :header-rows: 1
   :widths: 35 45 20

   * - Bedtools Command
     - GIQL Equivalent
     - Recipe
   * - ``intersect -a A -b B``
     - ``SELECT DISTINCT a.* FROM a, b WHERE a.pos INTERSECTS b.pos``
     - :ref:`intersect-basic`
   * - ``intersect -a A -b B -wa``
     - ``SELECT a.* FROM a, b WHERE a.pos INTERSECTS b.pos``
     - :ref:`intersect-wa`
   * - ``intersect -a A -b B -wb``
     - ``SELECT b.* FROM a, b WHERE a.pos INTERSECTS b.pos``
     - :ref:`intersect-wb`
   * - ``intersect -a A -b B -wa -wb``
     - ``SELECT a.*, b.* FROM a, b WHERE a.pos INTERSECTS b.pos``
     - :ref:`intersect-wawb`
   * - ``intersect -a A -b B -v``
     - ``SELECT a.* FROM a LEFT JOIN b ... WHERE b.chr IS NULL``
     - :ref:`intersect-v`
   * - ``intersect -a A -b B -u``
     - ``SELECT DISTINCT a.* FROM a JOIN b ...``
     - :ref:`intersect-u`
   * - ``intersect -a A -b B -c``
     - ``SELECT a.*, COUNT(b.name) ... GROUP BY ...``
     - :ref:`intersect-c`
   * - ``intersect -a A -b B -wo``
     - ``SELECT a.*, b.*, (overlap calculation) ...``
     - :ref:`intersect-wo`
   * - ``intersect -a A -b B -loj``
     - ``SELECT a.*, b.* FROM a LEFT JOIN b ...``
     - :ref:`intersect-loj`
   * - ``closest -a A -b B -k N``
     - ``CROSS JOIN LATERAL NEAREST(b, reference=a.pos, k=N)``
     - :ref:`closest-k`
   * - ``closest -a A -b B -d``
     - ``SELECT ..., DISTANCE(a.pos, b.pos) ...``
     - :ref:`closest-d`
   * - ``cluster -i A``
     - ``SELECT *, CLUSTER(interval) AS cluster_id FROM a``
     - :ref:`cluster-basic`
   * - ``cluster -i A -d N``
     - ``SELECT *, CLUSTER(interval, N) AS cluster_id FROM a``
     - :ref:`cluster-d`
   * - ``merge -i A``
     - ``SELECT MERGE(interval) FROM a``
     - :ref:`merge-basic`
   * - ``merge -i A -d N``
     - ``SELECT MERGE(interval, N) FROM a``
     - :ref:`merge-d`
   * - ``merge -i A -c 1 -o count``
     - ``SELECT MERGE(interval), COUNT(*) FROM a``
     - :ref:`merge-count`

bedtools intersect
------------------

.. _intersect-basic:

Default: Report overlaps between A and B
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Bedtools:**

.. code-block:: bash

   bedtools intersect -a file_a.bed -b file_b.bed

**GIQL:**

.. code-block:: python

   cursor = engine.execute("""
       SELECT DISTINCT a.*
       FROM features_a a, features_b b
       WHERE a.interval INTERSECTS b.interval
   """)

.. _intersect-wa:

``-wa``: Write original A entry for each overlap
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Bedtools:**

.. code-block:: bash

   bedtools intersect -a file_a.bed -b file_b.bed -wa

**GIQL:**

.. code-block:: python

   cursor = engine.execute("""
       SELECT a.*
       FROM features_a a, features_b b
       WHERE a.interval INTERSECTS b.interval
   """)

.. _intersect-wb:

``-wb``: Write original B entry for each overlap
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Bedtools:**

.. code-block:: bash

   bedtools intersect -a file_a.bed -b file_b.bed -wb

**GIQL:**

.. code-block:: python

   cursor = engine.execute("""
       SELECT b.*
       FROM features_a a, features_b b
       WHERE a.interval INTERSECTS b.interval
   """)

.. _intersect-wawb:

``-wa -wb``: Write both A and B entries
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Bedtools:**

.. code-block:: bash

   bedtools intersect -a file_a.bed -b file_b.bed -wa -wb

**GIQL:**

.. code-block:: python

   cursor = engine.execute("""
       SELECT a.*, b.*
       FROM features_a a, features_b b
       WHERE a.interval INTERSECTS b.interval
   """)

.. _intersect-v:

``-v``: Report A entries with NO overlap in B
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Bedtools:**

.. code-block:: bash

   bedtools intersect -a file_a.bed -b file_b.bed -v

**GIQL:**

.. code-block:: python

   cursor = engine.execute("""
       SELECT a.*
       FROM features_a a
       LEFT JOIN features_b b ON a.interval INTERSECTS b.interval
       WHERE b.chromosome IS NULL
   """)

.. _intersect-u:

``-u``: Report A entries with ANY overlap (unique)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Bedtools:**

.. code-block:: bash

   bedtools intersect -a file_a.bed -b file_b.bed -u

**GIQL:**

.. code-block:: python

   cursor = engine.execute("""
       SELECT DISTINCT a.*
       FROM features_a a
       INNER JOIN features_b b ON a.interval INTERSECTS b.interval
   """)

.. _intersect-c:

``-c``: Count B overlaps for each A feature
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Bedtools:**

.. code-block:: bash

   bedtools intersect -a file_a.bed -b file_b.bed -c

**GIQL:**

.. code-block:: python

   cursor = engine.execute("""
       SELECT a.*, COUNT(b.name) AS overlap_count
       FROM features_a a
       LEFT JOIN features_b b ON a.interval INTERSECTS b.interval
       GROUP BY a.chromosome, a.start_pos, a.end_pos, a.name, a.score, a.strand
   """)

.. _intersect-wo:

``-wo``: Write overlap amount in base pairs
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Bedtools:**

.. code-block:: bash

   bedtools intersect -a file_a.bed -b file_b.bed -wo

**GIQL:**

.. code-block:: python

   cursor = engine.execute("""
       SELECT
           a.*,
           b.*,
           (LEAST(a.end_pos, b.end_pos) - GREATEST(a.start_pos, b.start_pos)) AS overlap_bp
       FROM features_a a, features_b b
       WHERE a.interval INTERSECTS b.interval
   """)

.. _intersect-wao:

``-wao``: Write overlap amount for ALL A features
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Bedtools:**

.. code-block:: bash

   bedtools intersect -a file_a.bed -b file_b.bed -wao

**GIQL:**

.. code-block:: python

   cursor = engine.execute("""
       SELECT
           a.*,
           b.*,
           CASE
               WHEN b.chromosome IS NULL THEN 0
               ELSE LEAST(a.end_pos, b.end_pos) - GREATEST(a.start_pos, b.start_pos)
           END AS overlap_bp
       FROM features_a a
       LEFT JOIN features_b b ON a.interval INTERSECTS b.interval
   """)

.. _intersect-loj:

``-loj``: Left outer join
~~~~~~~~~~~~~~~~~~~~~~~~~

**Bedtools:**

.. code-block:: bash

   bedtools intersect -a file_a.bed -b file_b.bed -loj

**GIQL:**

.. code-block:: python

   cursor = engine.execute("""
       SELECT a.*, b.*
       FROM features_a a
       LEFT JOIN features_b b ON a.interval INTERSECTS b.interval
   """)

``-s``: Same strand overlaps only
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Bedtools:**

.. code-block:: bash

   bedtools intersect -a file_a.bed -b file_b.bed -s

**GIQL:**

.. code-block:: python

   cursor = engine.execute("""
       SELECT a.*
       FROM features_a a, features_b b
       WHERE a.interval INTERSECTS b.interval
         AND a.strand = b.strand
   """)

``-S``: Opposite strand overlaps only
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Bedtools:**

.. code-block:: bash

   bedtools intersect -a file_a.bed -b file_b.bed -S

**GIQL:**

.. code-block:: python

   cursor = engine.execute("""
       SELECT a.*
       FROM features_a a, features_b b
       WHERE a.interval INTERSECTS b.interval
         AND a.strand != b.strand
         AND a.strand IN ('+', '-')
         AND b.strand IN ('+', '-')
   """)

``-f``: Minimum overlap fraction of A
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Bedtools:**

.. code-block:: bash

   bedtools intersect -a file_a.bed -b file_b.bed -f 0.5

**GIQL:**

.. code-block:: python

   cursor = engine.execute("""
       SELECT a.*
       FROM features_a a, features_b b
       WHERE a.interval INTERSECTS b.interval
         AND (
             LEAST(a.end_pos, b.end_pos) - GREATEST(a.start_pos, b.start_pos)
         ) >= 0.5 * (a.end_pos - a.start_pos)
   """)

``-r``: Reciprocal overlap
~~~~~~~~~~~~~~~~~~~~~~~~~~

**Bedtools:**

.. code-block:: bash

   bedtools intersect -a file_a.bed -b file_b.bed -f 0.5 -r

**GIQL:**

.. code-block:: python

   cursor = engine.execute("""
       WITH overlap_calcs AS (
           SELECT
               a.*,
               (LEAST(a.end_pos, b.end_pos) - GREATEST(a.start_pos, b.start_pos)) AS overlap_bp,
               (a.end_pos - a.start_pos) AS a_length,
               (b.end_pos - b.start_pos) AS b_length
           FROM features_a a, features_b b
           WHERE a.interval INTERSECTS b.interval
       )
       SELECT chromosome, start_pos, end_pos, name, score, strand
       FROM overlap_calcs
       WHERE overlap_bp >= 0.5 * a_length
         AND overlap_bp >= 0.5 * b_length
   """)

bedtools closest
----------------

.. _closest-k:

``-k``: Find k nearest features
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Bedtools:**

.. code-block:: bash

   bedtools closest -a peaks.bed -b genes.bed -k 3

**GIQL:**

.. code-block:: python

   cursor = engine.execute("""
       SELECT
           peaks.name AS peak,
           nearest.name AS gene,
           nearest.distance
       FROM peaks
       CROSS JOIN LATERAL NEAREST(genes, reference=peaks.interval, k=3) AS nearest
       ORDER BY peaks.name, nearest.distance
   """)

.. _closest-d:

``-d``: Report distance
~~~~~~~~~~~~~~~~~~~~~~~

**Bedtools:**

.. code-block:: bash

   bedtools closest -a peaks.bed -b genes.bed -d

**GIQL:**

.. code-block:: python

   cursor = engine.execute("""
       SELECT
           a.name AS peak,
           b.name AS gene,
           DISTANCE(a.interval, b.interval) AS distance
       FROM peaks a
       CROSS JOIN genes b
       WHERE a.chromosome = b.chromosome
       ORDER BY a.name, distance
   """)

Or using NEAREST for just the closest:

.. code-block:: python

   cursor = engine.execute("""
       SELECT
           peaks.name AS peak,
           nearest.name AS gene,
           nearest.distance
       FROM peaks
       CROSS JOIN LATERAL NEAREST(genes, reference=peaks.interval, k=1) AS nearest
   """)

``-s``: Same strand only
~~~~~~~~~~~~~~~~~~~~~~~~

**Bedtools:**

.. code-block:: bash

   bedtools closest -a peaks.bed -b genes.bed -s -k 3

**GIQL:**

.. code-block:: python

   cursor = engine.execute("""
       SELECT
           peaks.name,
           nearest.name AS gene,
           nearest.distance
       FROM peaks
       CROSS JOIN LATERAL NEAREST(
           genes,
           reference=peaks.interval,
           k=3,
           stranded=true
       ) AS nearest
       ORDER BY peaks.name, nearest.distance
   """)

bedtools cluster
----------------

.. _cluster-basic:

Basic clustering
~~~~~~~~~~~~~~~~

**Bedtools:**

.. code-block:: bash

   bedtools cluster -i features.bed

**GIQL:**

.. code-block:: python

   cursor = engine.execute("""
       SELECT
           *,
           CLUSTER(interval) AS cluster_id
       FROM features
       ORDER BY chromosome, start_pos
   """)

.. _cluster-d:

``-d``: Cluster with distance parameter
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Bedtools:**

.. code-block:: bash

   bedtools cluster -i features.bed -d 1000

**GIQL:**

.. code-block:: python

   cursor = engine.execute("""
       SELECT
           *,
           CLUSTER(interval, 1000) AS cluster_id
       FROM features
       ORDER BY chromosome, start_pos
   """)

``-s``: Strand-specific clustering
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Bedtools:**

.. code-block:: bash

   bedtools cluster -i features.bed -s

**GIQL:**

.. code-block:: python

   cursor = engine.execute("""
       SELECT
           *,
           CLUSTER(interval, stranded=true) AS cluster_id
       FROM features
       ORDER BY chromosome, strand, start_pos
   """)

bedtools merge
--------------

.. _merge-basic:

Basic merge
~~~~~~~~~~~

**Bedtools:**

.. code-block:: bash

   bedtools merge -i features.bed

**GIQL:**

.. code-block:: python

   cursor = engine.execute("""
       SELECT MERGE(interval)
       FROM features
   """)

.. _merge-d:

``-d``: Merge with distance parameter
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Bedtools:**

.. code-block:: bash

   bedtools merge -i features.bed -d 1000

**GIQL:**

.. code-block:: python

   cursor = engine.execute("""
       SELECT MERGE(interval, 1000)
       FROM features
   """)

``-s``: Strand-specific merge
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Bedtools:**

.. code-block:: bash

   bedtools merge -i features.bed -s

**GIQL:**

.. code-block:: python

   cursor = engine.execute("""
       SELECT MERGE(interval, stranded=true)
       FROM features
   """)

.. _merge-count:

``-c -o count``: Count merged features
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Bedtools:**

.. code-block:: bash

   bedtools merge -i features.bed -c 1 -o count

**GIQL:**

.. code-block:: python

   cursor = engine.execute("""
       SELECT
           MERGE(interval),
           COUNT(*) AS feature_count
       FROM features
   """)

``-c -o mean``: Average score
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Bedtools:**

.. code-block:: bash

   bedtools merge -i features.bed -c 5 -o mean

**GIQL:**

.. code-block:: python

   cursor = engine.execute("""
       SELECT
           MERGE(interval),
           AVG(score) AS avg_score
       FROM features
   """)

``-c -o collapse``: Collect names
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Bedtools:**

.. code-block:: bash

   bedtools merge -i features.bed -c 4 -o collapse

**GIQL:**

.. code-block:: python

   cursor = engine.execute("""
       SELECT
           MERGE(interval),
           STRING_AGG(name, ',') AS feature_names
       FROM features
   """)

Key Differences from Bedtools
-----------------------------

1. **SQL-based syntax**: GIQL uses SQL syntax, which may be more familiar to
   users with database experience and allows integration with other SQL features.

2. **Explicit joins**: Instead of implicit A/B file relationships, GIQL uses
   explicit JOIN syntax, making the relationship between tables clearer.

3. **Flexible output**: SQL's SELECT clause gives you full control over which
   columns to return and how to format them.

4. **Built-in aggregation**: SQL's GROUP BY and aggregate functions (COUNT, AVG,
   SUM, etc.) are available directly, without needing separate post-processing.

5. **Database integration**: GIQL queries run against database tables, enabling
   integration with other data and persistence of results.

6. **Multi-backend support**: The same GIQL query can run on DuckDB, SQLite,
   or other supported backends without modification.
