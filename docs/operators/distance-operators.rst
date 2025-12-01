Distance and Proximity Operators
================================

Distance and proximity operators calculate genomic distances and find nearest features.
These operators are essential for proximity analysis, such as finding genes near
regulatory elements or variants near transcription start sites.

.. contents::
   :local:
   :depth: 2

.. _distance-operator:

DISTANCE
--------

Calculate the genomic distance between two intervals.

Description
~~~~~~~~~~~

The ``DISTANCE`` operator returns the number of base pairs separating two genomic
intervals. It follows standard genomic distance conventions:

- **Overlapping intervals**: Returns ``0``
- **Non-overlapping intervals**: Returns the gap in base pairs (positive integer)
- **Different chromosomes**: Returns ``NULL``

Syntax
~~~~~~

.. code-block:: sql

   DISTANCE(interval_a, interval_b)

Parameters
~~~~~~~~~~

**interval_a**
   A genomic column registered with the engine.

**interval_b**
   Another genomic column to measure distance to.

Return Value
~~~~~~~~~~~~

- ``0`` for overlapping intervals
- Positive integer (gap in base pairs) for non-overlapping same-chromosome intervals
- ``NULL`` for intervals on different chromosomes

Examples
~~~~~~~~

**Calculate Distances Between Features:**

Calculate distance between peaks and genes:

.. code-block:: python

   cursor = engine.execute("""
       SELECT
           p.name AS peak,
           g.name AS gene,
           DISTANCE(p.interval, g.interval) AS distance
       FROM peaks p
       CROSS JOIN genes g
       WHERE p.chromosome = g.chromosome
       ORDER BY p.name, distance
   """)

**Filter by Distance:**

Find features within 10kb of each other:

.. code-block:: python

   cursor = engine.execute("""
       SELECT a.name, b.name, DISTANCE(a.interval, b.interval) AS dist
       FROM features_a a
       CROSS JOIN features_b b
       WHERE a.chromosome = b.chromosome
         AND DISTANCE(a.interval, b.interval) <= 10000
   """)

**Identify Overlapping vs. Proximal:**

Distinguish between overlapping and nearby features:

.. code-block:: python

   cursor = engine.execute("""
       SELECT
           p.name,
           g.name,
           CASE
               WHEN DISTANCE(p.interval, g.interval) = 0 THEN 'overlapping'
               WHEN DISTANCE(p.interval, g.interval) <= 1000 THEN 'proximal'
               ELSE 'distant'
           END AS relationship
       FROM peaks p
       CROSS JOIN genes g
       WHERE p.chromosome = g.chromosome
   """)

Backend Compatibility
~~~~~~~~~~~~~~~~~~~~~

.. list-table::
   :header-rows: 1
   :widths: 20 20 60

   * - Backend
     - Support
     - Notes
   * - DuckDB
     - Full
     -
   * - SQLite
     - Full
     -
   * - PostgreSQL
     - Planned
     -

Performance Notes
~~~~~~~~~~~~~~~~~

- Always include ``WHERE a.chromosome = b.chromosome`` to avoid unnecessary
  cross-chromosome comparisons
- For large datasets, consider pre-filtering by region before calculating distances
- Create indexes on chromosome and position columns for better performance

Related Operators
~~~~~~~~~~~~~~~~~

- :ref:`NEAREST <nearest-operator>` - Find k-nearest features (uses distance internally)
- :ref:`INTERSECTS <intersects-operator>` - Alternative for checking overlap (returns boolean)

----

.. _nearest-operator:

NEAREST
-------

Find the k-nearest genomic features to a reference point or interval.

Description
~~~~~~~~~~~

The ``NEAREST`` operator performs k-nearest neighbor (k-NN) queries on genomic data.
It finds the closest features from a target table relative to a reference position,
supporting various filtering options including strand awareness and distance constraints.

This operator uses ``CROSS JOIN LATERAL`` syntax to efficiently find nearest neighbors
for each row in the driving table.

Syntax
~~~~~~

.. code-block:: sql

   -- Find k nearest features for each row
   SELECT *
   FROM source_table
   CROSS JOIN LATERAL NEAREST(
       target_table,
       reference=source_table.interval,
       k=5
   ) AS nearest

   -- With additional parameters
   NEAREST(
       target_table,
       reference=interval,
       k=5,
       max_distance=100000,
       stranded=true,
       signed=true
   )

   -- Standalone query with literal reference
   SELECT * FROM NEAREST(genes, reference='chr1:1000000-1001000', k=5)

Parameters
~~~~~~~~~~

**target_table**
   The table to search for nearest features.

**reference**
   The reference position to measure distances from. Can be a column reference
   (e.g., ``peaks.interval``) or a literal range (e.g., ``'chr1:1000-2000'``).

**k**
   The number of nearest neighbors to return. Default: ``1``.

**max_distance** *(optional)*
   Maximum distance threshold. Only features within this distance are returned.

**stranded** *(optional)*
   When ``true``, only consider features on the same strand. Default: ``false``.

**signed** *(optional)*
   When ``true``, return signed distances (negative = upstream, positive = downstream).
   Default: ``false``.

Return Value
~~~~~~~~~~~~

Returns rows from the target table with an additional ``distance`` column indicating
the distance to the reference position. Results are ordered by distance (closest first).

Examples
~~~~~~~~

**Find K Nearest Genes:**

Find the 3 nearest genes for each peak:

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

**Standalone Query:**

Find 5 nearest genes to a specific genomic location:

.. code-block:: python

   cursor = engine.execute("""
       SELECT gene_name, distance
       FROM NEAREST(genes, reference='chr1:1000000-1001000', k=5)
       ORDER BY distance
   """)

**Distance-Constrained Search:**

Find nearest features within 100kb:

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
           k=5,
           max_distance=100000
       ) AS nearest
       ORDER BY peaks.name, nearest.distance
   """)

**Strand-Specific Nearest Neighbors:**

Find nearest same-strand features:

.. code-block:: python

   cursor = engine.execute("""
       SELECT
           peaks.name,
           nearest.name AS gene,
           nearest.strand,
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

**Directional (Upstream/Downstream) Queries:**

Find upstream features using signed distances:

.. code-block:: python

   # Upstream features have negative distances
   cursor = engine.execute("""
       SELECT
           peaks.name,
           nearest.name AS gene,
           nearest.distance
       FROM peaks
       CROSS JOIN LATERAL NEAREST(
           genes,
           reference=peaks.interval,
           k=10,
           signed=true
       ) AS nearest
       WHERE nearest.distance < 0
       ORDER BY peaks.name, nearest.distance DESC
   """)

   # Downstream features have positive distances
   cursor = engine.execute("""
       SELECT
           peaks.name,
           nearest.name AS gene,
           nearest.distance
       FROM peaks
       CROSS JOIN LATERAL NEAREST(
           genes,
           reference=peaks.interval,
           k=10,
           signed=true
       ) AS nearest
       WHERE nearest.distance > 0
       ORDER BY peaks.name, nearest.distance
   """)

**Combined Parameters:**

Find nearby same-strand features within distance constraints:

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
           k=5,
           max_distance=50000,
           stranded=true,
           signed=true
       ) AS nearest
       WHERE nearest.distance BETWEEN -10000 AND 10000
       ORDER BY peaks.name, ABS(nearest.distance)
   """)

Backend Compatibility
~~~~~~~~~~~~~~~~~~~~~

.. list-table::
   :header-rows: 1
   :widths: 20 20 60

   * - Backend
     - Support
     - Notes
   * - DuckDB
     - Full
     - Efficient lateral join support
   * - SQLite
     - Partial
     - Works but slower for large k values
   * - PostgreSQL
     - Planned
     -

Performance Notes
~~~~~~~~~~~~~~~~~

- **Chromosome pre-filtering**: NEAREST automatically filters by chromosome for efficiency
- **Use max_distance**: Specifying a maximum distance reduces the search space significantly
- **Limit k**: Only request as many neighbors as you actually need
- **Create indexes**: Add indexes on ``(chromosome, start_pos, end_pos)`` for better performance

.. code-block:: python

   # Create indexes for better NEAREST performance
   engine.conn.execute("""
       CREATE INDEX idx_genes_position
       ON genes (chromosome, start_pos, end_pos)
   """)

Related Operators
~~~~~~~~~~~~~~~~~

- :ref:`DISTANCE <distance-operator>` - Calculate distance between specific pairs
- :ref:`INTERSECTS <intersects-operator>` - Find overlapping features (distance = 0)
