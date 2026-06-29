Distance and Neighbors
======================

Distance and proximity operators calculate genomic distances and find nearest features.
These operators are essential for proximity analysis, such as finding genes near
regulatory elements or variants near transcription start sites.

.. _distance-operator:

DISTANCE
--------

Calculate the genomic distance between two intervals.

Description
~~~~~~~~~~~

The ``DISTANCE`` operator returns the number of base pairs separating two genomic
intervals, matching ``bedtools closest -d`` semantics:

- **Overlapping intervals**: Returns ``0``
- **Book-ended (adjacent) intervals** (``A.end == B.start`` in half-open coordinates): Returns ``1``
- **Non-overlapping intervals**: Returns the half-open gap plus one (a raw gap of ``N`` bases reports ``N + 1``)
- **Different chromosomes**: Returns ``NULL``

Syntax
~~~~~~

.. code-block:: sql

   DISTANCE(interval_a, interval_b)

Parameters
~~~~~~~~~~

**interval_a**
   A genomic column.

**interval_b**
   Another genomic column to measure distance to.

Return Value
~~~~~~~~~~~~

- ``0`` for overlapping intervals
- ``1`` for book-ended (adjacent) intervals, matching ``bedtools closest -d``
- Positive integer (half-open gap ``+ 1``) for non-overlapping same-chromosome intervals
- ``NULL`` for intervals on different chromosomes

Examples
~~~~~~~~

**Calculate Distances Between Features:**

Calculate distance between peaks and genes:

.. code-block:: sql

   SELECT
       p.name AS peak,
       g.name AS gene,
       DISTANCE(p.interval, g.interval) AS distance
   FROM peaks p
   CROSS JOIN genes g
   WHERE p.chrom = g.chrom
   ORDER BY p.name, distance

**Filter by Distance:**

Find features within 10kb of each other:

.. code-block:: sql

   SELECT a.name, b.name, DISTANCE(a.interval, b.interval) AS dist
   FROM features_a a
   CROSS JOIN features_b b
   WHERE a.chrom = b.chrom
     AND DISTANCE(a.interval, b.interval) <= 10000

**Identify Overlapping vs. Proximal:**

Distinguish between overlapping and nearby features:

.. code-block:: sql

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
   WHERE p.chrom = g.chrom

Notes
~~~~~

- Always include ``WHERE a.chrom = b.chrom`` to avoid unnecessary
  cross-chromosome comparisons
- For large datasets, consider pre-filtering by region before calculating distances

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
       reference := source_table.interval,
       k := 5
   ) AS nearest

   -- With additional parameters
   NEAREST(
       target_table,
       reference := interval,
       k := 5,
       max_distance := 100000,
       stranded := true,
       signed := true
   )

   -- Standalone query with literal reference
   SELECT * FROM NEAREST(genes, reference := 'chr1:1000000-1001000', k := 5)

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

.. code-block:: sql

   SELECT
       peaks.name AS peak,
       nearest.name AS gene,
       nearest.distance
   FROM peaks
   CROSS JOIN LATERAL NEAREST(genes, reference := peaks.interval, k := 3) AS nearest
   ORDER BY peaks.name, nearest.distance

**Standalone Query:**

Find 5 nearest genes to a specific genomic location:

.. code-block:: sql

   SELECT gene_name, distance
   FROM NEAREST(genes, reference := 'chr1:1000000-1001000', k := 5)
   ORDER BY distance

**Distance-Constrained Search:**

Find nearest features within 100kb:

.. code-block:: sql

   SELECT
       peaks.name,
       nearest.name AS gene,
       nearest.distance
   FROM peaks
   CROSS JOIN LATERAL NEAREST(
       genes,
       reference := peaks.interval,
       k := 5,
       max_distance := 100000
   ) AS nearest
   ORDER BY peaks.name, nearest.distance

**Strand-Specific Nearest Neighbors:**

Find nearest same-strand features:

.. code-block:: sql

   SELECT
       peaks.name,
       nearest.name AS gene,
       nearest.strand,
       nearest.distance
   FROM peaks
   CROSS JOIN LATERAL NEAREST(
       genes,
       reference := peaks.interval,
       k := 3,
       stranded := true
   ) AS nearest
   ORDER BY peaks.name, nearest.distance

**Directional (Upstream/Downstream) Queries:**

Find upstream features using signed distances:

.. code-block:: sql

   -- Upstream features have negative distances
   SELECT
       peaks.name,
       nearest.name AS gene,
       nearest.distance
   FROM peaks
   CROSS JOIN LATERAL NEAREST(
       genes,
       reference := peaks.interval,
       k := 10,
       signed := true
   ) AS nearest
   WHERE nearest.distance < 0
   ORDER BY peaks.name, nearest.distance DESC

.. code-block:: sql

   -- Downstream features have positive distances
   SELECT
       peaks.name,
       nearest.name AS gene,
       nearest.distance
   FROM peaks
   CROSS JOIN LATERAL NEAREST(
       genes,
       reference := peaks.interval,
       k := 10,
       signed := true
   ) AS nearest
   WHERE nearest.distance > 0
   ORDER BY peaks.name, nearest.distance

**Combined Parameters:**

Find nearby same-strand features within distance constraints:

.. code-block:: sql

   SELECT
       peaks.name,
       nearest.name AS gene,
       nearest.distance
   FROM peaks
   CROSS JOIN LATERAL NEAREST(
       genes,
       reference := peaks.interval,
       k := 5,
       max_distance := 50000,
       stranded := true,
       signed := true
   ) AS nearest
   WHERE nearest.distance BETWEEN -10000 AND 10000
   ORDER BY peaks.name, ABS(nearest.distance)

Target support
~~~~~~~~~~~~~~

A correlated ``NEAREST`` (its reference is an outer-row column) runs on lateral-capable engines — DuckDB and the generic target — via a correlated ``LATERAL`` subquery, and on Apache DataFusion, which has no correlated-``LATERAL`` physical plan, via a decorrelated window-function rewrite. For an **explicitly-projected** query (one that selects named columns, e.g. ``SELECT a.start, b.start, b.distance``) the two forms return identical results: the ``(start, end)`` tiebreaker orders rows tied at the k-th distance the same way on every engine, deterministically whenever ``(start, end)`` distinguishes the tied candidates. A standalone ``NEAREST`` with a literal reference is uncorrelated and uses the same ordered, limited subquery on every target.

.. note::

   **Known limitation —** ``SELECT *`` **/** ``SELECT b.*`` **over a correlated NEAREST on DataFusion.** The decorrelated window-function rewrite needs its reference-key and rank columns (``__giql_x_rk_*``, ``__giql_x_rn``) visible on the rewritten join, so a ``SELECT *`` or ``SELECT b.*`` over a correlated NEAREST exposes those reserved internal columns on DataFusion — a different output schema than the LATERAL form emits on DuckDB. The cross-target identity claim above therefore holds for **explicitly-projected** queries only. Projecting named columns avoids the leak entirely. A query-level wrapper that projects the reserved columns away on the DataFusion path is tracked by `#160 <https://github.com/abdenlab/giql/issues/160>`_ (it depends on the query-level expander seam from #146).

Notes
~~~~~

- **Chromosome pre-filtering**: NEAREST automatically filters by chromosome for efficiency
- **Use max_distance**: Specifying a maximum distance reduces the search space significantly
- **Limit k**: Only request as many neighbors as you actually need

Related Operators
~~~~~~~~~~~~~~~~~

- :ref:`DISTANCE <distance-operator>` - Calculate distance between specific pairs
- :ref:`INTERSECTS <intersects-operator>` - Find overlapping features (distance = 0)
