Aggregation Operators
=====================

Aggregation operators combine and cluster genomic intervals. These operators are
essential for reducing complex interval data into summarized regions, such as
merging overlapping peaks or identifying clusters of related features.

.. contents::
   :local:
   :depth: 2

.. _cluster-operator:

CLUSTER
-------

Assign cluster IDs to overlapping or nearby genomic intervals.

Description
~~~~~~~~~~~

The ``CLUSTER`` operator assigns a unique cluster identifier to groups of intervals
that overlap or are within a specified distance of each other. Intervals in the same
cluster share a common cluster ID, while non-overlapping intervals receive different
IDs.

This is useful for:

- Grouping overlapping features
- Identifying regions of high feature density
- Preparing data for downstream merge operations

Syntax
~~~~~~

.. code-block:: sql

   -- Basic clustering (overlapping intervals)
   CLUSTER(interval) AS cluster_id

   -- Clustering with distance parameter
   CLUSTER(interval, distance) AS cluster_id

   -- Strand-specific clustering
   CLUSTER(interval, stranded=true) AS cluster_id

   -- Combined parameters
   CLUSTER(interval, distance, stranded=true) AS cluster_id

Parameters
~~~~~~~~~~

**interval**
   A genomic column registered with the engine.

**distance** *(optional)*
   Maximum gap between intervals to consider them part of the same cluster.
   Default: ``0`` (only overlapping intervals are clustered).

**stranded** *(optional)*
   When ``true``, only cluster intervals on the same strand. Default: ``false``.

Return Value
~~~~~~~~~~~~

Integer cluster ID. Intervals in the same cluster have the same ID.
IDs are assigned per-chromosome (and per-strand if ``stranded=true``).

Examples
~~~~~~~~

**Basic Clustering:**

Assign cluster IDs to overlapping intervals:

.. code-block:: python

   cursor = engine.execute("""
       SELECT
           *,
           CLUSTER(interval) AS cluster_id
       FROM features
       ORDER BY chromosome, start_pos
   """)

**Distance-Based Clustering:**

Cluster intervals within 1000bp of each other:

.. code-block:: python

   cursor = engine.execute("""
       SELECT
           *,
           CLUSTER(interval, 1000) AS cluster_id
       FROM features
       ORDER BY chromosome, start_pos
   """)

**Strand-Specific Clustering:**

Cluster intervals separately by strand:

.. code-block:: python

   cursor = engine.execute("""
       SELECT
           *,
           CLUSTER(interval, stranded=true) AS cluster_id
       FROM features
       ORDER BY chromosome, strand, start_pos
   """)

**Analyze Cluster Statistics:**

Count features per cluster:

.. code-block:: python

   cursor = engine.execute("""
       WITH clustered AS (
           SELECT
               *,
               CLUSTER(interval) AS cluster_id
           FROM features
       )
       SELECT
           chromosome,
           cluster_id,
           COUNT(*) AS feature_count,
           MIN(start_pos) AS cluster_start,
           MAX(end_pos) AS cluster_end
       FROM clustered
       GROUP BY chromosome, cluster_id
       ORDER BY chromosome, cluster_start
   """)

**Filter by Cluster Size:**

Find regions with multiple overlapping features:

.. code-block:: python

   cursor = engine.execute("""
       WITH clustered AS (
           SELECT
               *,
               CLUSTER(interval) AS cluster_id
           FROM features
       ),
       cluster_sizes AS (
           SELECT cluster_id, COUNT(*) AS size
           FROM clustered
           GROUP BY cluster_id
       )
       SELECT c.*
       FROM clustered c
       INNER JOIN cluster_sizes s ON c.cluster_id = s.cluster_id
       WHERE s.size >= 3
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
     - Efficient window function implementation
   * - SQLite
     - Full
     -
   * - PostgreSQL
     - Planned
     -

Performance Notes
~~~~~~~~~~~~~~~~~

- Data should be sorted by chromosome and position for efficient clustering
- For large datasets, consider partitioning by chromosome
- Cluster IDs are computed using window functions, which scale well

Related Operators
~~~~~~~~~~~~~~~~~

- :ref:`MERGE <merge-operator>` - Combine clustered intervals into single regions
- :ref:`INTERSECTS <intersects-operator>` - Test for overlap between specific pairs

----

.. _merge-operator:

MERGE
-----

Combine overlapping genomic intervals into unified regions.

Description
~~~~~~~~~~~

The ``MERGE`` operator combines overlapping (or nearby) intervals into single,
non-overlapping regions. This is useful for:

- Creating consensus regions from overlapping features
- Reducing redundant annotations
- Calculating total coverage

The operator works as an aggregate function, returning one row per merged region
with the unified coordinates.

Syntax
~~~~~~

.. code-block:: sql

   -- Basic merge
   SELECT MERGE(interval) FROM features

   -- Merge with distance parameter
   SELECT MERGE(interval, distance) FROM features

   -- Strand-specific merge
   SELECT MERGE(interval, stranded=true) FROM features

   -- Merge with additional aggregations
   SELECT
       MERGE(interval),
       COUNT(*) AS feature_count,
       AVG(score) AS avg_score
   FROM features

Parameters
~~~~~~~~~~

**interval**
   A genomic column registered with the engine.

**distance** *(optional)*
   Maximum gap between intervals to merge. Default: ``0`` (only overlapping
   intervals are merged).

**stranded** *(optional)*
   When ``true``, merge intervals separately by strand. Default: ``false``.

Return Value
~~~~~~~~~~~~

Returns merged interval coordinates:

- ``chromosome`` - Chromosome of the merged region
- ``start_pos`` - Start position of the merged region
- ``end_pos`` - End position of the merged region
- ``strand`` - Strand (if ``stranded=true``)

Examples
~~~~~~~~

**Basic Merge:**

Merge all overlapping intervals:

.. code-block:: python

   cursor = engine.execute("""
       SELECT MERGE(interval)
       FROM features
   """)

   # Returns: chromosome, start_pos, end_pos for each merged region

**Distance-Based Merge:**

Merge intervals within 1000bp of each other:

.. code-block:: python

   cursor = engine.execute("""
       SELECT MERGE(interval, 1000)
       FROM features
   """)

**Strand-Specific Merge:**

Merge intervals separately by strand:

.. code-block:: python

   cursor = engine.execute("""
       SELECT MERGE(interval, stranded=true)
       FROM features
   """)

**Merge with Feature Count:**

Count how many features were merged into each region:

.. code-block:: python

   cursor = engine.execute("""
       SELECT
           MERGE(interval),
           COUNT(*) AS feature_count
       FROM features
   """)

**Merge with Aggregations:**

Calculate statistics for merged regions:

.. code-block:: python

   cursor = engine.execute("""
       SELECT
           MERGE(interval),
           COUNT(*) AS feature_count,
           AVG(score) AS avg_score,
           MAX(score) AS max_score
       FROM features
   """)

**Collect Merged Feature Names:**

List the names of features that were merged:

.. code-block:: python

   cursor = engine.execute("""
       SELECT
           MERGE(interval),
           STRING_AGG(name, ',') AS feature_names
       FROM features
   """)

**Merge by Chromosome:**

Process each chromosome separately (explicit grouping):

.. code-block:: python

   cursor = engine.execute("""
       SELECT
           chromosome,
           MERGE(interval),
           COUNT(*) AS feature_count
       FROM features
       GROUP BY chromosome
       ORDER BY chromosome
   """)

**Calculate Total Coverage:**

Calculate the total base pairs covered after merging:

.. code-block:: python

   cursor = engine.execute("""
       WITH merged AS (
           SELECT MERGE(interval) AS merged_pos
           FROM features
       )
       SELECT SUM(end_pos - start_pos) AS total_coverage
       FROM merged
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

- MERGE is an aggregate operation that processes all matching rows
- For very large datasets, consider filtering by chromosome first
- The operation sorts data internally, so pre-sorting is not required

Related Operators
~~~~~~~~~~~~~~~~~

- :ref:`CLUSTER <cluster-operator>` - Assign cluster IDs without merging
- :ref:`INTERSECTS <intersects-operator>` - Test for overlap between specific pairs
