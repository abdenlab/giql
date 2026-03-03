Aggregation
===========

Aggregation operators combine and cluster genomic intervals. These operators are
essential for reducing complex interval data into summarized regions, such as
merging overlapping peaks or identifying clusters of related features.

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
   A genomic column.

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

.. code-block:: sql

   SELECT
       *,
       CLUSTER(interval) AS cluster_id
   FROM features
   ORDER BY chrom, start

**Distance-Based Clustering:**

Cluster intervals within 1000bp of each other:

.. code-block:: sql

   SELECT
       *,
       CLUSTER(interval, 1000) AS cluster_id
   FROM features
   ORDER BY chrom, start

**Strand-Specific Clustering:**

Cluster intervals separately by strand:

.. code-block:: sql

   SELECT
       *,
       CLUSTER(interval, stranded=true) AS cluster_id
   FROM features
   ORDER BY chrom, strand, start

**Analyze Cluster Statistics:**

Count features per cluster:

.. code-block:: sql

   WITH clustered AS (
       SELECT
           *,
           CLUSTER(interval) AS cluster_id
       FROM features
   )
   SELECT
       chrom,
       cluster_id,
       COUNT(*) AS feature_count,
       MIN(start) AS cluster_start,
       MAX(end) AS cluster_end
   FROM clustered
   GROUP BY chrom, cluster_id
   ORDER BY chrom, cluster_start

**Filter by Cluster Size:**

Find regions with multiple overlapping features:

.. code-block:: sql

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
   A genomic column.

**distance** *(optional)*
   Maximum gap between intervals to merge. Default: ``0`` (only overlapping
   intervals are merged).

**stranded** *(optional)*
   When ``true``, merge intervals separately by strand. Default: ``false``.

Return Value
~~~~~~~~~~~~

Returns merged interval coordinates:

- ``chrom`` - Chromosome of the merged region
- ``start`` - Start position of the merged region
- ``end`` - End position of the merged region
- ``strand`` - Strand (if ``stranded=true``)

Examples
~~~~~~~~

**Basic Merge:**

Merge all overlapping intervals:

.. code-block:: sql

   SELECT MERGE(interval)
   FROM features

   -- Returns: chrom, start, end for each merged region

**Distance-Based Merge:**

Merge intervals within 1000bp of each other:

.. code-block:: sql

   SELECT MERGE(interval, 1000)
   FROM features

**Strand-Specific Merge:**

Merge intervals separately by strand:

.. code-block:: sql

   SELECT MERGE(interval, stranded=true)
   FROM features

**Merge with Feature Count:**

Count how many features were merged into each region:

.. code-block:: sql

   SELECT
       MERGE(interval),
       COUNT(*) AS feature_count
   FROM features

**Merge with Aggregations:**

Calculate statistics for merged regions:

.. code-block:: sql

   SELECT
       MERGE(interval),
       COUNT(*) AS feature_count,
       AVG(score) AS avg_score,
       MAX(score) AS max_score
   FROM features

**Collect Merged Feature Names:**

List the names of features that were merged:

.. code-block:: sql

   SELECT
       MERGE(interval),
       STRING_AGG(name, ',') AS feature_names
   FROM features

**Merge by Chromosome:**

Process each chromosome separately (explicit grouping):

.. code-block:: sql

   SELECT
       chrom,
       MERGE(interval),
       COUNT(*) AS feature_count
   FROM features
   GROUP BY chrom
   ORDER BY chrom

**Calculate Total Coverage:**

Calculate the total base pairs covered after merging:

.. code-block:: sql

   WITH merged AS (
       SELECT MERGE(interval) AS merged_pos
       FROM features
   )
   SELECT SUM(end - start) AS total_coverage
   FROM merged

Notes
~~~~~

- MERGE is an aggregate operation that processes all matching rows
- The operation sorts data internally, so pre-sorting is not required

Related Operators
~~~~~~~~~~~~~~~~~

- :ref:`CLUSTER <cluster-operator>` - Assign cluster IDs without merging
- :ref:`INTERSECTS <intersects-operator>` - Test for overlap between specific pairs
