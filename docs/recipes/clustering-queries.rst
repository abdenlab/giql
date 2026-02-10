Clustering and Merging Queries
==============================

This section covers patterns for clustering overlapping intervals and merging
them into unified regions using GIQL's aggregation operators.

.. contents::
   :local:
   :depth: 2

Basic Clustering
----------------

Assign Cluster IDs
~~~~~~~~~~~~~~~~~~

Assign unique cluster IDs to groups of overlapping intervals:

.. code-block:: sql

   SELECT
       *,
       CLUSTER(interval) AS cluster_id
   FROM features
   ORDER BY chrom, start

**Use case:** Group overlapping peaks or annotations for downstream analysis.

View Cluster Assignments
~~~~~~~~~~~~~~~~~~~~~~~~

See which features belong to which cluster:

.. code-block:: sql

   SELECT
       cluster_id,
       chrom,
       name,
       start,
       end
   FROM (
       SELECT *, CLUSTER(interval) AS cluster_id
       FROM features
   )
   ORDER BY cluster_id, start

**Use case:** Inspect clustering results to understand feature groupings.

Distance-Based Clustering
-------------------------

Cluster with Gap Tolerance
~~~~~~~~~~~~~~~~~~~~~~~~~~

Cluster intervals that are within a specified distance of each other:

.. code-block:: sql

   SELECT
       *,
       CLUSTER(interval, 1000) AS cluster_id
   FROM features
   ORDER BY chrom, start

**Use case:** Group nearby features even if they don't directly overlap
(e.g., cluster peaks within 1kb of each other).

Variable Distance Thresholds
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Experiment with different clustering distances:

.. code-block:: sql

   -- Tight clustering (overlapping only)
   SELECT *, CLUSTER(interval, 0) AS tight_cluster FROM features

   -- Medium clustering (within 500bp)
   SELECT *, CLUSTER(interval, 500) AS medium_cluster FROM features

   -- Loose clustering (within 5kb)
   SELECT *, CLUSTER(interval, 5000) AS loose_cluster FROM features

**Use case:** Compare clustering at different resolutions for sensitivity analysis.

Strand-Specific Clustering
--------------------------

Cluster by Strand
~~~~~~~~~~~~~~~~~

Cluster intervals separately for each strand:

.. code-block:: sql

   SELECT
       *,
       CLUSTER(interval, stranded=true) AS cluster_id
   FROM features
   ORDER BY chrom, strand, start

**Use case:** Maintain strand separation when clustering transcripts or
strand-specific regulatory elements.

Strand-Specific with Distance
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Combine strand awareness with distance tolerance:

.. code-block:: sql

   SELECT
       *,
       CLUSTER(interval, 1000, stranded=true) AS cluster_id
   FROM features
   ORDER BY chrom, strand, start

**Use case:** Cluster nearby same-strand features while keeping opposite
strands separate.

Cluster Statistics
------------------

Count Features per Cluster
~~~~~~~~~~~~~~~~~~~~~~~~~~

Calculate how many features are in each cluster:

.. code-block:: sql

   WITH clustered AS (
       SELECT *, CLUSTER(interval) AS cluster_id
       FROM features
   )
   SELECT
       cluster_id,
       chrom,
       COUNT(*) AS feature_count,
       MIN(start) AS cluster_start,
       MAX(end) AS cluster_end
   FROM clustered
   GROUP BY cluster_id, chrom
   ORDER BY chrom, cluster_start

**Use case:** Identify cluster sizes and boundaries.

Filter by Cluster Size
~~~~~~~~~~~~~~~~~~~~~~

Find clusters with a minimum number of features:

.. code-block:: sql

   WITH clustered AS (
       SELECT *, CLUSTER(interval) AS cluster_id
       FROM features
   ),
   cluster_sizes AS (
       SELECT cluster_id, COUNT(*) AS size
       FROM clustered
       GROUP BY cluster_id
   )
   SELECT c.*
   FROM clustered c
   JOIN cluster_sizes s ON c.cluster_id = s.cluster_id
   WHERE s.size >= 3
   ORDER BY c.cluster_id, c.start

**Use case:** Focus on regions with multiple overlapping features (hotspots).

Cluster Summary Statistics
~~~~~~~~~~~~~~~~~~~~~~~~~~

Calculate statistics for each cluster:

.. code-block:: sql

   WITH clustered AS (
       SELECT *, CLUSTER(interval) AS cluster_id
       FROM features
   )
   SELECT
       cluster_id,
       chrom,
       COUNT(*) AS feature_count,
       MIN(start) AS cluster_start,
       MAX(end) AS cluster_end,
       MAX(end) - MIN(start) AS cluster_span,
       AVG(score) AS avg_score,
       MAX(score) AS max_score
   FROM clustered
   GROUP BY cluster_id, chrom
   ORDER BY feature_count DESC

**Use case:** Rank clusters by size, span, or aggregate scores.

Basic Merging
-------------

Merge Overlapping Intervals
~~~~~~~~~~~~~~~~~~~~~~~~~~~

Combine overlapping intervals into unified regions:

.. code-block:: sql

   SELECT MERGE(interval)
   FROM features

**Use case:** Create non-overlapping consensus regions from redundant annotations.

Merge with Distance
~~~~~~~~~~~~~~~~~~~

Merge intervals within a specified distance:

.. code-block:: sql

   SELECT MERGE(interval, 1000)
   FROM features

**Use case:** Create broader regions by joining nearby features.

Strand-Specific Merge
~~~~~~~~~~~~~~~~~~~~~

Merge intervals separately by strand:

.. code-block:: sql

   SELECT MERGE(interval, stranded=true)
   FROM features

**Use case:** Create strand-aware consensus regions.

Merge with Aggregations
-----------------------

Count Merged Features
~~~~~~~~~~~~~~~~~~~~~

Track how many features were merged into each region:

.. code-block:: sql

   SELECT
       MERGE(interval),
       COUNT(*) AS feature_count
   FROM features

**Use case:** Understand the complexity of each merged region.

Aggregate Scores
~~~~~~~~~~~~~~~~

Calculate statistics for merged regions:

.. code-block:: sql

   SELECT
       MERGE(interval),
       COUNT(*) AS feature_count,
       AVG(score) AS avg_score,
       MAX(score) AS max_score,
       SUM(score) AS total_score
   FROM features

**Use case:** Summarize signal intensity across merged regions.

Collect Feature Names
~~~~~~~~~~~~~~~~~~~~~

List the names of features that were merged:

.. code-block:: sql

   SELECT
       MERGE(interval),
       STRING_AGG(name, ',') AS merged_features
   FROM features

**Use case:** Track provenance of merged regions.

Coverage Calculations
---------------------

Total Base Pair Coverage
~~~~~~~~~~~~~~~~~~~~~~~~

Calculate total genomic coverage after merging:

.. code-block:: sql

   WITH merged AS (
       SELECT MERGE(interval)
       FROM features
   )
   SELECT SUM(end - start) AS total_coverage_bp
   FROM merged

**Use case:** Calculate the total genome fraction covered by features.

Coverage per Chromosome
~~~~~~~~~~~~~~~~~~~~~~~

Calculate coverage for each chromosome:

.. code-block:: sql

   WITH merged AS (
       SELECT MERGE(interval)
       FROM features
   )
   SELECT
       chrom,
       COUNT(*) AS region_count,
       SUM(end - start) AS coverage_bp
   FROM merged
   GROUP BY chrom
   ORDER BY chrom

**Use case:** Compare feature density across chromosomes.

Coverage Reduction
~~~~~~~~~~~~~~~~~~

Compare raw vs merged coverage:

.. code-block:: sql

   WITH raw_stats AS (
       SELECT
           COUNT(*) AS raw_count,
           SUM(end - start) AS raw_bp
       FROM features
   ),
   merged_stats AS (
       SELECT
           COUNT(*) AS merged_count,
           SUM(end - start) AS merged_bp
       FROM (SELECT MERGE(interval) FROM features)
   )
   SELECT
       raw_count,
       merged_count,
       raw_bp,
       merged_bp,
       ROUND(100.0 * merged_bp / raw_bp, 2) AS coverage_retained_pct
   FROM raw_stats, merged_stats

**Use case:** Quantify the redundancy in your feature set.

Advanced Patterns
-----------------

Cluster Then Merge
~~~~~~~~~~~~~~~~~~

First cluster features, then analyze each cluster:

.. code-block:: sql

   WITH clustered AS (
       SELECT *, CLUSTER(interval) AS cluster_id
       FROM features
   )
   SELECT
       cluster_id,
       MIN(chrom) AS chrom,
       MIN(start) AS start,
       MAX(end) AS end,
       COUNT(*) AS feature_count,
       STRING_AGG(name, ',') AS features
   FROM clustered
   GROUP BY cluster_id
   ORDER BY chrom, start

**Use case:** Alternative to MERGE that preserves cluster identifiers.

Hierarchical Clustering
~~~~~~~~~~~~~~~~~~~~~~~

Apply multiple clustering levels:

.. code-block:: sql

   WITH level1 AS (
       SELECT *, CLUSTER(interval, 0) AS cluster_l1
       FROM features
   ),
   level2 AS (
       SELECT *, CLUSTER(interval, 1000) AS cluster_l2
       FROM level1
   )
   SELECT
       cluster_l1,
       cluster_l2,
       chrom,
       name,
       start,
       end
   FROM level2
   ORDER BY cluster_l2, cluster_l1, start

**Use case:** Analyze feature relationships at multiple scales.
