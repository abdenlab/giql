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
   CLUSTER(interval, stranded := true) AS cluster_id

   -- Predicate-gated clustering (run-length encoding on a column)
   CLUSTER(interval, predicate := depth = PREV(depth)) AS cluster_id

   -- Combined parameters
   CLUSTER(interval, distance, stranded := true) AS cluster_id

Parameters
~~~~~~~~~~

**interval**
   A genomic column.

**distance** *(optional)*
   Maximum gap between intervals to consider them part of the same cluster.
   Default: ``0`` (only overlapping intervals are clustered).

**stranded** *(optional)*
   When ``true``, only cluster intervals on the same strand. Default: ``false``.

**predicate** *(optional)*
   A boolean expression evaluated between each interval and its sorted
   predecessor. When supplied, the cluster-boundary condition becomes
   **adjacent AND predicate**: an interval stays in the current cluster only
   when it is within ``distance`` of its predecessor *and* the predicate holds
   between the two. A change in the predicate forces a new cluster, so an
   equality predicate yields a run-length encoding of the input sequence.
   Omitting the predicate preserves the default adjacency-only behavior.

   Bare column references resolve to the *current* interval; the predecessor's
   value of a column is referenced with ``PREV(column)``
   (e.g. ``depth = PREV(depth)``). The predicate composes with ``distance`` and
   ``stranded`` and is evaluated under the operator's existing per-chromosome
   (and per-strand) partition and start-position order.

   Two constraints apply:

   - **References existing columns only.** The predicate *gates* merging on
     columns already present on the input rows; it does not synthesize a
     statistic. Coverage depth, for example, must already be a column on the
     rows (typically produced upstream by :ref:`DISJOIN <disjoin-operator>` and
     aggregation).
   - **Pairwise only, with single-linkage drift.** The predicate compares each
     interval to its immediate sorted predecessor (everything ``LAG`` can
     express). Whole-cluster conditions are out of scope. When the predicate is
     not an equivalence relation (e.g. ``ABS(score - PREV(score)) < 5``),
     consecutive pairs may each satisfy it while the cluster's extremes do not
     — the same single-linkage behavior that ``distance``-based clustering
     already exhibits.

Return Value
~~~~~~~~~~~~

Integer cluster ID. Intervals in the same cluster have the same ID.
IDs are assigned per-chromosome (and per-strand if ``stranded := true``).

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
       CLUSTER(interval, stranded := true) AS cluster_id
   FROM features
   ORDER BY chrom, strand, start

**Predicate-Gated Clustering:**

Cut adjacent intervals into clusters wherever a column's value changes
(run-length encoding). ``PREV(column)`` references the predecessor row's value:

.. code-block:: sql

   SELECT
       *,
       CLUSTER(interval, predicate := depth = PREV(depth)) AS cluster_id
   FROM features
   ORDER BY chrom, start

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

.. note::

   **Synthesized flag hidden under** ``SELECT *``. A ``SELECT *, CLUSTER(...)`` query materializes an internal ``__giql_is_new_cluster`` flag in a subquery, so the outer star is emitted as ``SELECT * EXCEPT (__giql_is_new_cluster)`` to keep that helper column out of the result (#184). ``* EXCEPT`` is a DataFusion-family extension: the generic and ``datafusion`` dialects emit it, while ``duckdb`` spells the exclusion ``EXCLUDE``. Transpile with ``dialect="duckdb"`` to execute on DuckDB — the portable generic ``* EXCEPT`` form is not DuckDB-runnable. A qualified ``SELECT t.*, CLUSTER(...)`` receives the same treatment: because ``CLUSTER`` runs over a single relation, the qualifier is dropped and the outer star is emitted as the same bare ``* EXCEPT (__giql_is_new_cluster)`` (#185). An explicitly-projected ``CLUSTER`` (no star) surfaces no helper column and needs no exclusion.

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
   SELECT MERGE(interval, stranded := true) FROM features

   -- Predicate-gated merge (merge only equal-valued adjacent runs)
   SELECT MERGE(interval, predicate := depth = PREV(depth)) FROM features

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

**predicate** *(optional)*
   A boolean expression that further restricts which adjacent intervals are
   merged. ``MERGE`` decomposes into :ref:`CLUSTER <cluster-operator>` plus a
   ``GROUP BY`` over the cluster id, so it inherits predicate-aware boundaries
   directly — see the :ref:`CLUSTER predicate <cluster-operator>` description
   for the full semantics, the ``PREV(column)`` convention, the
   references-existing-columns-only constraint, and the pairwise-only /
   single-linkage caveat. Omitting the predicate preserves the default
   adjacency-only merge.

Return Value
~~~~~~~~~~~~

Returns merged interval coordinates:

- ``chrom`` - Chromosome of the merged region
- ``start`` - Start position of the merged region
- ``end`` - End position of the merged region
- ``strand`` - Strand (if ``stranded := true``)

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

   SELECT MERGE(interval, stranded := true)
   FROM features

**Predicate-Gated Merge (coverage depth):**

Merge only adjacent intervals that share the same coverage depth, reconstructing
a re-clustered, depth-segmented partition from per-breakpoint segments produced
by :ref:`DISJOIN <disjoin-operator>` and aggregation:

.. code-block:: sql

   SELECT MERGE(interval, predicate := depth = PREV(depth))
   FROM (
       SELECT disjoin_chrom AS chrom,
              disjoin_start AS start,
              disjoin_end AS end,
              COUNT(*) AS depth
       FROM DISJOIN(features)
       GROUP BY disjoin_chrom, disjoin_start, disjoin_end
   ) AS segments

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

.. note::

   CLUSTER and MERGE cannot be combined in a single ``SELECT`` — MERGE
   aggregates rows away while CLUSTER is a per-row window over those same rows,
   so no single query expresses both. Transpiling ``SELECT MERGE(interval),
   CLUSTER(interval) FROM features`` raises a ``ValueError``. Use them in
   separate queries instead — for example, CLUSTER over a subquery, or MERGE
   over one.

.. note::

   MERGE cannot be projected alongside a star. MERGE aggregates rows into one
   row per merged region, so a ``SELECT *, MERGE(...)`` or ``SELECT t.*,
   MERGE(...)`` has no coherent per-row meaning — the star names the
   pre-aggregation input columns of a relation MERGE has already collapsed and
   grouped away. Transpiling either shape raises a ``ValueError`` rather than
   emitting non-executable SQL (a bare ``*`` re-surfaces non-grouped columns
   under the synthesized ``GROUP BY``; a qualified ``rel.*`` dangles an alias the
   aggregation no longer exposes). Drop the star, or project only grouping
   columns and aggregates (e.g. ``COUNT(*)``) alongside ``MERGE`` — not raw input
   columns, which are neither grouped nor aggregated. This is unlike
   :ref:`CLUSTER <cluster-operator>`, a per-row window over which a star *is*
   meaningful and supported.

Related Operators
~~~~~~~~~~~~~~~~~

- :ref:`CLUSTER <cluster-operator>` - Assign cluster IDs without merging
- :ref:`INTERSECTS <intersects-operator>` - Test for overlap between specific pairs
