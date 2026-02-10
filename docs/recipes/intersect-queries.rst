Intersection Queries
====================

This section covers common patterns for finding overlapping genomic features
using GIQL's spatial operators.

.. contents::
   :local:
   :depth: 2

Finding Overlapping Features
----------------------------

Basic Overlap Query
~~~~~~~~~~~~~~~~~~~

Find all features in table A that overlap with any feature in table B:

.. code-block:: sql

   SELECT DISTINCT a.*
   FROM features_a a, features_b b
   WHERE a.interval INTERSECTS b.interval

**Use case:** Identify variants that fall within gene regions.

Get All Overlap Pairs
~~~~~~~~~~~~~~~~~~~~~

Return every pair of overlapping features (may produce duplicates if one
feature overlaps multiple others):

.. code-block:: sql

   SELECT a.*, b.*
   FROM features_a a, features_b b
   WHERE a.interval INTERSECTS b.interval

**Use case:** Generate a full overlap matrix for downstream analysis.

Query Against a Specific Region
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Find features overlapping a literal genomic range:

.. code-block:: sql

   SELECT * FROM variants
   WHERE interval INTERSECTS 'chr1:1000000-2000000'

**Use case:** Extract all data for a specific chromosomal region.

Filtering by Overlap
--------------------

Excluding Overlaps
~~~~~~~~~~~~~~~~~~

Find features in A that do NOT overlap with any feature in B:

.. code-block:: sql

   SELECT a.*
   FROM features_a a
   LEFT JOIN features_b b ON a.interval INTERSECTS b.interval
   WHERE b.chrom IS NULL

**Use case:** Find regulatory regions that don't overlap with known genes,
or identify variants outside of exonic regions.

Features with Any Overlap (Unique)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Return each feature from A only once, regardless of how many B features it overlaps:

.. code-block:: sql

   SELECT DISTINCT a.*
   FROM features_a a
   INNER JOIN features_b b ON a.interval INTERSECTS b.interval

**Use case:** Get a deduplicated list of features that have at least one overlap.

Counting Overlaps
-----------------

Count Overlapping Features
~~~~~~~~~~~~~~~~~~~~~~~~~~

Count how many B features each A feature overlaps:

.. code-block:: sql

   SELECT a.*, COUNT(b.name) AS overlap_count
   FROM features_a a
   LEFT JOIN features_b b ON a.interval INTERSECTS b.interval
   GROUP BY a.chrom, a.start, a.end, a.name, a.score, a.strand

**Use case:** Calculate how many enhancers each gene overlaps with,
or count variants per feature.

Filter by Overlap Count
~~~~~~~~~~~~~~~~~~~~~~~

Find features that overlap at least N other features:

.. code-block:: sql

   SELECT a.*
   FROM features_a a
   INNER JOIN features_b b ON a.interval INTERSECTS b.interval
   GROUP BY a.chrom, a.start, a.end, a.name, a.score, a.strand
   HAVING COUNT(*) >= 3

**Use case:** Identify hotspot regions with high feature density.

Strand-Specific Operations
--------------------------

Same-Strand Overlaps
~~~~~~~~~~~~~~~~~~~~

Find overlapping features on the same strand:

.. code-block:: sql

   SELECT a.*, b.name AS b_name
   FROM features_a a, features_b b
   WHERE a.interval INTERSECTS b.interval
     AND a.strand = b.strand

**Use case:** Find sense-strand overlaps for transcript analysis.

Opposite-Strand Overlaps
~~~~~~~~~~~~~~~~~~~~~~~~

Find overlapping features on opposite strands:

.. code-block:: sql

   SELECT a.*, b.name AS b_name
   FROM features_a a, features_b b
   WHERE a.interval INTERSECTS b.interval
     AND a.strand != b.strand
     AND a.strand IN ('+', '-')
     AND b.strand IN ('+', '-')

**Use case:** Identify antisense overlaps or convergent transcription.

Overlap Fraction Requirements
-----------------------------

Minimum Overlap Fraction of A
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Find overlaps where at least 50% of feature A is covered:

.. code-block:: sql

   SELECT a.*
   FROM features_a a, features_b b
   WHERE a.interval INTERSECTS b.interval
     AND (
         LEAST(a.end, b.end) - GREATEST(a.start, b.start)
     ) >= 0.5 * (a.end - a.start)

**Use case:** Ensure substantial overlap rather than just touching edges.

Minimum Overlap Fraction of B
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Find overlaps where at least 50% of feature B is covered:

.. code-block:: sql

   SELECT a.*
   FROM features_a a, features_b b
   WHERE a.interval INTERSECTS b.interval
     AND (
         LEAST(a.end, b.end) - GREATEST(a.start, b.start)
     ) >= 0.5 * (b.end - b.start)

**Use case:** Find features that substantially cover smaller annotations.

Reciprocal Overlap
~~~~~~~~~~~~~~~~~~

Require both features to have at least 50% mutual overlap:

.. code-block:: sql

   WITH overlap_calcs AS (
       SELECT
           a.*,
           b.name AS b_name,
           (LEAST(a.end, b.end) - GREATEST(a.start, b.start)) AS overlap_bp,
           (a.end - a.start) AS a_length,
           (b.end - b.start) AS b_length
       FROM features_a a, features_b b
       WHERE a.interval INTERSECTS b.interval
   )
   SELECT *
   FROM overlap_calcs
   WHERE overlap_bp >= 0.5 * a_length
     AND overlap_bp >= 0.5 * b_length

**Use case:** Find high-confidence overlaps where features mutually cover each other.

Join Patterns
-------------

Left Outer Join
~~~~~~~~~~~~~~~

Report all features from A, with B information where available:

.. code-block:: sql

   SELECT a.*, b.name AS overlapping_feature
   FROM features_a a
   LEFT JOIN features_b b ON a.interval INTERSECTS b.interval

**Use case:** Annotate features with overlap information while keeping all records.

Calculate Overlap Amount
~~~~~~~~~~~~~~~~~~~~~~~~

Return the overlap size in base pairs:

.. code-block:: sql

   SELECT
       a.*,
       b.name AS b_name,
       (LEAST(a.end, b.end) - GREATEST(a.start, b.start)) AS overlap_bp
   FROM features_a a, features_b b
   WHERE a.interval INTERSECTS b.interval

**Use case:** Quantify the extent of each overlap.

Overlap with NULL Handling
~~~~~~~~~~~~~~~~~~~~~~~~~~

Report overlap amount for all A features, with 0 for non-overlapping:

.. code-block:: sql

   SELECT
       a.*,
       b.name AS b_name,
       CASE
           WHEN b.chrom IS NULL THEN 0
           ELSE LEAST(a.end, b.end) - GREATEST(a.start, b.start)
       END AS overlap_bp
   FROM features_a a
   LEFT JOIN features_b b ON a.interval INTERSECTS b.interval

**Use case:** Create a complete overlap report including non-overlapping features.

Multi-Table Operations
----------------------

Union Multiple Sources
~~~~~~~~~~~~~~~~~~~~~~

Intersect A with features from multiple B tables:

.. code-block:: sql

   WITH all_b_features AS (
       SELECT * FROM features_b1
       UNION ALL
       SELECT * FROM features_b2
       UNION ALL
       SELECT * FROM features_b3
   )
   SELECT DISTINCT a.*
   FROM features_a a
   INNER JOIN all_b_features b ON a.interval INTERSECTS b.interval

**Use case:** Find features overlapping any region from multiple annotation sources.

Track Overlap Source
~~~~~~~~~~~~~~~~~~~~

Know which source table each overlap came from:

.. code-block:: sql

   WITH all_b_features AS (
       SELECT *, 'source1' AS source FROM features_b1
       UNION ALL
       SELECT *, 'source2' AS source FROM features_b2
       UNION ALL
       SELECT *, 'source3' AS source FROM features_b3
   )
   SELECT a.*, b.name AS overlap_name, b.source
   FROM features_a a
   INNER JOIN all_b_features b ON a.interval INTERSECTS b.interval

**Use case:** Track which annotation database each overlap originated from.

Complex Filtering
-----------------

Overlap with Quality Filters
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Combine spatial and attribute filters:

.. code-block:: sql

   SELECT v.*, g.name AS gene_name
   FROM variants v
   INNER JOIN genes g ON v.interval INTERSECTS g.interval
   WHERE v.quality >= 30
     AND g.biotype = 'protein_coding'
   ORDER BY v.chrom, v.start

**Use case:** Find high-quality variants in protein-coding genes.

Specific Target Genes
~~~~~~~~~~~~~~~~~~~~~

Find overlaps with a specific set of genes:

.. code-block:: sql

   SELECT v.*, g.name AS gene_name
   FROM variants v
   INNER JOIN genes g ON v.interval INTERSECTS g.interval
   WHERE g.name IN ('BRCA1', 'BRCA2', 'TP53', 'EGFR')
   ORDER BY g.name, v.start

**Use case:** Extract variants in clinically relevant genes.
