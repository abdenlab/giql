Advanced Queries
================

This section covers advanced query patterns including multi-range matching,
complex filtering, aggregate statistics, and multi-table workflows.

.. contents::
   :local:
   :depth: 2

Multi-Range Matching
--------------------

Match Any of Multiple Regions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Find features overlapping any of several regions of interest:

.. code-block:: python

   cursor = engine.execute("""
       SELECT * FROM variants
       WHERE interval INTERSECTS ANY(
           'chr1:1000000-2000000',
           'chr1:5000000-6000000',
           'chr2:1000000-3000000'
       )
   """)

**Use case:** Query multiple regions of interest in a single statement.

Match All of Multiple Points
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Find features containing all specified positions:

.. code-block:: python

   cursor = engine.execute("""
       SELECT * FROM genes
       WHERE interval CONTAINS ALL(
           'chr1:1500',
           'chr1:1600',
           'chr1:1700'
       )
   """)

**Use case:** Find genes spanning a set of SNP positions.

Exclude Multiple Regions
~~~~~~~~~~~~~~~~~~~~~~~~

Find features that don't overlap any blacklisted region:

.. code-block:: python

   cursor = engine.execute("""
       SELECT * FROM peaks
       WHERE NOT interval INTERSECTS ANY(
           'chr1:120000000-125000000',   -- Centromere region
           'chr1:140000000-142000000',   -- Known artifact
           'chrM:1-16569'                -- Mitochondrial
       )
   """)

**Use case:** Filter out features in problematic genomic regions.

Combine ANY and ALL
~~~~~~~~~~~~~~~~~~~

Complex multi-range logic:

.. code-block:: python

   cursor = engine.execute("""
       SELECT * FROM features
       WHERE interval INTERSECTS ANY('chr1:1000-2000', 'chr1:5000-6000')
         AND interval CONTAINS ALL('chr1:1100', 'chr1:1200')
   """)

**Use case:** Find features matching complex spatial criteria.

Complex Filtering
-----------------

Multi-Attribute Filtering
~~~~~~~~~~~~~~~~~~~~~~~~~

Combine spatial and attribute filters:

.. code-block:: python

   cursor = engine.execute("""
       SELECT v.*, g.name AS gene_name, g.biotype
       FROM variants v
       INNER JOIN genes g ON v.interval INTERSECTS g.interval
       WHERE v.quality >= 30
         AND v.filter = 'PASS'
         AND v.allele_frequency > 0.01
         AND g.biotype = 'protein_coding'
       ORDER BY v.chromosome, v.start_pos
   """)

**Use case:** Extract high-quality variants in protein-coding genes.

Target Gene Lists
~~~~~~~~~~~~~~~~~

Filter to specific genes of interest:

.. code-block:: python

   cursor = engine.execute("""
       SELECT v.*, g.name AS gene_name
       FROM variants v
       INNER JOIN genes g ON v.interval INTERSECTS g.interval
       WHERE g.name IN (
           'BRCA1', 'BRCA2', 'TP53', 'EGFR', 'KRAS',
           'BRAF', 'PIK3CA', 'PTEN', 'APC', 'ATM'
       )
       ORDER BY g.name, v.start_pos
   """)

**Use case:** Extract variants in clinically actionable genes.

Conditional Logic
~~~~~~~~~~~~~~~~~

Apply different criteria based on feature type:

.. code-block:: python

   cursor = engine.execute("""
       SELECT v.*, g.name, g.biotype,
           CASE
               WHEN g.biotype = 'protein_coding' THEN 'coding'
               WHEN g.biotype LIKE '%RNA%' THEN 'noncoding_RNA'
               ELSE 'other'
           END AS gene_category
       FROM variants v
       INNER JOIN genes g ON v.interval INTERSECTS g.interval
       WHERE CASE
           WHEN g.biotype = 'protein_coding' THEN v.quality >= 30
           ELSE v.quality >= 20
       END
   """)

**Use case:** Apply different quality thresholds based on genomic context.

Aggregate Statistics
--------------------

Per-Chromosome Statistics
~~~~~~~~~~~~~~~~~~~~~~~~~

Calculate summary statistics by chromosome:

.. code-block:: python

   cursor = engine.execute("""
       SELECT
           a.chromosome,
           COUNT(DISTINCT a.name) AS total_features,
           COUNT(b.name) AS total_overlaps,
           COUNT(DISTINCT CASE WHEN b.name IS NOT NULL THEN a.name END) AS features_with_overlap
       FROM features_a a
       LEFT JOIN features_b b ON a.interval INTERSECTS b.interval
       GROUP BY a.chromosome
       ORDER BY a.chromosome
   """)

**Use case:** Compare feature distribution across chromosomes.

Overlap Statistics
~~~~~~~~~~~~~~~~~~

Calculate overlap metrics:

.. code-block:: python

   cursor = engine.execute("""
       SELECT
           a.chromosome,
           COUNT(*) AS overlap_count,
           AVG(LEAST(a.end_pos, b.end_pos) - GREATEST(a.start_pos, b.start_pos)) AS avg_overlap_bp,
           SUM(LEAST(a.end_pos, b.end_pos) - GREATEST(a.start_pos, b.start_pos)) AS total_overlap_bp
       FROM features_a a
       INNER JOIN features_b b ON a.interval INTERSECTS b.interval
       GROUP BY a.chromosome
       ORDER BY a.chromosome
   """)

**Use case:** Quantify overlap patterns across the genome.

Feature Size Distribution
~~~~~~~~~~~~~~~~~~~~~~~~~

Analyze feature sizes by category:

.. code-block:: python

   cursor = engine.execute("""
       SELECT
           biotype,
           COUNT(*) AS count,
           AVG(end_pos - start_pos) AS avg_length,
           MIN(end_pos - start_pos) AS min_length,
           MAX(end_pos - start_pos) AS max_length
       FROM genes
       GROUP BY biotype
       ORDER BY count DESC
   """)

**Use case:** Compare size distributions across feature types.

Multi-Table Workflows
---------------------

Three-Way Intersection
~~~~~~~~~~~~~~~~~~~~~~

Find features overlapping in all three tables:

.. code-block:: python

   cursor = engine.execute("""
       SELECT DISTINCT a.*
       FROM features_a a
       INNER JOIN features_b b ON a.interval INTERSECTS b.interval
       INNER JOIN features_c c ON a.interval INTERSECTS c.interval
   """)

**Use case:** Find consensus regions across multiple datasets.

Hierarchical Annotations
~~~~~~~~~~~~~~~~~~~~~~~~

Join multiple annotation levels:

.. code-block:: python

   cursor = engine.execute("""
       SELECT
           v.name AS variant,
           e.name AS exon,
           t.name AS transcript,
           g.name AS gene
       FROM variants v
       INNER JOIN exons e ON v.interval INTERSECTS e.interval
       INNER JOIN transcripts t ON e.interval WITHIN t.interval
       INNER JOIN genes g ON t.interval WITHIN g.interval
   """)

**Use case:** Build hierarchical annotations for variants.

Union with Deduplication
~~~~~~~~~~~~~~~~~~~~~~~~

Combine features from multiple sources:

.. code-block:: python

   cursor = engine.execute("""
       WITH all_peaks AS (
           SELECT *, 'chip_seq' AS source FROM chip_peaks
           UNION ALL
           SELECT *, 'atac_seq' AS source FROM atac_peaks
           UNION ALL
           SELECT *, 'dnase_seq' AS source FROM dnase_peaks
       )
       SELECT
           chromosome,
           start_pos,
           end_pos,
           STRING_AGG(DISTINCT source, ',') AS sources,
           COUNT(DISTINCT source) AS source_count
       FROM all_peaks
       GROUP BY chromosome, start_pos, end_pos
       HAVING COUNT(DISTINCT source) >= 2
   """)

**Use case:** Find regulatory regions supported by multiple assays.

Subqueries and CTEs
-------------------

Filtered Subquery
~~~~~~~~~~~~~~~~~

Use subqueries to pre-filter data:

.. code-block:: python

   cursor = engine.execute("""
       SELECT v.*
       FROM variants v
       WHERE v.interval INTERSECTS ANY(
           SELECT position FROM genes WHERE biotype = 'protein_coding'
       )
   """)

**Use case:** Intersect with dynamically filtered reference data.

.. note::

   Subquery support depends on the target database backend.

Chained CTEs
~~~~~~~~~~~~

Build complex analyses with Common Table Expressions:

.. code-block:: python

   cursor = engine.execute("""
       WITH
       -- Step 1: Find high-quality variants
       hq_variants AS (
           SELECT * FROM variants
           WHERE quality >= 30 AND filter = 'PASS'
       ),
       -- Step 2: Annotate with genes
       annotated AS (
           SELECT v.*, g.name AS gene_name, g.biotype
           FROM hq_variants v
           LEFT JOIN genes g ON v.interval INTERSECTS g.interval
       ),
       -- Step 3: Summarize by gene
       gene_summary AS (
           SELECT
               gene_name,
               biotype,
               COUNT(*) AS variant_count
           FROM annotated
           WHERE gene_name IS NOT NULL
           GROUP BY gene_name, biotype
       )
       SELECT * FROM gene_summary
       ORDER BY variant_count DESC
       LIMIT 20
   """)

**Use case:** Build multi-step analysis pipelines in a single query.

Window Functions
----------------

Rank Overlaps
~~~~~~~~~~~~~

Rank features by their overlap characteristics:

.. code-block:: python

   cursor = engine.execute("""
       SELECT
           a.name,
           a.chromosome,
           a.start_pos,
           overlap_count,
           RANK() OVER (ORDER BY overlap_count DESC) AS rank
       FROM (
           SELECT a.*, COUNT(b.name) AS overlap_count
           FROM features_a a
           LEFT JOIN features_b b ON a.interval INTERSECTS b.interval
           GROUP BY a.chromosome, a.start_pos, a.end_pos, a.name, a.score, a.strand
       ) a
   """)

**Use case:** Identify features with the most overlaps.

Running Totals
~~~~~~~~~~~~~~

Calculate cumulative coverage:

.. code-block:: python

   cursor = engine.execute("""
       SELECT
           chromosome,
           start_pos,
           end_pos,
           end_pos - start_pos AS length,
           SUM(end_pos - start_pos) OVER (
               PARTITION BY chromosome
               ORDER BY start_pos
           ) AS cumulative_bp
       FROM features
       ORDER BY chromosome, start_pos
   """)

**Use case:** Track cumulative coverage along each chromosome.

Debugging and Optimization
--------------------------

View Generated SQL
~~~~~~~~~~~~~~~~~~

Use transpile() to see the SQL GIQL generates:

.. code-block:: python

   sql = engine.transpile("""
       SELECT * FROM variants
       WHERE interval INTERSECTS 'chr1:1000-2000'
   """)
   print(sql)
   # See the actual SQL that will be executed

**Use case:** Debug queries or understand GIQL's translation.

Verbose Mode
~~~~~~~~~~~~

Enable detailed logging:

.. code-block:: python

   with GIQLEngine(target_dialect="duckdb", verbose=True) as engine:
       # All queries will print transpilation details
       cursor = engine.execute("""
           SELECT * FROM variants
           WHERE interval INTERSECTS 'chr1:1000-2000'
       """)

**Use case:** Diagnose query translation issues.

Explain Query Plan
~~~~~~~~~~~~~~~~~~

Analyze query execution:

.. code-block:: python

   # First transpile to get the SQL
   sql = engine.transpile("""
       SELECT v.*, g.name
       FROM variants v
       JOIN genes g ON v.interval INTERSECTS g.interval
   """)

   # Then use database-native EXPLAIN
   cursor = engine.execute(f"EXPLAIN {sql}")
   for row in cursor:
       print(row)

**Use case:** Optimize slow queries by examining execution plans.
