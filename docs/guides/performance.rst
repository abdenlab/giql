Performance Guide
=================

This guide covers strategies for optimizing GIQL query performance, including
indexing, query patterns, and backend-specific optimizations.

.. contents::
   :local:
   :depth: 2

Understanding Query Performance
-------------------------------

How GIQL Queries Execute
~~~~~~~~~~~~~~~~~~~~~~~~

When you execute a GIQL query:

1. GIQL parses the query and identifies genomic operators
2. Operators are expanded into standard SQL predicates
3. The SQL is sent to the database backend
4. The database executes the query using its optimizer

Performance depends on both the generated SQL and how the database executes it.

Common Performance Bottlenecks
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

- **Full table scans**: No indexes to speed up filtering
- **Cartesian products**: Large cross joins without early filtering
- **Missing chromosome filters**: Comparing features across all chromosomes
- **Inefficient join order**: Small tables should drive joins

Indexing Strategies
-------------------

Creating Indexes
~~~~~~~~~~~~~~~~

Create indexes on genomic columns for faster queries:

.. code-block:: python

   # DuckDB
   engine.conn.execute("""
       CREATE INDEX idx_features_position
       ON features (chromosome, start_pos, end_pos)
   """)

   # SQLite
   engine.conn.execute("""
       CREATE INDEX idx_features_position
       ON features (chromosome, start_pos, end_pos)
   """)

Recommended Index Patterns
~~~~~~~~~~~~~~~~~~~~~~~~~~

**For single-table queries (filtering):**

.. code-block:: sql

   CREATE INDEX idx_table_position ON table_name (chromosome, start_pos, end_pos)

**For join queries:**

.. code-block:: sql

   -- Index both tables involved in joins
   CREATE INDEX idx_variants_position ON variants (chromosome, start_pos, end_pos)
   CREATE INDEX idx_genes_position ON genes (chromosome, start_pos, end_pos)

**For strand-specific queries:**

.. code-block:: sql

   CREATE INDEX idx_features_strand ON features (chromosome, strand, start_pos, end_pos)

When to Create Indexes
~~~~~~~~~~~~~~~~~~~~~~

Create indexes when:

- Tables have more than ~10,000 rows
- You're running repeated queries on the same tables
- Join queries are slow
- Filtering by genomic position is common

Skip indexes when:

- Tables are small
- You're doing one-time analysis
- Full table scans are acceptable

Query Optimization Patterns
---------------------------

Pre-filter by Chromosome
~~~~~~~~~~~~~~~~~~~~~~~~

Always include chromosome filtering when joining tables:

.. code-block:: python

   # Good: Explicit chromosome filter
   cursor = engine.execute("""
       SELECT a.*, b.name
       FROM features_a a
       JOIN features_b b ON a.position INTERSECTS b.position
       WHERE a.chromosome = 'chr1'
   """)

   # Also good: Cross-chromosome join with implicit filtering
   # GIQL handles this, but explicit is clearer
   cursor = engine.execute("""
       SELECT a.*, b.name
       FROM features_a a
       JOIN features_b b ON a.position INTERSECTS b.position
         AND a.chromosome = b.chromosome
   """)

Use Selective Filters Early
~~~~~~~~~~~~~~~~~~~~~~~~~~~

Apply selective filters before joins:

.. code-block:: python

   # Good: Filter before joining
   cursor = engine.execute("""
       WITH filtered_variants AS (
           SELECT * FROM variants
           WHERE quality >= 30 AND filter = 'PASS'
       )
       SELECT f.*, g.name
       FROM filtered_variants f
       JOIN genes g ON f.position INTERSECTS g.position
   """)

   # Less efficient: Filter after joining
   cursor = engine.execute("""
       SELECT v.*, g.name
       FROM variants v
       JOIN genes g ON v.position INTERSECTS g.position
       WHERE v.quality >= 30 AND v.filter = 'PASS'
   """)

Limit Result Sets
~~~~~~~~~~~~~~~~~

Use LIMIT for exploratory queries:

.. code-block:: python

   # Good: Limit results during exploration
   cursor = engine.execute("""
       SELECT * FROM variants
       WHERE position INTERSECTS 'chr1:1000000-2000000'
       LIMIT 100
   """)

Use DISTINCT Wisely
~~~~~~~~~~~~~~~~~~~

DISTINCT can be expensive. Only use when necessary:

.. code-block:: python

   # Only use DISTINCT when you actually need unique rows
   cursor = engine.execute("""
       SELECT DISTINCT a.*
       FROM features_a a
       JOIN features_b b ON a.position INTERSECTS b.position
   """)

   # If you just need to check existence, use EXISTS instead
   cursor = engine.execute("""
       SELECT a.*
       FROM features_a a
       WHERE EXISTS (
           SELECT 1 FROM features_b b
           WHERE a.position INTERSECTS b.position
       )
   """)

NEAREST Query Optimization
--------------------------

Optimizing K-NN Queries
~~~~~~~~~~~~~~~~~~~~~~~

The NEAREST operator can be expensive for large datasets. Optimize with:

**1. Use max_distance to limit search space:**

.. code-block:: python

   # Good: Constrained search
   cursor = engine.execute("""
       SELECT peaks.name, nearest.name, nearest.distance
       FROM peaks
       CROSS JOIN LATERAL NEAREST(
           genes,
           reference=peaks.position,
           k=5,
           max_distance=100000   -- Only search within 100kb
       ) AS nearest
   """)

**2. Request only the k you need:**

.. code-block:: python

   # Good: Request exactly what you need
   NEAREST(genes, reference=peaks.position, k=3)

   # Wasteful: Request more than needed
   NEAREST(genes, reference=peaks.position, k=100)

**3. Index the target table:**

.. code-block:: sql

   CREATE INDEX idx_genes_position ON genes (chromosome, start_pos, end_pos)

Merge and Cluster Optimization
------------------------------

Efficient Clustering
~~~~~~~~~~~~~~~~~~~~

For large datasets, consider pre-sorting:

.. code-block:: python

   # Pre-sort data for clustering
   cursor = engine.execute("""
       WITH sorted AS (
           SELECT * FROM features
           ORDER BY chromosome, start_pos
       )
       SELECT *, CLUSTER(position) AS cluster_id
       FROM sorted
   """)

Efficient Merging
~~~~~~~~~~~~~~~~~

Filter before merging to reduce data volume:

.. code-block:: python

   # Good: Filter first, then merge
   cursor = engine.execute("""
       WITH filtered AS (
           SELECT * FROM features
           WHERE score >= 10
       )
       SELECT MERGE(position), COUNT(*) AS count
       FROM filtered
   """)

Analyzing Query Performance
---------------------------

Using EXPLAIN
~~~~~~~~~~~~~

Analyze query execution plans:

.. code-block:: python

   # Get the transpiled SQL
   sql = engine.transpile("""
       SELECT a.*, b.name
       FROM variants a
       JOIN genes b ON a.position INTERSECTS b.position
   """)

   # Analyze the execution plan
   cursor = engine.execute(f"EXPLAIN {sql}")
   for row in cursor:
       print(row)

   # DuckDB also supports EXPLAIN ANALYZE for actual timing
   cursor = engine.execute(f"EXPLAIN ANALYZE {sql}")

Timing Queries
~~~~~~~~~~~~~~

Measure query execution time:

.. code-block:: python

   import time

   start = time.time()
   cursor = engine.execute("""
       SELECT * FROM variants
       WHERE position INTERSECTS 'chr1:1000000-2000000'
   """)
   results = cursor.fetchall()
   elapsed = time.time() - start

   print(f"Query returned {len(results)} rows in {elapsed:.2f} seconds")

Backend-Specific Tips
---------------------

DuckDB Optimizations
~~~~~~~~~~~~~~~~~~~~

**Use columnar strengths:**

DuckDB is columnar, so queries that select few columns are faster:

.. code-block:: python

   # Faster: Select only needed columns
   cursor = engine.execute("""
       SELECT chromosome, start_pos, end_pos, name
       FROM features
       WHERE position INTERSECTS 'chr1:1000-2000'
   """)

   # Slower: Select all columns
   cursor = engine.execute("""
       SELECT *
       FROM features
       WHERE position INTERSECTS 'chr1:1000-2000'
   """)

**Parallel execution:**

DuckDB automatically parallelizes queries. For very large datasets,
ensure you're not limiting parallelism.

SQLite Optimizations
~~~~~~~~~~~~~~~~~~~~

**Use covering indexes:**

.. code-block:: sql

   -- Include commonly selected columns in the index
   CREATE INDEX idx_features_covering
   ON features (chromosome, start_pos, end_pos, name, score)

**Analyze tables:**

.. code-block:: python

   # Help SQLite's query planner
   engine.conn.execute("ANALYZE features")

Memory Management
-----------------

Streaming Results
~~~~~~~~~~~~~~~~~

For large result sets, iterate instead of fetching all:

.. code-block:: python

   # Good: Stream results
   cursor = engine.execute("SELECT * FROM large_table")
   for row in cursor:
       process(row)

   # Memory-intensive: Fetch all at once
   cursor = engine.execute("SELECT * FROM large_table")
   all_rows = cursor.fetchall()  # Loads everything into memory

Batch Processing
~~~~~~~~~~~~~~~~

Process large datasets in batches:

.. code-block:: python

   chromosomes = ['chr1', 'chr2', 'chr3', ...]  # All chromosomes

   for chrom in chromosomes:
       cursor = engine.execute(f"""
           SELECT * FROM features
           WHERE chromosome = '{chrom}'
             AND position INTERSECTS '{chrom}:1-1000000'
       """)
       process_chromosome(cursor)

Performance Checklist
---------------------

Before running large queries, check:

.. code-block:: text

   □ Indexes created on genomic columns
   □ Chromosome filtering included in joins
   □ Selective filters applied early
   □ LIMIT used for exploration
   □ Only necessary columns selected
   □ NEAREST queries use max_distance
   □ Results streamed instead of fetched all at once

Quick Wins
~~~~~~~~~~

1. **Add indexes** - Usually the biggest performance improvement
2. **Filter by chromosome** - Reduces join complexity significantly
3. **Use max_distance with NEAREST** - Limits search space
4. **Stream results** - Reduces memory pressure
5. **Use DuckDB** - Generally faster for analytical queries
