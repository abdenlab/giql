Performance
-----------

When you use GIQL:

1. GIQL parses the query and identifies genomic operators
2. Operators are expanded into SQL predicates
3. You execute the SQL on your database backend or analytics engine
4. The system optimizes the query and executes it

Performance depends on both the generated SQL and how the target data system
plans, optimizes, and executes it. Some common performance bottlenecks include:

- **Full table scans**: No indexes to speed up filtering
- **Cartesian products**: Large cross joins without early filtering
- **Missing chromosome filters**: Comparing features across all chromosomes
- **Inefficient join order**: Small tables should drive joins

Streaming
---------

Analytics engines like DuckDB and Polars support streaming data sources
in sequences of small "record batches", enabling parallel processing and 
out-of-core workflows on files that may be much larger than memory.

For delimited text files, you can use native APIs:

**DuckDB:**

.. code-block:: python

   import duckdb
   from giql import transpile

   conn = duckdb.connect()

   # DuckDB can query CSV/TSV files directly
   conn.execute("""
       CREATE VIEW peaks AS
       SELECT * FROM read_csv('peaks.bed', delim='\t',
           columns={'chrom': 'VARCHAR', 'start': 'INTEGER',
                    'end': 'INTEGER', 'name': 'VARCHAR',
                    'score': 'INTEGER', 'strand': 'VARCHAR'})
   """)

   sql = transpile(
       "SELECT * FROM peaks WHERE interval INTERSECTS 'chr1:1000-2000'",
       tables=["peaks"],
   )

   df = conn.execute(sql).fetchdf()

**Polars:**

.. code-block:: python

   import polars as pl
   from giql import transpile

   lf = pl.scan_csv("peaks.bed", separator="\t",
       new_columns=["chrom", "start", "end", "name", "score", "strand"])

   sql = transpile(
       "SELECT * FROM peaks WHERE interval INTERSECTS 'chr1:1000-2000'",
       tables=["peaks"],
   )

   ctx = pl.SQLContext(peaks=lf)
   df = ctx.execute(sql).collect()

For specialized NGS formats, you can supply streaming data using the
`oxbow <https://github.com/abdenlab/oxbow>`_ package:

.. code-block:: python

   import duckdb
   import oxbow as ox
   from giql import transpile

   conn = duckdb.connect()

   # Load a streaming data source as a DuckDB relation
   peaks = ox.read_bed("peaks.bed").to_duckdb(conn, "peaks")

   sql = transpile(
       "SELECT * FROM peaks WHERE interval INTERSECTS 'chr1:1000-2000'",
       tables=["peaks"],
   )

   df = conn.execute(sql).fetchdf()

.. code-block:: python

   import polars as pl
   import oxbow as ox
   from giql import transpile

   # Load a streaming data source as a Polars LazyFrame
   lf = ox.read_bed("peaks.bed").pl(lazy=True)

   sql = transpile(
      """SELECT *, CLUSTER(interval) AS cluster_id
         FROM features
         ORDER BY chrom, start
      """,
      tables=["peaks"],
   )
   ctx = pl.SQLContext(peaks=lf)
   ctx.execute(sql).sink_parquet("peaks_clustered.parquet")

Indexing
--------

If your data source is a database table, you can create indexes on 
genomic columns for faster queries:

.. code-block:: sql

   -- DuckDB or SQLite
   CREATE INDEX idx_features_position
   ON features (chrom, start, "end")

**For single-table queries (filtering):**

.. code-block:: sql

   CREATE INDEX idx_table_position ON table_name (chrom, start, "end")

**For join queries:**

.. code-block:: sql

   -- Index both tables involved in joins
   CREATE INDEX idx_variants_position ON variants (chrom, start, "end")
   CREATE INDEX idx_genes_position ON genes (chrom, start, "end")

**For strand-specific queries:**

.. code-block:: sql

   CREATE INDEX idx_features_strand ON features (chrom, strand, start, "end")

Create indexes when:

- Tables are very large
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

.. code-block:: sql

   -- Good: Explicit chromosome filter
   SELECT a.*, b.name
   FROM features_a a
   JOIN features_b b ON a.interval INTERSECTS b.interval
   WHERE a.chrom = 'chr1'

Use Selective Filters Early
~~~~~~~~~~~~~~~~~~~~~~~~~~~

Apply selective filters before joins:

.. code-block:: sql

   -- Good: Filter before joining
   WITH filtered_variants AS (
       SELECT * FROM variants
       WHERE quality >= 30 AND filter = 'PASS'
   )
   SELECT f.*, g.name
   FROM filtered_variants f
   JOIN genes g ON f.interval INTERSECTS g.interval

Limit Result Sets
~~~~~~~~~~~~~~~~~

Use LIMIT for exploratory queries:

.. code-block:: sql

   SELECT * FROM variants
   WHERE interval INTERSECTS 'chr1:1000000-2000000'
   LIMIT 100

Use DISTINCT Wisely
~~~~~~~~~~~~~~~~~~~

DISTINCT can be expensive. Only use when necessary:

.. code-block:: sql

   -- If you just need to check existence, use EXISTS instead
   SELECT a.*
   FROM features_a a
   WHERE EXISTS (
       SELECT 1 FROM features_b b
       WHERE a.interval INTERSECTS b.interval
   )

Optimizing K-NN Queries
~~~~~~~~~~~~~~~~~~~~~~~

The NEAREST operator can be expensive for large datasets. Optimize with:

**1. Use max_distance to limit search space:**

.. code-block:: sql

   SELECT peaks.name, nearest.name, nearest.distance
   FROM peaks
   CROSS JOIN LATERAL NEAREST(
       genes,
       reference=peaks.interval,
       k=5,
       max_distance=100000   -- Only search within 100kb
   ) AS nearest

**2. Request only the k you need:**

.. code-block:: sql

   -- Good: Request exactly what you need
   NEAREST(genes, reference=peaks.interval, k=3)

   -- Wasteful: Request more than needed
   NEAREST(genes, reference=peaks.interval, k=100)

**3. Index the target table:**

.. code-block:: sql

   CREATE INDEX idx_genes_position ON genes (chrom, start, "end")

Efficient Clustering
~~~~~~~~~~~~~~~~~~~~

For large datasets, consider pre-sorting:

.. code-block:: sql

   WITH sorted AS (
       SELECT * FROM features
       ORDER BY chrom, start
   )
   SELECT *, CLUSTER(interval) AS cluster_id
   FROM sorted

Efficient Merging
~~~~~~~~~~~~~~~~~

Filter before merging to reduce data volume:

.. code-block:: sql

   WITH filtered AS (
       SELECT * FROM features
       WHERE score >= 10
   )
   SELECT MERGE(interval), COUNT(*) AS count
   FROM filtered

Analyzing Query Performance
---------------------------

Using EXPLAIN
~~~~~~~~~~~~~

Analyze query execution plans by running EXPLAIN on the transpiled SQL:

.. code-block:: python

   from giql import transpile

   sql = transpile(
       """
       SELECT a.*, b.name
       FROM variants a
       JOIN genes b ON a.interval INTERSECTS b.interval
       """,
       tables=["variants", "genes"],
   )

   # Run EXPLAIN on your database connection
   # conn.execute(f"EXPLAIN {sql}")
   # DuckDB also supports EXPLAIN ANALYZE for actual timing

Backend-Specific Tips
---------------------

DuckDB Optimizations
~~~~~~~~~~~~~~~~~~~~

**Use columnar strengths:**

DuckDB is columnar, so queries that select few columns are faster:

.. code-block:: sql

   -- Faster: Select only needed columns
   SELECT chrom, start, "end", name
   FROM features
   WHERE interval INTERSECTS 'chr1:1000-2000'

**Parallel execution:**

DuckDB automatically parallelizes queries. For very large datasets,
ensure you're not limiting parallelism.

SQLite Optimizations
~~~~~~~~~~~~~~~~~~~~

**Use covering indexes:**

.. code-block:: sql

   -- Include commonly selected columns in the index
   CREATE INDEX idx_features_covering
   ON features (chrom, start, "end", name, score)

**Analyze tables:**

.. code-block:: sql

   ANALYZE features
