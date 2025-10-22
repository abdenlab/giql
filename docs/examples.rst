Examples
========

This section provides examples of common GIQL usage patterns.

Replicating Bedtools with GIQL
-------------------------------

These examples show how to replicate common ``bedtools intersect`` operations using GIQL.

Setup
~~~~~

.. code-block:: python

   from giql import GIQLEngine

   # Load two genomic datasets
   engine = GIQLEngine(target_dialect="duckdb")
   engine.load_csv("features_a", "file_a.bed")
   engine.load_csv("features_b", "file_b.bed")

   # Register schemas
   for table in ["features_a", "features_b"]:
       engine.register_table_schema(
           table,
           {
               "chromosome": "VARCHAR",
               "start_pos": "BIGINT",
               "end_pos": "BIGINT",
               "name": "VARCHAR",
               "score": "FLOAT",
               "strand": "VARCHAR",
           },
           genomic_column="position",
       )

Basic Intersect Operations
~~~~~~~~~~~~~~~~~~~~~~~~~~~

Default: Report overlaps between A and B
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Bedtools:**

.. code-block:: bash

   bedtools intersect -a file_a.bed -b file_b.bed

**GIQL:**

.. code-block:: python

   cursor = engine.execute("""
       SELECT DISTINCT a.*
       FROM features_a a, features_b b
       WHERE a.position INTERSECTS b.position
   """)

``-wa``: Write original A entry for each overlap
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Bedtools:**

.. code-block:: bash

   bedtools intersect -a file_a.bed -b file_b.bed -wa

**GIQL:**

.. code-block:: python

   cursor = engine.execute("""
       SELECT a.*
       FROM features_a a, features_b b
       WHERE a.position INTERSECTS b.position
   """)

``-wb``: Write original B entry for each overlap
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Bedtools:**

.. code-block:: bash

   bedtools intersect -a file_a.bed -b file_b.bed -wb

**GIQL:**

.. code-block:: python

   cursor = engine.execute("""
       SELECT b.*
       FROM features_a a, features_b b
       WHERE a.position INTERSECTS b.position
   """)

``-wa -wb``: Write both A and B entries
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Bedtools:**

.. code-block:: bash

   bedtools intersect -a file_a.bed -b file_b.bed -wa -wb

**GIQL:**

.. code-block:: python

   cursor = engine.execute("""
       SELECT a.*, b.*
       FROM features_a a, features_b b
       WHERE a.position INTERSECTS b.position
   """)

Filtering Operations
~~~~~~~~~~~~~~~~~~~~

``-v``: Report A entries with NO overlap in B
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Bedtools:**

.. code-block:: bash

   bedtools intersect -a file_a.bed -b file_b.bed -v

**GIQL:**

.. code-block:: python

   cursor = engine.execute("""
       SELECT a.*
       FROM features_a a
       LEFT JOIN features_b b ON a.position INTERSECTS b.position
       WHERE b.chromosome IS NULL
   """)

``-u``: Report A entries with ANY overlap in B (unique)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Bedtools:**

.. code-block:: bash

   bedtools intersect -a file_a.bed -b file_b.bed -u

**GIQL:**

.. code-block:: python

   cursor = engine.execute("""
       SELECT DISTINCT a.*
       FROM features_a a
       INNER JOIN features_b b ON a.position INTERSECTS b.position
   """)

Counting Operations
~~~~~~~~~~~~~~~~~~~

``-c``: Count B overlaps for each A feature
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Bedtools:**

.. code-block:: bash

   bedtools intersect -a file_a.bed -b file_b.bed -c

**GIQL:**

.. code-block:: python

   cursor = engine.execute("""
       SELECT a.*, COUNT(b.name) as overlap_count
       FROM features_a a
       LEFT JOIN features_b b ON a.position INTERSECTS b.position
       GROUP BY a.chromosome, a.start_pos, a.end_pos, a.name, a.score, a.strand
   """)

Strand-Specific Operations
~~~~~~~~~~~~~~~~~~~~~~~~~~~

``-s``: Same strand overlaps only
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Bedtools:**

.. code-block:: bash

   bedtools intersect -a file_a.bed -b file_b.bed -s

**GIQL:**

.. code-block:: python

   cursor = engine.execute("""
       SELECT a.*
       FROM features_a a, features_b b
       WHERE a.position INTERSECTS b.position
         AND a.strand = b.strand
   """)

``-S``: Opposite strand overlaps only
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Bedtools:**

.. code-block:: bash

   bedtools intersect -a file_a.bed -b file_b.bed -S

**GIQL:**

.. code-block:: python

   cursor = engine.execute("""
       SELECT a.*
       FROM features_a a, features_b b
       WHERE a.position INTERSECTS b.position
         AND a.strand != b.strand
         AND a.strand IN ('+', '-')
         AND b.strand IN ('+', '-')
   """)

Overlap Fraction Requirements
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

``-f``: Minimum overlap fraction of A
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Bedtools:**

.. code-block:: bash

   bedtools intersect -a file_a.bed -b file_b.bed -f 0.5

**GIQL:**

.. code-block:: python

   cursor = engine.execute("""
       SELECT a.*
       FROM features_a a, features_b b
       WHERE a.position INTERSECTS b.position
         AND (
             LEAST(a.end_pos, b.end_pos) - GREATEST(a.start_pos, b.start_pos)
         ) >= 0.5 * (a.end_pos - a.start_pos)
   """)

``-F``: Minimum overlap fraction of B
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Bedtools:**

.. code-block:: bash

   bedtools intersect -a file_a.bed -b file_b.bed -F 0.5

**GIQL:**

.. code-block:: python

   cursor = engine.execute("""
       SELECT a.*
       FROM features_a a, features_b b
       WHERE a.position INTERSECTS b.position
         AND (
             LEAST(a.end_pos, b.end_pos) - GREATEST(a.start_pos, b.start_pos)
         ) >= 0.5 * (b.end_pos - b.start_pos)
   """)

``-r``: Reciprocal overlap (both A and B must meet fraction)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Bedtools:**

.. code-block:: bash

   bedtools intersect -a file_a.bed -b file_b.bed -f 0.5 -r

**GIQL:**

.. code-block:: python

   cursor = engine.execute("""
       WITH overlap_calcs AS (
           SELECT
               a.*,
               (LEAST(a.end_pos, b.end_pos) - GREATEST(a.start_pos, b.start_pos)) as overlap_bp,
               (a.end_pos - a.start_pos) as a_length,
               (b.end_pos - b.start_pos) as b_length
           FROM features_a a, features_b b
           WHERE a.position INTERSECTS b.position
       )
       SELECT chromosome, start_pos, end_pos, name, score, strand
       FROM overlap_calcs
       WHERE overlap_bp >= 0.5 * a_length
         AND overlap_bp >= 0.5 * b_length
   """)

Join Operations
~~~~~~~~~~~~~~~

``-loj``: Left outer join (report all A, with NULL for non-overlapping)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Bedtools:**

.. code-block:: bash

   bedtools intersect -a file_a.bed -b file_b.bed -loj

**GIQL:**

.. code-block:: python

   cursor = engine.execute("""
       SELECT a.*, b.*
       FROM features_a a
       LEFT JOIN features_b b ON a.position INTERSECTS b.position
   """)

``-wo``: Write overlap amount in base pairs
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Bedtools:**

.. code-block:: bash

   bedtools intersect -a file_a.bed -b file_b.bed -wo

**GIQL:**

.. code-block:: python

   cursor = engine.execute("""
       SELECT
           a.*,
           b.*,
           (LEAST(a.end_pos, b.end_pos) - GREATEST(a.start_pos, b.start_pos)) as overlap_bp
       FROM features_a a, features_b b
       WHERE a.position INTERSECTS b.position
   """)

``-wao``: Write overlap amount for ALL A features
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Bedtools:**

.. code-block:: bash

   bedtools intersect -a file_a.bed -b file_b.bed -wao

**GIQL:**

.. code-block:: python

   cursor = engine.execute("""
       SELECT
           a.*,
           b.*,
           CASE
               WHEN b.chromosome IS NULL THEN 0
               ELSE LEAST(a.end_pos, b.end_pos) - GREATEST(a.start_pos, b.start_pos)
           END as overlap_bp
       FROM features_a a
       LEFT JOIN features_b b ON a.position INTERSECTS b.position
   """)

Multiple Database Files
~~~~~~~~~~~~~~~~~~~~~~~~

Intersect with multiple B files
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Bedtools:**

.. code-block:: bash

   bedtools intersect -a file_a.bed -b file1.bed file2.bed file3.bed

**GIQL:**

.. code-block:: python

   # Load multiple files
   engine.load_csv("features_b1", "file1.bed")
   engine.load_csv("features_b2", "file2.bed")
   engine.load_csv("features_b3", "file3.bed")

   # Register schemas for each...

   # Query using UNION to combine all B features
   cursor = engine.execute("""
       WITH all_b_features AS (
           SELECT * FROM features_b1
           UNION ALL
           SELECT * FROM features_b2
           UNION ALL
           SELECT * FROM features_b3
       )
       SELECT DISTINCT a.*
       FROM features_a a
       INNER JOIN all_b_features b ON a.position INTERSECTS b.position
   """)

Bedtools Advanced Queries
~~~~~~~~~~~~~~~~~~~~~~~~~~

Find features in A that overlap at least 2 features in B
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Bedtools:**

.. code-block:: bash

   bedtools intersect -a file_a.bed -b file_b.bed -c | awk '$NF >= 2'

**GIQL:**

.. code-block:: python

   cursor = engine.execute("""
       SELECT a.*
       FROM features_a a
       INNER JOIN features_b b ON a.position INTERSECTS b.position
       GROUP BY a.chromosome, a.start_pos, a.end_pos, a.name, a.score, a.strand
       HAVING COUNT(*) >= 2
   """)

Complex filtering: Overlaps specific genes with quality threshold
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: python

   cursor = engine.execute("""
       SELECT v.*, g.name as gene_name
       FROM variants v
       INNER JOIN genes g ON v.position INTERSECTS g.position
       WHERE v.quality >= 30
         AND g.name IN ('BRCA1', 'BRCA2', 'TP53')
         AND v.chromosome = g.chromosome
       ORDER BY v.chromosome, v.start_pos
   """)

Aggregate statistics per chromosome
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: python

   cursor = engine.execute("""
       SELECT
           a.chromosome,
           COUNT(DISTINCT a.name) as total_features,
           COUNT(b.name) as total_overlaps,
           AVG(LEAST(a.end_pos, b.end_pos) - GREATEST(a.start_pos, b.start_pos)) as avg_overlap_bp
       FROM features_a a
       LEFT JOIN features_b b ON a.position INTERSECTS b.position
       GROUP BY a.chromosome
       ORDER BY a.chromosome
   """)

Advanced Queries
----------------

Multiple Range Matching
~~~~~~~~~~~~~~~~~~~~~~~~

Find variants in multiple genomic regions:

.. code-block:: python

   cursor = engine.execute("""
       SELECT * FROM variants
       WHERE position INTERSECTS ANY(
           'chr1:1000-2000',
           'chr1:5000-6000',
           'chr2:1000-3000'
       )
   """)

Complex Filtering
~~~~~~~~~~~~~~~~~

Find high-quality variants in specific genes:

.. code-block:: python

   cursor = engine.execute("""
       SELECT v.*, g.name as gene_name
       FROM variants v
       INNER JOIN genes g ON v.position INTERSECTS g.position
       WHERE v.quality >= 30
         AND g.name IN ('BRCA1', 'BRCA2', 'TP53')
       ORDER BY v.chromosome, v.start_pos
   """)

Aggregate Statistics
~~~~~~~~~~~~~~~~~~~~

Calculate per-chromosome statistics:

.. code-block:: python

   cursor = engine.execute("""
       SELECT
           a.chromosome,
           COUNT(DISTINCT a.name) as total_features,
           COUNT(b.name) as total_overlaps,
           AVG(LEAST(a.end_pos, b.end_pos) - GREATEST(a.start_pos, b.start_pos)) as avg_overlap_bp
       FROM features_a a
       LEFT JOIN features_b b ON a.position INTERSECTS b.position
       GROUP BY a.chromosome
       ORDER BY a.chromosome
   """)

Common Patterns
---------------

Find Non-Overlapping Features
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   cursor = engine.execute("""
       SELECT a.*
       FROM features_a a
       LEFT JOIN features_b b ON a.position INTERSECTS b.position
       WHERE b.chromosome IS NULL
   """)

Nearest Neighbors
~~~~~~~~~~~~~~~~~

Find features within a distance window:

.. code-block:: python

   cursor = engine.execute("""
       WITH expanded AS (
           SELECT
               id,
               chromosome,
               start_pos - 1000 as search_start,
               end_pos + 1000 as search_end
           FROM features_a
       )
       SELECT a.*, b.*
       FROM expanded a
       JOIN features_b b
           ON b.chromosome = a.chromosome
           AND b.start_pos < a.search_end
           AND b.end_pos > a.search_start
   """)

Debugging and Transpilation
----------------------------

Understanding Generated SQL
~~~~~~~~~~~~~~~~~~~~~~~~~~~

Use the ``transpile()`` method to see what SQL is generated for your GIQL queries:

.. code-block:: python

   from giql import GIQLEngine

   with GIQLEngine(target_dialect="duckdb") as engine:
       # Register schema
       engine.register_table_schema(
           "variants",
           {
               "chromosome": "VARCHAR",
               "start_pos": "BIGINT",
               "end_pos": "BIGINT",
           },
           genomic_column="position",
       )

       # Transpile without executing
       sql = engine.transpile("""
           SELECT * FROM variants
           WHERE position INTERSECTS 'chr1:1000-2000'
       """)

       print(sql)
       # Shows: SELECT * FROM variants WHERE chromosome = 'chr1' AND start_pos < 2000 AND end_pos > 1000

Comparing Dialects
~~~~~~~~~~~~~~~~~~

See how different dialects handle the same query:

.. code-block:: python

   from giql import GIQLEngine

   query = """
       SELECT * FROM variants
       WHERE position INTERSECTS ANY('chr1:1000-2000', 'chr2:5000-6000')
   """

   # DuckDB
   with GIQLEngine(target_dialect="duckdb") as engine:
       engine.register_table_schema("variants", {
           "chromosome": "VARCHAR",
           "start_pos": "BIGINT",
           "end_pos": "BIGINT",
       }, genomic_column="position")
       duckdb_sql = engine.transpile(query)
       print("DuckDB SQL:")
       print(duckdb_sql)

   # SQLite
   with GIQLEngine(target_dialect="sqlite") as engine:
       engine.register_table_schema("variants", {
           "chromosome": "VARCHAR",
           "start_pos": "BIGINT",
           "end_pos": "BIGINT",
       }, genomic_column="position")
       sqlite_sql = engine.transpile(query)
       print("\nSQLite SQL:")
       print(sqlite_sql)

Verbose Mode for Debugging
~~~~~~~~~~~~~~~~~~~~~~~~~~~

Enable verbose mode to see detailed transpilation information:

.. code-block:: python

   from giql import GIQLEngine

   with GIQLEngine(target_dialect="duckdb", verbose=True) as engine:
       engine.register_table_schema("variants", {
           "chromosome": "VARCHAR",
           "start_pos": "BIGINT",
           "end_pos": "BIGINT",
       }, genomic_column="position")

       # This will print detailed transpilation steps
       sql = engine.transpile("""
           SELECT * FROM variants
           WHERE position INTERSECTS 'chr1:1000-2000'
       """)

       # Also shows transpilation details when executing
       cursor = engine.execute("""
           SELECT * FROM variants
           WHERE position INTERSECTS 'chr1:1000-2000'
       """)

Using Transpiled SQL Externally
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The transpiled SQL can be used with external tools or libraries:

.. code-block:: python

   from giql import GIQLEngine
   import duckdb

   # Get transpiled SQL
   with GIQLEngine(target_dialect="duckdb") as engine:
       engine.register_table_schema("variants", {
           "chromosome": "VARCHAR",
           "start_pos": "BIGINT",
           "end_pos": "BIGINT",
       }, genomic_column="position")

       sql = engine.transpile("""
           SELECT * FROM variants
           WHERE position INTERSECTS 'chr1:1000-2000'
       """)

   # Execute with external DuckDB connection
   conn = duckdb.connect("my_database.duckdb")
   result = conn.execute(sql).fetchall()
   conn.close()
