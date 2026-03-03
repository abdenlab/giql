Transpilation
=============

The ``giql`` Python package transpiles GIQL into SQL.

How it works
------------

When you do this:

.. code-block:: python

   from giql import transpile

   sql = transpile(
       "SELECT * FROM variants WHERE interval INTERSECTS 'chr1:1000-2000'",
       tables=["variants"],
   )

   print(sql)

The transpiler performs three main steps:

1. **Parses** the GIQL query into an abstract syntax tree (AST) to identify GIQL-specific operators
2. **Transforms** genomic operators into SQL predicates and Common Table Expressions (CTEs), and replace genomic pseudo-columns with actual column references
3. **Generates** SQL output from the modified AST

The result is a standard SQL query that can be consumed by an execution engine that is not genome-aware.

.. code-block:: sql

   SELECT * FROM variants
   WHERE "chrom" = 'chr1' AND "start" < 2000 AND "end" > 1000

Notably, the transpiler expands logical genomic range columns into physical column comparisons.

The ``Table`` configuration of ``"variants"`` tells GIQL which physical columns correspond to 
the logical ``interval`` column. The above example simply maps to the default column names: 
``chrom``, ``start``, ``end``.


Examples
--------

Each GIQL operator expands to specific SQL patterns.

**INTERSECTS** expands to range overlap checks:

.. tab-set::

   .. tab-item:: GIQL

      .. code-block:: sql

         a.interval INTERSECTS b.interval

   .. tab-item:: SQL

      .. code-block:: sql

         a."chrom" = b."chrom"
         AND a."start" < b."end"
         AND a."end" > b."start"

**CONTAINS** expands to containment checks:

.. tab-set::

   .. tab-item:: GIQL

      .. code-block:: sql

         a.interval CONTAINS b.interval

   .. tab-item:: SQL

      .. code-block:: sql

         a."chrom" = b."chrom"
         AND a."start" <= b."start"
         AND a."end" >= b."end"

**DISTANCE** expands to gap calculations:

.. tab-set::

   .. tab-item:: GIQL

      .. code-block:: sql

         DISTANCE(a.interval, b.interval)

   .. tab-item:: SQL

      .. code-block:: sql

         CASE
             WHEN a."chrom" != b."chrom" THEN NULL
             WHEN a."end" <= b."start" THEN b."start" - a."end"
             WHEN b."end" <= a."start" THEN a."start" - b."end"
             ELSE 0
         END

**Intersection joins** expand to inequality joins:

.. tab-set::

   .. tab-item:: GIQL

      .. code-block:: sql

         SELECT v.*, g.name AS gene_name
         FROM variants v
         JOIN genes g ON v.interval INTERSECTS g.interval
         WHERE v.quality >= 30

   .. tab-item:: SQL

      .. code-block:: sql

         SELECT v.*, g.name AS gene_name
         FROM variants AS v
         JOIN genes AS g
           ON v."chrom" = g."chrom"
           AND v."start" < g."end"
           AND v."end" > g."start"
         WHERE v.quality >= 30

**NEAREST** expands to lateral subqueries:

.. tab-set::

   .. tab-item:: GIQL

      .. code-block:: sql

         SELECT peaks.name, nearest.name, nearest.distance
         FROM peaks
         CROSS JOIN LATERAL NEAREST(
             genes, reference=peaks.interval, k=5
         ) AS nearest

   .. tab-item:: SQL

      .. code-block:: sql

         SELECT peaks.name, nearest.name, nearest.distance
         FROM peaks
         CROSS JOIN LATERAL (
             SELECT
                 genes.*,
                 CASE
                     WHEN peaks."chrom" != genes."chrom" THEN NULL
                     WHEN peaks."start" < genes."end"
                      AND peaks."end" > genes."start" THEN 0
                     WHEN peaks."end" <= genes."start"
                       THEN genes."start" - peaks."end"
                     ELSE peaks."start" - genes."end"
                 END AS distance
             FROM genes
             WHERE peaks."chrom" = genes."chrom"
             ORDER BY ABS(
                 CASE
                     WHEN peaks."chrom" != genes."chrom" THEN NULL
                     WHEN peaks."start" < genes."end"
                      AND peaks."end" > genes."start" THEN 0
                     WHEN peaks."end" <= genes."start"
                       THEN genes."start" - peaks."end"
                     ELSE peaks."start" - genes."end"
                 END
             )
             LIMIT 5
         ) AS nearest

**MERGE** expands to window-function-based clustering:

.. tab-set::

   .. tab-item:: GIQL

      .. code-block:: sql

         SELECT MERGE(interval), COUNT(*) AS count
         FROM features

   .. tab-item:: SQL

      .. code-block:: sql

         SELECT
             "chrom",
             MIN("start") AS start,
             MAX("end") AS end,
             COUNT(*) AS count
         FROM (
             SELECT
                 *,
                 SUM(is_new_cluster) OVER (
                     PARTITION BY "chrom"
                     ORDER BY "start" NULLS LAST
                 ) AS __giql_cluster_id
             FROM (
                 SELECT
                     *,
                     CASE
                         WHEN LAG("end") OVER (
                             PARTITION BY "chrom"
                             ORDER BY "start" NULLS LAST
                         ) >= "start" THEN 0
                         ELSE 1
                     END AS is_new_cluster
                 FROM features
             ) AS lag_calc
         ) AS clustered
         GROUP BY chrom, __giql_cluster_id
         ORDER BY "chrom" NULLS LAST, "start" NULLS LAST
