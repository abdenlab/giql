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
       reference := peaks.interval,
       k := 5,
       max_distance := 100000   -- Only search within 100kb
   ) AS nearest

**2. Request only the k you need:**

.. code-block:: sql

   -- Good: Request exactly what you need
   NEAREST(genes, reference := peaks.interval, k := 3)

   -- Wasteful: Request more than needed
   NEAREST(genes, reference := peaks.interval, k := 100)

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

DuckDB IEJoin Dialect
~~~~~~~~~~~~~~~~~~~~~

For column-to-column ``INTERSECTS`` joins (INNER, SEMI, or ANTI) on
DuckDB, the ``dialect="duckdb"`` opt-in transpiles the join into a
per-chromosome dynamic-SQL pattern. Each per-chromosome subquery contains
only inequality predicates, which DuckDB plans through its range-join
family (``IE_JOIN`` or ``PIECEWISE_MERGE_JOIN``) — avoiding the
row-inflation and deduplication cost of the default binned equi-join plan
(the ``dialect=None`` path, controlled by ``intersects_bin_size``, which
expands each interval into ``UNNEST``-generated bins and deduplicates
the resulting matches).

.. code-block:: python

   from giql import transpile

   sql = transpile(
       """
       SELECT a.chrom, a.start, b.start
       FROM peaks a
       JOIN genes b ON a.interval INTERSECTS b.interval
       """,
       tables=["peaks", "genes"],
       dialect="duckdb",
   )
   # sql is a two-statement script — see the next example for execution.

The output is a multi-statement string of the form::

   SET VARIABLE __giql_iejoin_<token> = COALESCE((... string_agg per chromosome ...), '<empty schema>');
   SELECT ... FROM query(getvariable('__giql_iejoin_<token>')) AS __giql_iejoin_wrapper

The ``<token>`` is a per-call ``uuid4().hex`` (128 bits) suffix so the
``SET VARIABLE`` name is collision-resistant even when outputs from many
``transpile()`` calls are interleaved in a single DuckDB session
(session variables are global session state). The outer SELECT's
wrapper-relation alias is constant (``__giql_iejoin_wrapper``) because
it isn't user-visible and doesn't need to vary per call.

Inside the per-chromosome string-builder, chromosome names are emitted
via ``replace(chrom, '''', '''''')`` wrapped in single quotes, so
chromosome identifiers containing literal single quotes interpolate
safely into the dynamic SQL.

Because the output is a two-statement script that depends on session
state, execute it through a single DuckDB connection's ``.execute()`` —
DuckDB ≥1.4 (the version pinned by GIQL's development and CI
environment) accepts multi-statement strings in one call and returns
the result of the final statement. SQLAlchemy's ``text()``,
``pandas.read_sql_query``, and similar wrappers that split or rewrite
the string may drop the ``SET VARIABLE`` and produce empty or NULL
results.

.. code-block:: python

   import duckdb
   from giql import transpile

   conn = duckdb.connect()
   conn.execute("CREATE TABLE peaks (chrom VARCHAR, \"start\" INTEGER, \"end\" INTEGER)")
   conn.execute("CREATE TABLE genes (chrom VARCHAR, \"start\" INTEGER, \"end\" INTEGER)")
   sql = transpile(
       """
       SELECT a.chrom, a.start, b.start
       FROM peaks a
       JOIN genes b ON a.interval INTERSECTS b.interval
       """,
       tables=["peaks", "genes"],
       dialect="duckdb",
   )
   rows = conn.execute(sql).fetchall()

The dialect rewrites the whole query, so the supported shape is
``SELECT <qualified projections> FROM <base table> {INNER|SEMI|ANTI}
JOIN <base table> {ON|WHERE} <one column-to-column INTERSECTS>`` — the
JOIN and the column-to-column ``INTERSECTS`` are both required
(literal-range ``INTERSECTS`` against a single table falls through to
the default predicate generator).

**Join variants.** ``INNER JOIN`` is the default. ``SEMI JOIN`` returns
distinct left rows that have at least one overlapping match; ``ANTI
JOIN`` returns left rows with no overlapping match. Both restrict the
outer SELECT to left-side projections — any reference to the right
side (``b.col`` / ``b.*`` / aggregate over ``b.*``) raises
``ValueError``. ANTI uses a left-distinct chromosome partition rather
than the chromosome INTERSECT used by INNER / SEMI, so left rows on
chromosomes absent from the right table are preserved.

DuckDB plans INNER through the IE_JOIN / PIECEWISE_MERGE_JOIN
sort-merge family. SEMI and ANTI with inequality predicates plan
through ``BLOCKWISE_NL_JOIN`` instead — still avoiding the UNNEST
row-inflation of the binned plan, but not the IE_JOIN sort-merge fast
path. Expect substantial speedups vs. the binned plan for SEMI / ANTI
where the chromosome partition already filters most pairs; INNER gets
the largest speedup.

On top of the core shape the dialect also absorbs several common
decorations:

- **Outer modifiers.** ``ORDER BY`` / ``LIMIT`` / ``OFFSET`` /
  ``DISTINCT`` on the outer query are appended to the outer SELECT;
  column references in ``ORDER BY`` are rewritten to the wrapper
  relation's inner aliases.
- **Aggregates.** ``GROUP BY`` / ``HAVING`` with aggregate functions
  (``COUNT``, ``SUM``, ``MIN``, ``MAX``, ``AVG``, ``COUNT(DISTINCT
  a.col)``, ...) are appended to the outer SELECT. Aggregate arguments
  must be table-qualified (``COUNT(*)`` and ``COUNT(a.col)`` are both
  accepted).
- **Extra JOIN / WHERE predicates.** Additional non-INTERSECTS
  predicates ANDed onto the join ON or WHERE (e.g. ``a.score >
  b.score`` or ``WHERE a.score > 100``) are inlined into each
  per-chromosome subquery's ON, so DuckDB filters them inside each
  IEJoin candidate set. **Limitations:** the dialect peels extra
  predicates only across top-level ``AND`` connectives, so predicates
  that wrap the ``INTERSECTS`` in ``OR`` / ``NOT`` / parentheses still
  fall back. Predicates whose AND-tree residual contains a subquery,
  aggregate, or window function fall back too (DuckDB forbids window
  functions inside ``JOIN ON``, and subqueries / aggregates inside an
  IEJoin candidate set would either break the planner or produce
  semantically wrong results). If you hit the
  WHERE-INTERSECTS-plus-extra-JOIN-ON-predicate shape described in
  `#94 <https://github.com/abdenlab/giql/issues/94>`_ on the default
  ``dialect=None`` path, ``dialect="duckdb"`` is one workaround — the
  dialect inlines extras directly into each per-chromosome subquery
  and is unaffected by that bug.
- **Star projections (``a.*`` / ``b.*``).** Expanded to the genomic
  columns declared by the corresponding :class:`Table` config (chrom
  / start / end, plus strand when set). This narrows the result schema
  relative to the binned plan, which delegates star expansion to
  DuckDB and projects every base-table column. Users with additional
  non-genomic columns should list them explicitly.

The dialect splits unsupported shapes into two buckets. Soft-fallback
shapes route to the binned plan automatically and return correct
results; within those shapes the ``dialect`` kwarg is safe to set
without risk of silent incorrectness. Hard-error shapes (enumerated
further below) raise ``ValueError`` at transpile time — the dialect
deliberately refuses them rather than silently producing the wrong
SQL.

The soft-fallback shapes are:

- **Outer joins.** ``LEFT`` / ``RIGHT`` / ``FULL`` ``JOIN ... ON ...
  INTERSECTS ...`` falls back to the binned plan. The same applies
  when the join keeps its side modifier and the ``INTERSECTS`` lives in
  the top-level ``WHERE`` (e.g. ``LEFT JOIN ... ON TRUE WHERE
  a.interval INTERSECTS b.interval``).
- **NATURAL joins.** ``NATURAL JOIN`` falls back because the dialect
  cannot enumerate shared columns at transpile time (only the
  registered ``chrom`` / ``start`` / ``end`` / ``strand`` columns are
  known).
- **USING joins.** Single-column ``USING(<chrom_col>)`` admits (the
  per-chromosome partition is exactly the equi-join). Multi-column
  ``USING`` and ``USING(<non-chrom-col>)`` fall back; inline support
  is a documented follow-up.
- **Multiple INTERSECTS predicates.** Queries with more than one
  column-to-column ``INTERSECTS`` (anywhere in the AST, including
  subqueries) fall back.
- **INTERSECTS references unrecognized join sides.** A column-to-column
  INTERSECTS whose operands don't reference the FROM table's alias,
  or where the second operand's alias doesn't resolve to a registered
  JOIN target, falls back to the binned plan.
- **Self-joins.** Joining a table to itself falls back.
- **OR / NOT / paren-wrapped extra predicates.** See the Extra JOIN /
  WHERE predicates note above.
- **Subquery, aggregate, or window-function extra predicates.**
  Predicates whose AND-tree residual contains a subquery, an aggregate
  function, or a window function (``ROW_NUMBER() OVER (...)`` etc.)
  fall back rather than inline into the per-chromosome subquery's
  ``ON`` (DuckDB forbids window functions in ``ON`` clauses, and
  aggregates / subqueries inside an IEJoin candidate set would either
  break the planner or produce semantically wrong results).
- **Subquery inside GROUP BY / HAVING / ORDER BY.** A modifier clause
  containing a nested subquery (including ``EXISTS (SELECT ...)``)
  falls back to the binned plan, because the dialect's modifier-ref
  rewriter is not scope-aware and would corrupt the subquery's
  column scope. Subqueries hidden inside aggregate arguments in the
  SELECT list fall back for the same reason.
- **Top-level WITH clauses.** A query with a leading ``WITH`` falls
  back so the CTE survives.
- **More than two tables.** Any extra FROM-clause table (auxiliary
  comma-join or non-INTERSECTS ``JOIN``) falls back.
- **Non-base-table operands.** A subquery, CTE, or GIQL table-function
  (e.g. ``DISJOIN(genes)``) used as a join operand falls back.
- **DISTINCT ON (...).** Only plain ``DISTINCT`` is absorbed; the
  ``DISTINCT ON`` variant falls back.

Query-shape ``ValueError`` raises (the dialect deliberately refuses
these rather than silently miscompile):

- **Bare** ``SELECT *`` **and unqualified columns.** The dialect path
  needs to know which side each output column comes from. Use ``a.*``,
  ``b.*``, ``a.col``, or ``a.col AS x``.
- **Expression-form projections.** ``a.start + 1`` and other non-column
  expressions (other than aggregates) raise. Project the underlying
  column and compute in a wrapping query around the dialect's output,
  or omit ``dialect="duckdb"`` to use the binned plan.
- **Unknown table qualifier.** ``c.col`` when the only sides are ``a``
  and ``b``.
- **Unknown alias in GROUP BY / HAVING / ORDER BY.** A modifier-clause
  column qualified with a table alias outside the join's two sides
  (``c.col``) raises with a ``cannot resolve … in <CLAUSE>`` message
  naming the expected aliases.
- **Right-side reference in SEMI / ANTI joins.** ``b.col`` / ``b.*`` /
  aggregate over ``b.*`` in the outer SELECT, GROUP BY, HAVING, ORDER
  BY, or aggregate argument under ``SEMI JOIN`` / ``ANTI JOIN`` raises;
  SEMI / ANTI return left-side columns only.
- **Aggregate of unqualified column.** ``SUM(score)`` raises; use
  ``SUM(a.score)`` or ``SUM(b.score)``.
- **Star projection with a user alias.** ``a.* AS x`` raises (most
  engines reject star-with-alias outright); either drop the alias and
  use ``a.*``, or list the columns explicitly with per-column aliases.
- **Window aggregates and FILTER clauses.** ``SUM(a.score) OVER
  (...)`` and ``SUM(a.score) FILTER (WHERE ...)`` raise with a
  dedicated error; rewrite the projection or omit ``dialect="duckdb"``.
- **Aggregates inside expressions.** ``COUNT(*) * 2`` and similar
  arithmetic-over-aggregate projections raise; project the aggregate
  on its own and do arithmetic in a wrapping query.
- **Scalar subqueries in the SELECT list.** ``(SELECT SUM(x) FROM
  other_table) AS s`` raises; rewrite without the subquery or omit
  ``dialect="duckdb"``.
- **Catalog/schema-qualified extra predicates.** ``mycat.myschema.a.col``
  in an extra predicate raises; the alias rewriter would leave the
  catalog/schema qualifier intact in the inner subquery, where it
  references a different relation than intended.
- **Unqualified or unknown-aliased extra predicates.** ``WHERE strand
  = '+'`` (unqualified) or ``WHERE c.score > 0`` (alias outside the
  join's two sides) raise with an actionable message naming the
  expected aliases. The binned plan would either silently bind the
  column wrong or defer the error to the downstream engine, so the
  dialect surfaces the user mistake at transpile time instead.

Argument-validation ``ValueError`` raises (these fire before any AST
inspection and are independent of the query shape):

- **Mutually exclusive with** ``intersects_bin_size``. The two
  parameters cannot be combined; pass one or the other.
- **Unknown dialect string.** Anything other than ``"duckdb"`` or
  ``None`` raises with the offending value echoed.

Literal-range ``INTERSECTS`` (e.g. ``WHERE interval INTERSECTS
'chr1:1000-2000'``) is single-table and has no column-to-column join
to rewrite, so the dialect declines and the standard range-predicate
emission handles it identically to the ``dialect=None`` path.

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
