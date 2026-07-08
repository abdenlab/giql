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

By default (``dialect=None`` / ``"datafusion"``) a column-to-column
``INTERSECTS`` join emits the **naive overlap predicate** — a plain ``ON
a.chrom = b.chrom AND a.start < b.end AND b.start < a.end`` condition — and
lets the engine's own optimizer plan it. On DuckDB and DataFusion this
becomes a hash join keyed on ``chrom`` with the two position inequalities
as a residual join filter, correct for both inner and outer joins with no
special handling. This is the lowest-common-denominator plan: standard SQL
that every target runs.

For column-to-column ``INTERSECTS`` joins (INNER, SEMI, or ANTI) on
DuckDB, the ``dialect="duckdb"`` opt-in instead transpiles the join into a
per-chromosome dynamic-SQL pattern. Each per-chromosome subquery contains
only inequality predicates, which DuckDB plans through its range-join
family (``IE_JOIN`` or ``PIECEWISE_MERGE_JOIN``). Shapes this path declines
fall through to the naive predicate above.

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
``SELECT <qualified non-star projections> FROM <base table>
{INNER|SEMI|ANTI} JOIN <base table> {ON|WHERE} <one column-to-column
INTERSECTS>`` — the
JOIN and the column-to-column ``INTERSECTS`` are both required
(literal-range ``INTERSECTS`` against a single table falls through to
the default predicate generator).

**Join variants.** ``INNER JOIN`` is the default. ``SEMI JOIN`` returns
distinct left rows that have at least one overlapping match; ``ANTI
JOIN`` returns left rows with no overlapping match. Both restrict the
outer SELECT to left-side projections — any ``b.col`` reference (or an
aggregate over ``b.col``) raises ``ValueError``; a right-side ``b.*``
declines to the naive-predicate plan with every other star (which then
rejects the out-of-scope right table at bind time). ANTI uses a
left-distinct chromosome partition rather than the chromosome INTERSECT
used by INNER / SEMI, so left rows on chromosomes absent from the right
table are preserved.

DuckDB plans all three variants through its ``IE_JOIN`` /
``PIECEWISE_MERGE_JOIN`` sort-merge range-join family. INNER emits a
per-chromosome inner join directly. SEMI and ANTI are emitted as a
correlated ``WHERE EXISTS`` / ``WHERE NOT EXISTS`` subquery over the
per-chromosome right partition, because DuckDB plans a bare ``SEMI
JOIN`` / ``ANTI JOIN`` inequality overlap as a quadratic
``BLOCKWISE_NL_JOIN`` (a nested loop) rather than the range-join fast
path — the correlated form reaches the fast path while preserving exact
semantics (#208). Expect large speedups over the naive predicate for all
three variants at scale.

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

The dialect splits unsupported shapes into two buckets. Soft-fallback
shapes route to the naive-predicate plan automatically and return correct
results; within those shapes the ``dialect`` kwarg is safe to set
without risk of silent incorrectness. Hard-error shapes (enumerated
further below) raise ``ValueError`` at transpile time — the dialect
deliberately refuses them rather than silently producing the wrong
SQL.

The soft-fallback shapes are:

- **Outer joins.** ``LEFT`` / ``RIGHT`` / ``FULL`` ``JOIN ... ON ...
  INTERSECTS ...`` falls back to the naive-predicate plan. The same applies
  when the join keeps its side modifier and the ``INTERSECTS`` lives in
  the top-level ``WHERE`` (e.g. ``LEFT JOIN ... ON TRUE WHERE
  a.interval INTERSECTS b.interval``).
- **SEMI / ANTI join with the INTERSECTS in the WHERE.** A ``SEMI`` /
  ``ANTI`` join whose column-to-column ``INTERSECTS`` sits in the
  top-level ``WHERE`` rather than its own ``ON`` (e.g. ``ANTI JOIN ... ON
  TRUE WHERE a.interval INTERSECTS b.interval``) falls back. The right
  table is out of scope in the ``WHERE`` after a left-only join, so the
  reference plans reject it with a binder error; the dialect declines so
  it surfaces that same error instead of relocating the predicate into
  the join and inventing anti/semi-overlap results. The idiomatic form
  with the ``INTERSECTS`` in the join ``ON`` is unaffected.
- **Star projections.** Any star in the top-level SELECT list — bare
  ``*``, ``a.*``, ``b.*``, or ``a.* AS x`` — falls back. Building the
  dialect's per-chromosome dynamic SQL requires enumerating the star's
  columns at transpile time, but only the registered genomic columns
  (``chrom`` / ``start`` / ``end`` / ``strand``) are known — arbitrary
  user columns are invisible. The naive-predicate plan expands the real
  star against DuckDB's live schema, so declining keeps the projection
  identical across backends rather than silently narrowing it.
- **Projections the IEJoin cannot rebuild.** The dialect's projection
  rebuild handles a qualified column (``a.col``) or a plain aggregate
  over qualified columns (``SUM(a.score)`` / ``COUNT(*)``). Any other
  SELECT-list shape the naive plan compiles for free — an expression
  (``a.start + 1``), a window aggregate (``SUM(a.score) OVER (...)``), a
  ``FILTER`` clause, a scalar subquery, an aggregate nested in an
  expression (``COUNT(*) * 2``), a bare literal, or a star nested in an
  aggregate argument (``COUNT(a.*)`` / ``MIN(COLUMNS(*))``) — falls back,
  so ``dialect="duckdb"`` stays consistent with every other backend
  instead of hard-erroring or, for the aggregate-nested star,
  miscompiling. A projection whose column the rebuild cannot attribute to
  a join side (unqualified, an unknown table, or the right side under a
  SEMI / ANTI join) raises at transpile time instead (see the hard-error
  list below) — the unknown-table and right-side cases are rejected by the
  naive plan too, and a bare unqualified column, though it may be
  naive-valid when unambiguous, has no side the dialect can infer without a
  live schema.
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
  JOIN target, falls back to the naive-predicate plan.
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
  falls back to the naive-predicate plan, because the dialect's modifier-ref
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

- **Unqualified columns.** A bare column reference (``SELECT score``)
  raises — the dialect path needs to know which side each output column
  comes from, and has no live schema to infer it (so an unambiguous
  unqualified column that the naive plan would accept still raises here).
  Use ``a.col`` or ``a.col AS x``. (Bare ``SELECT *`` does not raise; like
  every star it falls back to the naive-predicate plan — see the
  soft-fallback list above. An expression / window aggregate / ``FILTER``
  clause over an *unqualified* column raises the same "requires qualified
  projections" error — once qualified, the wrapper itself falls back to the
  naive plan rather than erroring. A scalar subquery always falls back, even
  over an unqualified inner column, since that column resolves against the
  subquery's own scope.)
- **Unknown table qualifier.** ``c.col`` when the only sides are ``a``
  and ``b``.
- **Unknown alias in GROUP BY / HAVING / ORDER BY.** A modifier-clause
  column qualified with a table alias outside the join's two sides
  (``c.col``) raises with a ``cannot resolve … in <CLAUSE>`` message
  naming the expected aliases.
- **Right-side reference in SEMI / ANTI joins.** ``b.col`` or an
  aggregate over ``b.col`` in the outer SELECT, GROUP BY, HAVING, ORDER
  BY, or aggregate argument under ``SEMI JOIN`` / ``ANTI JOIN`` raises;
  SEMI / ANTI return left-side columns only. (A right-side ``b.*`` falls
  back to the naive-predicate plan with every other star — see the
  soft-fallback list above.)
- **Aggregate of unqualified column.** ``SUM(score)`` raises; use
  ``SUM(a.score)`` or ``SUM(b.score)``. (A qualified aggregate the IEJoin
  cannot rebuild — ``SUM(a.score) OVER (...)``, ``SUM(a.score) FILTER
  (WHERE ...)``, ``COUNT(*) * 2``, ``COUNT(a.*)`` — falls back to the
  naive plan instead; see the soft-fallback list above.)
- **Catalog/schema-qualified extra predicates.** ``mycat.myschema.a.col``
  in an extra predicate raises; the alias rewriter would leave the
  catalog/schema qualifier intact in the inner subquery, where it
  references a different relation than intended.
- **Unqualified or unknown-aliased extra predicates.** ``WHERE strand
  = '+'`` (unqualified) or ``WHERE c.score > 0`` (alias outside the
  join's two sides) raise with an actionable message naming the
  expected aliases. The naive-predicate plan would either silently bind the
  column wrong or defer the error to the downstream engine, so the
  dialect surfaces the user mistake at transpile time instead.

The one argument-validation ``ValueError`` raise (it fires before any AST
inspection and is independent of the query shape):

- **Unknown dialect string.** A dialect that resolves to neither a built-in
  target (``None``, ``"duckdb"``, ``"datafusion"``) nor a custom target
  registered on the plugin hub raises with the offending value echoed.

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
