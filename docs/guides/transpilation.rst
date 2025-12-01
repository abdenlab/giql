Transpilation Guide
===================

GIQL works by transpiling genomic queries into standard SQL. This guide explains
how transpilation works, how to debug query generation, and how to use transpiled
SQL with external tools.

.. contents::
   :local:
   :depth: 2

How Transpilation Works
-----------------------

The Transpilation Process
~~~~~~~~~~~~~~~~~~~~~~~~~

When you write a GIQL query:

.. code-block:: sql

   SELECT * FROM variants WHERE position INTERSECTS 'chr1:1000-2000'

GIQL performs these steps:

1. **Parse**: Parse the SQL to identify GIQL-specific operators
2. **Expand**: Replace genomic operators with standard SQL predicates
3. **Generate**: Produce SQL for the target database dialect

The result is standard SQL:

.. code-block:: sql

   SELECT * FROM variants
   WHERE chromosome = 'chr1' AND start_pos < 2000 AND end_pos > 1000

Operator Expansion
~~~~~~~~~~~~~~~~~~

Each GIQL operator expands to specific SQL patterns:

**INTERSECTS** expands to range overlap checks:

.. code-block:: sql

   -- GIQL
   a.position INTERSECTS b.position

   -- SQL (same chromosome, overlapping ranges)
   a.chromosome = b.chromosome
   AND a.start_pos < b.end_pos
   AND a.end_pos > b.start_pos

**CONTAINS** expands to containment checks:

.. code-block:: sql

   -- GIQL
   a.position CONTAINS b.position

   -- SQL
   a.chromosome = b.chromosome
   AND a.start_pos <= b.start_pos
   AND a.end_pos >= b.end_pos

**DISTANCE** expands to gap calculations:

.. code-block:: sql

   -- GIQL
   DISTANCE(a.position, b.position)

   -- SQL (simplified)
   CASE
       WHEN a.chromosome != b.chromosome THEN NULL
       WHEN a.end_pos <= b.start_pos THEN b.start_pos - a.end_pos
       WHEN b.end_pos <= a.start_pos THEN a.start_pos - b.end_pos
       ELSE 0
   END

Using the Transpile Method
--------------------------

Basic Transpilation
~~~~~~~~~~~~~~~~~~~

Use ``transpile()`` to see generated SQL without executing:

.. code-block:: python

   from giql import GIQLEngine

   with GIQLEngine(target_dialect="duckdb") as engine:
       engine.register_table_schema(
           "variants",
           {
               "chromosome": "VARCHAR",
               "start_pos": "BIGINT",
               "end_pos": "BIGINT",
           },
           genomic_column="position",
       )

       sql = engine.transpile("""
           SELECT * FROM variants
           WHERE position INTERSECTS 'chr1:1000-2000'
       """)

       print(sql)
       # Output: SELECT * FROM variants
       #         WHERE chromosome = 'chr1' AND start_pos < 2000 AND end_pos > 1000

Transpiling Complex Queries
~~~~~~~~~~~~~~~~~~~~~~~~~~~

Transpilation works with all GIQL features:

.. code-block:: python

   # Join query
   sql = engine.transpile("""
       SELECT v.*, g.name AS gene_name
       FROM variants v
       JOIN genes g ON v.position INTERSECTS g.position
       WHERE v.quality >= 30
   """)
   print(sql)

   # NEAREST query
   sql = engine.transpile("""
       SELECT peaks.name, nearest.name, nearest.distance
       FROM peaks
       CROSS JOIN LATERAL NEAREST(genes, reference=peaks.position, k=5) AS nearest
   """)
   print(sql)

   # Aggregation query
   sql = engine.transpile("""
       SELECT MERGE(position), COUNT(*) AS count
       FROM features
   """)
   print(sql)

Debugging with Transpilation
----------------------------

Understanding Query Expansion
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Use transpilation to understand what GIQL does:

.. code-block:: python

   # See how ANY quantifier expands
   sql = engine.transpile("""
       SELECT * FROM variants
       WHERE position INTERSECTS ANY('chr1:1000-2000', 'chr2:5000-6000')
   """)
   print(sql)
   # Shows the OR conditions for each range

   # See how join conditions expand
   sql = engine.transpile("""
       SELECT a.*, b.name
       FROM features_a a
       JOIN features_b b ON a.position INTERSECTS b.position
   """)
   print(sql)
   # Shows the full range comparison predicates

Verbose Mode
~~~~~~~~~~~~

Enable verbose mode for detailed transpilation information:

.. code-block:: python

   with GIQLEngine(target_dialect="duckdb", verbose=True) as engine:
       engine.register_table_schema("variants", {...}, genomic_column="position")

       # Transpilation will print detailed information
       sql = engine.transpile("""
           SELECT * FROM variants
           WHERE position INTERSECTS 'chr1:1000-2000'
       """)

       # Execution also shows transpilation details
       cursor = engine.execute("""
           SELECT * FROM variants
           WHERE position INTERSECTS 'chr1:1000-2000'
       """)

Troubleshooting Transpilation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Query not expanding correctly:**

.. code-block:: python

   # Check that schema is registered
   sql = engine.transpile("SELECT * FROM variants WHERE position INTERSECTS 'chr1:1000-2000'")
   if "position INTERSECTS" in sql:
       print("Schema not registered for 'variants' table")

**Wrong column names in output:**

.. code-block:: python

   # Verify column mapping
   engine.register_table_schema(
       "variants",
       {...},
       genomic_column="position",
       chromosome_column="chrom",      # Check these match your table
       start_column="start",
       end_column="end",
   )

Comparing Dialects
------------------

Same Query, Different SQL
~~~~~~~~~~~~~~~~~~~~~~~~~

See how the same query translates for different backends:

.. code-block:: python

   query = """
       SELECT * FROM variants
       WHERE position INTERSECTS 'chr1:1000-2000'
   """

   schema = {
       "chromosome": "VARCHAR",
       "start_pos": "BIGINT",
       "end_pos": "BIGINT",
   }

   # DuckDB
   with GIQLEngine(target_dialect="duckdb") as engine:
       engine.register_table_schema("variants", schema, genomic_column="position")
       print("DuckDB SQL:")
       print(engine.transpile(query))
       print()

   # SQLite
   with GIQLEngine(target_dialect="sqlite") as engine:
       engine.register_table_schema("variants", schema, genomic_column="position")
       print("SQLite SQL:")
       print(engine.transpile(query))

Dialect-Specific Differences
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Some queries may generate different SQL for different dialects:

- String functions may use different names
- Type casting syntax may vary
- Window function support may differ

GIQL handles these differences automatically, but understanding them helps
when debugging or integrating with external tools.

Using Transpiled SQL Externally
-------------------------------

With External Database Connections
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Use transpiled SQL with your own database connections:

.. code-block:: python

   import duckdb

   # Generate SQL using GIQL
   with GIQLEngine(target_dialect="duckdb") as giql_engine:
       giql_engine.register_table_schema("variants", {...}, genomic_column="position")
       sql = giql_engine.transpile("""
           SELECT * FROM variants
           WHERE position INTERSECTS 'chr1:1000-2000'
       """)

   # Execute with external connection
   conn = duckdb.connect("my_database.duckdb")
   result = conn.execute(sql).fetchall()
   conn.close()

With ORMs and Query Builders
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Integrate transpiled SQL with SQLAlchemy or other ORMs:

.. code-block:: python

   from sqlalchemy import create_engine, text

   # Generate SQL
   with GIQLEngine(target_dialect="duckdb") as giql_engine:
       giql_engine.register_table_schema("variants", {...}, genomic_column="position")
       sql = giql_engine.transpile("""
           SELECT * FROM variants
           WHERE position INTERSECTS 'chr1:1000-2000'
       """)

   # Execute with SQLAlchemy
   sa_engine = create_engine("duckdb:///my_database.duckdb")
   with sa_engine.connect() as conn:
       result = conn.execute(text(sql))
       for row in result:
           print(row)

Building SQL Pipelines
~~~~~~~~~~~~~~~~~~~~~~

Use transpilation in data pipelines:

.. code-block:: python

   def build_intersection_query(table_a, table_b, region):
       """Generate SQL for intersection query."""
       with GIQLEngine(target_dialect="duckdb") as engine:
           engine.register_table_schema(table_a, {...}, genomic_column="position")
           engine.register_table_schema(table_b, {...}, genomic_column="position")

           return engine.transpile(f"""
               SELECT a.*, b.name
               FROM {table_a} a
               JOIN {table_b} b ON a.position INTERSECTS b.position
               WHERE a.position INTERSECTS '{region}'
           """)

   # Use in pipeline
   sql = build_intersection_query("variants", "genes", "chr1:1000000-2000000")
   # Execute sql with your preferred method

Saving Queries
~~~~~~~~~~~~~~

Save transpiled SQL for documentation or reuse:

.. code-block:: python

   # Generate and save SQL
   with GIQLEngine(target_dialect="duckdb") as engine:
       engine.register_table_schema("variants", {...}, genomic_column="position")

       sql = engine.transpile("""
           SELECT * FROM variants
           WHERE position INTERSECTS 'chr1:1000-2000'
       """)

       with open("query.sql", "w") as f:
           f.write(sql)

   # Later, execute saved SQL
   with open("query.sql") as f:
       sql = f.read()

   conn = duckdb.connect("database.duckdb")
   result = conn.execute(sql).fetchall()

Advanced Transpilation
----------------------

Parameterized Queries
~~~~~~~~~~~~~~~~~~~~~

Build queries with parameters:

.. code-block:: python

   def query_region(engine, chrom, start, end):
       """Query a parameterized region."""
       region = f"{chrom}:{start}-{end}"
       return engine.execute(f"""
           SELECT * FROM variants
           WHERE position INTERSECTS '{region}'
       """)

   # Use with different regions
   cursor = query_region(engine, "chr1", 1000000, 2000000)
   cursor = query_region(engine, "chr2", 5000000, 6000000)

Dynamic Query Building
~~~~~~~~~~~~~~~~~~~~~~

Build queries programmatically:

.. code-block:: python

   def build_multi_table_query(tables, target_region):
       """Build a query that unions results from multiple tables."""
       union_parts = []
       for table in tables:
           union_parts.append(f"""
               SELECT *, '{table}' AS source FROM {table}
               WHERE position INTERSECTS '{target_region}'
           """)

       query = " UNION ALL ".join(union_parts)
       return engine.transpile(query)

Inspecting the AST
~~~~~~~~~~~~~~~~~~

For advanced debugging, you can inspect the parsed query:

.. code-block:: python

   # GIQL uses sqlglot internally
   # The transpiled SQL shows the final result
   sql = engine.transpile("SELECT * FROM variants WHERE position INTERSECTS 'chr1:1000-2000'")

   # For deep debugging, examine the generated SQL structure
   print(sql)
