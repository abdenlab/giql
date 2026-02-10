Execution
=========

How to use transpiled SQL
-------------------------

You can write queries in the GIQL dialect and execute them on any SQL-92 
compliant database or analytics engine, without needing native GIQL support.

With external database connections
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Use transpiled SQL with your own database connections:

.. code-block:: python

   import duckdb
   from giql import transpile

   sql = transpile(
       """
       SELECT * FROM variants
       WHERE interval INTERSECTS 'chr1:1000-2000'
       """,
       tables=["variants"],
   )

   conn = duckdb.connect("my_database.duckdb")
   result = conn.execute(sql).fetchall()
   conn.close()

With ORMs and query builders
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Integrate transpiled SQL with SQLAlchemy or other ORMs:

.. code-block:: python

   from sqlalchemy import create_engine, text
   from giql import transpile

   sql = transpile(
       """
       SELECT * FROM variants
       WHERE interval INTERSECTS 'chr1:1000-2000'
       """,
       tables=["variants"],
   )

   engine = create_engine("duckdb:///my_database.duckdb")
   with engine.connect() as conn:
       result = conn.execute(text(sql))
       for row in result:
           print(row)

Building SQL pipelines
~~~~~~~~~~~~~~~~~~~~~~

Use transpilation in data pipelines:

.. code-block:: python

   from giql import transpile

   def build_intersection_query(table_a, table_b, region):
       """Generate SQL for intersection query."""
       return transpile(
           f"""
           SELECT a.*, b.name
           FROM {table_a} a
           JOIN {table_b} b ON a.interval INTERSECTS b.interval
           WHERE a.interval INTERSECTS '{region}'
           """,
           tables=[table_a, table_b],
       )

   # Use in pipeline
   sql = build_intersection_query("variants", "genes", "chr1:1000000-2000000")
   # Execute sql with your preferred method

Saving queries
~~~~~~~~~~~~~~

Save transpiled SQL for documentation or reuse:

.. code-block:: python

   import duckdb
   from giql import transpile

   sql = transpile(
       """
       SELECT * FROM variants
       WHERE interval INTERSECTS 'chr1:1000-2000'
       """,
       tables=["variants"],
   )

   with open("query.sql", "w") as f:
       f.write(sql)

   # Later, execute saved SQL
   with open("query.sql") as f:
       sql = f.read()

   conn = duckdb.connect("database.duckdb")
   result = conn.execute(sql).fetchall()

Parameterized queries
~~~~~~~~~~~~~~~~~~~~~

Build queries with parameters:

.. code-block:: python

   from giql import transpile

   def query_region(chrom, start, end):
       """Transpile a parameterized region query."""
       region = f"{chrom}:{start}-{end}"
       return transpile(
           f"""
           SELECT * FROM variants
           WHERE interval INTERSECTS '{region}'
           """,
           tables=["variants"],
       )

   # Use with different regions
   sql = query_region("chr1", 1000000, 2000000)
   sql = query_region("chr2", 5000000, 6000000)

Dynamic query building
~~~~~~~~~~~~~~~~~~~~~~

Build queries programmatically:

.. code-block:: python

   from giql import transpile

   def build_multi_table_query(tables, target_region):
       """Build a query that unions results from multiple tables."""
       union_parts = []
       for table in tables:
           union_parts.append(f"""
               SELECT *, '{table}' AS source FROM {table}
               WHERE interval INTERSECTS '{target_region}'
           """)

       query = " UNION ALL ".join(union_parts)
       return transpile(query, tables=list(tables))
