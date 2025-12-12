GIQL - Genomic Interval Query Language
======================================

**GIQL** is a SQL dialect for genomic range queries with multi-database support.

Genomic analysis often requires repetitive, complex SQL patterns to express simple operations like finding overlapping intervals or merging features. GIQL extends SQL with dedicated operators for these common tasks, so you can declaratively express *what* you want to compute without getting lost in SQL boilerplate. GIQL queries read naturally, even without SQL expertise - this clarity makes your analysis code easier to review and share. Best of all, GIQL queries work across DuckDB, SQLite, PostgreSQL, and other databases, so you're never locked into a specific engine and can choose the tool that fits your use case. Finally, GIQL operators follow established conventions from tools like bedtools, so the semantics are familiar and predictable.

.. toctree::
   :maxdepth: 2
   :caption: Getting Started

   quickstart

.. toctree::
   :maxdepth: 2
   :caption: Operator Reference

   operators/index

.. toctree::
   :maxdepth: 2
   :caption: Guides

   guides/index

.. toctree::
   :maxdepth: 2
   :caption: Recipes

   recipes/index

.. toctree::
   :maxdepth: 2
   :caption: Reference

   reference/syntax-reference
   api/index

Quick Start
-----------

Install GIQL:

.. code-block:: bash

   pip install giql

Basic usage:

.. code-block:: python

   from giql import GIQLEngine

   # Create engine with DuckDB backend
   with GIQLEngine(target_dialect="duckdb") as engine:
       # Load genomic data
       engine.load_csv("variants", "variants.csv")
       engine.register_table_schema(
           "variants",
           {
               "id": "INTEGER",
               "chromosome": "VARCHAR",
               "start_pos": "BIGINT",
               "end_pos": "BIGINT",
           },
           genomic_column="interval",
       )

       # Query with genomic operators (returns cursor for streaming)
       cursor = engine.execute("""
           SELECT * FROM variants
           WHERE interval INTERSECTS 'chr1:1000-2000'
       """)

       # Process results
       for row in cursor:
           print(row)

       # Or just transpile to SQL without executing
       sql = engine.transpile("""
           SELECT * FROM variants
           WHERE interval INTERSECTS 'chr1:1000-2000'
       """)
       print(sql)  # See the generated SQL

Features
--------

* **SQL-based**: Familiar SQL syntax with genomic extensions
* **Multi-backend**: Works with DuckDB, SQLite, and more
* **Spatial operators**: INTERSECTS, CONTAINS, WITHIN, DISTANCE, NEAREST
* **Aggregation operators**: CLUSTER, MERGE for combining intervals
* **Set quantifiers**: ANY, ALL for multi-range queries
* **Column-to-column joins**: Join tables on genomic position
* **Transpilation**: Convert GIQL to standard SQL for debugging or external use

Operators at a Glance
---------------------

**Spatial Relationships:**

.. code-block:: sql

   -- Find overlapping features
   WHERE interval INTERSECTS 'chr1:1000-2000'

   -- Find containing/contained features
   WHERE gene.interval CONTAINS variant.interval

**Distance and Proximity:**

.. code-block:: sql

   -- Calculate distance between intervals
   SELECT DISTANCE(a.interval, b.interval) AS dist

   -- Find k-nearest neighbors
   FROM peaks CROSS JOIN LATERAL NEAREST(genes, reference=peaks.interval, k=5)

**Aggregation:**

.. code-block:: sql

   -- Cluster overlapping intervals
   SELECT *, CLUSTER(interval) AS cluster_id FROM features

   -- Merge overlapping intervals
   SELECT MERGE(interval) FROM features

**Set Quantifiers:**

.. code-block:: sql

   -- Match any of multiple regions
   WHERE interval INTERSECTS ANY('chr1:1000-2000', 'chr2:5000-6000')

See :doc:`operators/index` for complete operator documentation.

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
