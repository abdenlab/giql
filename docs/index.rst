GIQL - Genomic Interval Query Language
======================================

**GIQL** is a SQL dialect for genomic range queries with multi-database support.

GIQL extends SQL with spatial operators for genomic interval queries. It transpiles
to standard SQL that works across multiple database backends including DuckDB and SQLite.

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

   reference/operator-matrix
   reference/syntax-reference
   reference/changelog
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
           genomic_column="position",
       )

       # Query with genomic operators (returns cursor for streaming)
       cursor = engine.execute("""
           SELECT * FROM variants
           WHERE position INTERSECTS 'chr1:1000-2000'
       """)

       # Process results
       for row in cursor:
           print(row)

       # Or just transpile to SQL without executing
       sql = engine.transpile("""
           SELECT * FROM variants
           WHERE position INTERSECTS 'chr1:1000-2000'
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
   WHERE position INTERSECTS 'chr1:1000-2000'

   -- Find containing/contained features
   WHERE gene.position CONTAINS variant.position

**Distance and Proximity:**

.. code-block:: sql

   -- Calculate distance between intervals
   SELECT DISTANCE(a.position, b.position) AS dist

   -- Find k-nearest neighbors
   FROM peaks CROSS JOIN LATERAL NEAREST(genes, reference=peaks.position, k=5)

**Aggregation:**

.. code-block:: sql

   -- Cluster overlapping intervals
   SELECT *, CLUSTER(position) AS cluster_id FROM features

   -- Merge overlapping intervals
   SELECT MERGE(position) FROM features

**Set Quantifiers:**

.. code-block:: sql

   -- Match any of multiple regions
   WHERE position INTERSECTS ANY('chr1:1000-2000', 'chr2:5000-6000')

See :doc:`operators/index` for complete operator documentation.

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
