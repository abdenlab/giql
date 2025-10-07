GIQL - Genomic Interval Query Language
======================================

**GIQL** is a SQL dialect for genomic range queries with multi-database support.

GIQL extends SQL with spatial operators for genomic interval queries. It transpiles
to standard SQL that works across multiple database backends including DuckDB and SQLite.

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   quickstart
   api/index
   examples

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

       # Query with genomic operators
       result = engine.query("""
           SELECT * FROM variants
           WHERE position INTERSECTS 'chr1:1000-2000'
       """)

Features
--------

* **SQL-based**: Familiar SQL syntax with genomic extensions
* **Multi-backend**: Works with DuckDB, SQLite, and more
* **Spatial operators**: INTERSECTS, CONTAINS, WITHIN
* **Set quantifiers**: ANY, ALL for multi-range queries
* **Column-to-column joins**: Join tables on genomic position
* **Zero-copy**: Efficient in-memory operations with DuckDB

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
