Genomic Interval Query Language (GIQL)
======================================

.. toctree::
   :hidden:

   GIQL <self>
   quickstart

**GIQL** (pronounced "JEE-quel", rhymes with "sequel") is an extended SQL dialect that allows you to declaratively express genomic interval operations.

See the :doc:`quickstart` to get started.
   
Dialect
-------
GIQL extends the SQL query language with dedicated constructs for 
common tasks, allowing you to declare *what* you want to compute rather 
than how. Whether you're filtering variants by genomic region, finding 
overlapping features, or calculating distances between intervals, GIQL 
makes these operations intuitive and portable.

See the :doc:`GIQL dialect <dialect/index>` docs.

.. toctree::
   :hidden:
   :caption: Dialect

   Overview <dialect/index>
   dialect/spatial-operators
   dialect/distance-operators
   dialect/aggregation-operators
   dialect/quantifiers

Transpilation
-------------
The ``giql`` package *transpiles* queries written in GIQL to regular SQL
for use in existing database systems (like SQLite or PostgreSQL), data
warehouses, or analytics engines (like Polars and DuckDB).

See the :doc:`GIQL transpiler <transpilation/index>` docs.

.. toctree::
   :hidden:
   :caption: Transpilation

   Overview <transpilation/index>
   transpilation/schema-mapping
   transpilation/execution
   transpilation/performance
   transpilation/api-reference


Learn more
----------
See the following :doc:`recipes <recipes/index>` to learn how to use GIQL effectively.

.. toctree::
   :hidden:
   :caption: Recipes

   Overview <recipes/index>
   recipes/intersect
   recipes/distance
   recipes/clustering
   recipes/advanced
   recipes/bedtools-migration


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
