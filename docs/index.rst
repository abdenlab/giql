Genomic Interval Query Language (GIQL)
======================================

**GIQL** is an extended SQL dialect that allows you to declaratively express genomic interval operations.

Dialect
-------
GIQL extends the SQL query language with dedicated constructs for these 
common tasks, allowing you to declare *what* you want to compute rather 
than how. Whether you're filtering variants by genomic region, finding 
overlapping features, or calculating distances between intervals, GIQL 
makes these operations intuitive and portable.

.. toctree::
   :maxdepth: 1
   :caption: Dialect

   dialect/index
   dialect/syntax-reference

Transpilation
-------------
The ``giql`` package *transpiles* queries written in GIQL to regular SQL
for use in existing database systems and analytics engines.

.. toctree::
   :maxdepth: 1
   :caption: Transpilation

   transpilation/index
   transpilation/execution
   transpilation/api-reference


Learn more
----------
See the following guides to learn how to use GIQL effectively:

.. toctree::
   :maxdepth: 1
   :caption: Guides and Recipes

   guides/quickstart
   guides/index
   recipes/index


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
