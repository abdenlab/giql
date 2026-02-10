Recipes
=======

This section provides practical examples and patterns for common genomic analysis tasks
using GIQL. Each recipe focuses on a specific use case with ready-to-use query patterns.

.. contents::
   :local:
   :depth: 1

Getting Started with Recipes
----------------------------

All recipes show GIQL queries that you can transpile and execute on your database.
Setup:

.. code-block:: python

   from giql import transpile

   # Transpile any GIQL query to SQL
   sql = transpile(
       "... GIQL query from the recipes below ...",
       tables=["features_a", "features_b"],
   )

   # Then execute the SQL on your database connection
   # e.g., conn.execute(sql)

Recipe Categories
-----------------

:doc:`intersect-queries`
   Finding overlapping features, filtering by overlap, counting overlaps,
   strand-specific operations, and join patterns.

:doc:`distance-queries`
   Calculating distances between features, finding nearest neighbors,
   distance-constrained searches, and directional queries.

:doc:`clustering-queries`
   Clustering overlapping intervals, distance-based clustering,
   merging intervals, and aggregating cluster statistics.

:doc:`advanced-queries`
   Multi-range matching, complex filtering with joins, aggregate statistics,
   window expansions, and multi-table queries.

Coming from Bedtools?
---------------------

If you're familiar with bedtools and want to replicate specific commands in GIQL,
see the :doc:`bedtools-migration` guide for a complete mapping of bedtools
operations to GIQL equivalents.

.. toctree::
   :maxdepth: 2
   :hidden:

   intersect-queries
   distance-queries
   clustering-queries
   advanced-queries
   bedtools-migration
