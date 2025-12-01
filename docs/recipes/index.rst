Recipes
=======

This section provides practical examples and patterns for common genomic analysis tasks
using GIQL. Each recipe focuses on a specific use case with ready-to-use query patterns.

.. contents::
   :local:
   :depth: 1

Getting Started with Recipes
----------------------------

All recipes assume you have set up a GIQL engine and registered your table schemas:

.. code-block:: python

   from giql import GIQLEngine

   with GIQLEngine(target_dialect="duckdb") as engine:
       # Load your data
       engine.load_csv("features_a", "file_a.bed")
       engine.load_csv("features_b", "file_b.bed")

       # Register schemas with genomic column mapping
       for table in ["features_a", "features_b"]:
           engine.register_table_schema(
               table,
               {
                   "chromosome": "VARCHAR",
                   "start_pos": "BIGINT",
                   "end_pos": "BIGINT",
                   "name": "VARCHAR",
                   "score": "FLOAT",
                   "strand": "VARCHAR",
               },
               genomic_column="position",
           )

       # Now run queries from the recipes below
       cursor = engine.execute("...")

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
