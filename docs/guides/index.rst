Guides
======

Task-oriented guides for working with GIQL. These guides cover common workflows
and best practices for using GIQL effectively.

.. toctree::
   :maxdepth: 2

   schema-mapping
   multi-backend
   performance
   transpilation

Guide Overview
--------------

:doc:`schema-mapping`
   Learn how to configure GIQL to work with your genomic data, including
   registering table schemas and mapping logical genomic columns.

:doc:`multi-backend`
   Understand GIQL's multi-database support and how to work with different
   backends like DuckDB, SQLite, and PostgreSQL.

:doc:`performance`
   Optimize your GIQL queries for better performance with indexing strategies,
   query patterns, and backend-specific tips.

:doc:`transpilation`
   Understand how GIQL translates queries to SQL, debug query generation,
   and integrate transpiled SQL with external tools.
