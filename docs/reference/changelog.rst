Changelog
=========

All notable changes to GIQL are documented here.

The format is based on `Keep a Changelog <https://keepachangelog.com/>`_,
and this project adheres to `Semantic Versioning <https://semver.org/>`_.

.. contents::
   :local:
   :depth: 1

Unreleased
----------

*Changes in development, not yet released.*

Added
~~~~~

- Comprehensive documentation with operator reference
- Recipe-based examples for common patterns
- Bedtools migration guide
- Multi-backend support guide
- Performance optimization guide

Changed
~~~~~~~

- Documentation restructured with operator-first organization

0.1.0 - Initial Release
-----------------------

*Initial release of GIQL.*

Added
~~~~~

**Core Features:**

- SQL dialect for genomic interval queries
- Transpilation to standard SQL
- Multi-database backend support (DuckDB, SQLite)

**Spatial Operators:**

- ``INTERSECTS`` - Test range overlap
- ``CONTAINS`` - Test containment
- ``WITHIN`` - Test if range is within another

**Distance Operators:**

- ``DISTANCE`` - Calculate genomic distance
- ``NEAREST`` - K-nearest neighbor queries

**Aggregation Operators:**

- ``CLUSTER`` - Assign cluster IDs to overlapping intervals
- ``MERGE`` - Combine overlapping intervals

**Set Quantifiers:**

- ``ANY`` - Match any of multiple ranges
- ``ALL`` - Match all of multiple ranges

**API:**

- ``GIQLEngine`` - Main engine class
- ``execute()`` - Execute GIQL queries
- ``transpile()`` - Convert GIQL to SQL
- ``register_table_schema()`` - Register table schemas
- ``load_csv()`` - Load CSV files

Version History
---------------

.. list-table::
   :header-rows: 1
   :widths: 20 20 60

   * - Version
     - Date
     - Highlights
   * - 0.1.0
     - TBD
     - Initial release with core operators
