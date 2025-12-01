GIQL Operators
==============

GIQL extends SQL with operators specifically designed for genomic interval queries.
These operators enable powerful spatial reasoning over genomic coordinates without
requiring complex SQL expressions.

Operators are organized by functionality:

.. contents::
   :local:
   :depth: 1

Spatial Relationship Operators
------------------------------

Test positional relationships between genomic ranges.

.. list-table::
   :header-rows: 1
   :widths: 20 50 30

   * - Operator
     - Description
     - Example
   * - :ref:`INTERSECTS <intersects-operator>`
     - Returns true when ranges overlap by at least one base pair
     - ``position INTERSECTS 'chr1:1000-2000'``
   * - :ref:`CONTAINS <contains-operator>`
     - Returns true when one range fully contains another
     - ``position CONTAINS 'chr1:1500'``
   * - :ref:`WITHIN <within-operator>`
     - Returns true when one range is fully within another
     - ``position WITHIN 'chr1:1000-5000'``

See :doc:`spatial-operators` for detailed documentation.

Distance and Proximity Operators
--------------------------------

Calculate distances and find nearest features.

.. list-table::
   :header-rows: 1
   :widths: 20 50 30

   * - Operator
     - Description
     - Example
   * - :ref:`DISTANCE <distance-operator>`
     - Calculate genomic distance between two intervals
     - ``DISTANCE(a.position, b.position)``
   * - :ref:`NEAREST <nearest-operator>`
     - Find k-nearest genomic features
     - ``NEAREST(genes, reference=peaks.position, k=5)``

See :doc:`distance-operators` for detailed documentation.

Aggregation Operators
---------------------

Combine and cluster genomic intervals.

.. list-table::
   :header-rows: 1
   :widths: 20 50 30

   * - Operator
     - Description
     - Example
   * - :ref:`CLUSTER <cluster-operator>`
     - Assign cluster IDs to overlapping intervals
     - ``CLUSTER(position) AS cluster_id``
   * - :ref:`MERGE <merge-operator>`
     - Combine overlapping intervals into unified regions
     - ``SELECT MERGE(position) FROM features``

See :doc:`aggregation-operators` for detailed documentation.

Set Quantifiers
---------------

Apply operators to multiple ranges simultaneously.

.. list-table::
   :header-rows: 1
   :widths: 20 50 30

   * - Quantifier
     - Description
     - Example
   * - :ref:`ANY <any-quantifier>`
     - Match if condition holds for any of the specified ranges
     - ``position INTERSECTS ANY('chr1:1000-2000', 'chr2:5000-6000')``
   * - :ref:`ALL <all-quantifier>`
     - Match if condition holds for all of the specified ranges
     - ``position CONTAINS ALL('chr1:1500', 'chr1:1600')``

See :doc:`quantifiers` for detailed documentation.

Operator Compatibility
----------------------

All operators work across supported database backends (DuckDB, SQLite, with PostgreSQL planned).
Each operator page includes a compatibility table showing backend support status.

.. toctree::
   :maxdepth: 2
   :hidden:

   spatial-operators
   distance-operators
   aggregation-operators
   quantifiers
