Disjoining Intervals
====================

This section covers patterns for splitting intervals at breakpoints using
GIQL's ``DISJOIN`` operator -- partitioning a set into non-overlapping
segments and re-tiling features against a reference grid.

Partition a Set of Intervals
----------------------------

Build a Disjoint Partition
~~~~~~~~~~~~~~~~~~~~~~~~~~

Split a set of intervals into the maximal set of non-overlapping
sub-intervals defined by its own breakpoints:

.. code-block:: sql

   SELECT DISTINCT disjoin_chrom, disjoin_start, disjoin_end
   FROM DISJOIN(features)
   ORDER BY disjoin_chrom, disjoin_start

**Use case:** Given ChIP-seq peak calls pooled from several samples, produce a
non-overlapping segment track -- the equivalent of Bioconductor's ``disjoin()``
-- so downstream signal or count aggregates never double-count a base covered
by two overlapping sample peaks.

Track Each Segment's Parent
~~~~~~~~~~~~~~~~~~~~~~~~~~~

Keep the parent feature alongside each sub-interval:

.. code-block:: sql

   SELECT name, disjoin_start, disjoin_end
   FROM DISJOIN(features)
   ORDER BY name, disjoin_start

**Use case:** See how each original feature was fragmented and which segment
came from which parent.

Split Against a Reference
-------------------------

Split Features Against a Mask
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Cut target features at the boundaries of a reference (mask) set, keeping only
the pieces the mask covers:

.. code-block:: sql

   SELECT name, disjoin_start, disjoin_end
   FROM DISJOIN(features, reference := mask)

**Use case:** Restrict ATAC-seq or gene features to a set of callable (or
otherwise interesting) mask regions, splitting them at the mask boundaries so
no reported piece straddles a mask edge.

Re-tile Against a Uniform Grid
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Pass a generated set of fixed-width bins as the reference:

.. code-block:: sql

   WITH bins AS (
       SELECT 'chr1' AS chrom, x AS start, x + 1000 AS "end"
       FROM range(0, 250000000, 1000) AS t(x)
   )
   SELECT * FROM DISJOIN(features, reference := bins)

**Use case:** Break features onto a uniform coordinate grid -- for example to
build a fixed-width binned coverage matrix -- so each piece falls within a
single bin.

.. note::

   ``range()`` is DuckDB-specific syntax for generating the bin grid; other
   engines need their own generator. The grid must also span every chromosome
   present in ``features``, or features on an uncovered chromosome are dropped.

Coming from Bedtools?
---------------------

``DISJOIN`` has no direct bedtools equivalent -- no single bedtools command
splits intervals at breakpoints. See the :doc:`bedtools-migration` guide for
the GIQL operators that do map onto bedtools commands.
