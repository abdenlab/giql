Disjoining Intervals
====================

This section covers patterns for splitting intervals at breakpoints using
GIQL's ``DISJOIN`` operator -- partitioning a set into non-overlapping
segments and re-tiling features against a reference grid.

Partition a Set of Intervals
-----------------------------

Build a Disjoint Partition
~~~~~~~~~~~~~~~~~~~~~~~~~~~

Split a set of intervals into the maximal set of non-overlapping
sub-intervals defined by its own breakpoints:

.. code-block:: sql

   SELECT DISTINCT disjoin_chrom, disjoin_start, disjoin_end
   FROM DISJOIN(features)
   ORDER BY disjoin_chrom, disjoin_start

**Use case:** Produce a non-overlapping segment track -- the equivalent of
Bioconductor's ``disjoin()`` -- that downstream queries can aggregate without
double-counting overlaps.

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

**Use case:** Restrict features to mask regions while splitting them at mask
boundaries so no piece straddles a mask edge.

Re-tile Against a Uniform Grid
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Pass a generated set of fixed-width bins as the reference:

.. code-block:: sql

   WITH bins AS (
       SELECT 'chr1' AS chrom, x AS start, x + 1000 AS "end"
       FROM range(0, 250000000, 1000) AS t(x)
   )
   SELECT * FROM DISJOIN(features, reference := bins)

**Use case:** Break features onto a uniform coordinate grid so each piece
falls within a single bin.

Coming from Bedtools?
---------------------

``DISJOIN`` has no single bedtools equivalent. The self-mode partition is
closest to the breakpoints implied by ``bedtools merge`` re-split at every
input boundary; splitting against a reference is closest to ``bedtools
intersect`` combined with the reference's boundaries. See the
:doc:`bedtools-migration` guide for related operations.
