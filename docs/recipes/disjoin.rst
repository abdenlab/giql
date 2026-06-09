Disjoining Intervals
====================

This section covers patterns for splitting intervals at breakpoints using
GIQL's ``DISJOIN`` operator -- partitioning a set into non-overlapping
segments and re-tiling features against a reference grid.

See :doc:`../dialect/set-operators` for the operator's reference -- syntax,
parameter semantics, and the reserved output column / identifier conventions
the recipes below assume.

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

Disjoin a Filtered or Canonicalised Subset
------------------------------------------

Disjoin Only the Rows You Want
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Build the input to ``DISJOIN`` inline with a ``WITH`` clause, skipping the
need to register the filtered subset as a separate table:

.. code-block:: sql

   WITH filtered AS (
       SELECT chrom, "start", "end", name
       FROM features
       WHERE chrom = 'chr1'
   )
   SELECT DISTINCT disjoin_chrom, disjoin_start, disjoin_end
   FROM DISJOIN(filtered)
   ORDER BY disjoin_start

**Use case:** Partition only one chromosome (or any other inline filter) and
let the operator see exactly that subset, without registering a derived table
or materialising it on disk.

Canonicalise a Non-Default Schema
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

A CTE target is assumed to expose canonical 0-based half-open ``chrom`` /
``start`` / ``end`` columns. When the underlying table uses custom column
names or a non-default coordinate system, alias and convert inline:

.. code-block:: sql

   WITH canonical AS (
       SELECT seqid AS chrom, lo - 1 AS "start", hi AS "end", name
       FROM one_based_features
   )
   SELECT DISTINCT disjoin_chrom, disjoin_start, disjoin_end
   FROM DISJOIN(canonical)

Additional columns (here, ``name``) may be projected through the CTE and
will pass through to the DISJOIN output rows alongside the canonical three.

**Use case:** Use a 1-based-closed source table without registering a custom
``Table`` config -- the user takes responsibility for the canonical layout
the CTE-target contract requires.

.. warning::

   The CTE-target contract is unvalidated. If the CTE does not actually
   expose 0-based half-open ``chrom`` / ``start`` / ``end`` columns -- wrong
   names, wrong coordinate system, wrong types -- results will be either
   silently wrong or surfaced as an opaque engine error ("column not
   found"). There is no GIQL-level diagnostic.

Chain DISJOIN Through a CTE
~~~~~~~~~~~~~~~~~~~~~~~~~~~

``DISJOIN`` emits its appended ``disjoin_*`` columns in the *target*
table's coordinate system and passes the target's own ``chrom`` /
``start`` / ``end`` (and any other columns) through unchanged. When the
target is a registered non-canonical table, those passthrough columns
carry the table's native encoding into the CTE -- and the CTE-target
contract on a second ``DISJOIN`` assumes canonical 0-based half-open.
Canonicalise inside the CTE before re-disjoining:

.. code-block:: sql

   WITH segments AS (
       SELECT
           disjoin_chrom AS chrom,
           disjoin_start - 1 AS "start",          -- 1-based -> 0-based
           disjoin_end AS "end",                  -- closed -> half-open
           name
       FROM DISJOIN(one_based_features)
   )
   SELECT DISTINCT disjoin_chrom, disjoin_start, disjoin_end
   FROM DISJOIN(segments)

**Use case:** Two-stage partitioning -- e.g. disjoin features within
sample, then disjoin across samples -- when the source table is not
0-based half-open. Project only the canonicalised ``chrom`` / ``start`` /
``end`` (plus any passthrough columns) through the CTE; do **not**
re-export the raw ``start`` / ``end`` from the inner table, which still
carry the original encoding.

.. warning::

   If the inner ``DISJOIN``'s target is canonical 0-based half-open --
   either a CTE itself, or a registered table whose ``coordinate_system``
   / ``interval_type`` resolve to canonical -- chained ``DISJOIN`` is
   safe without any canonicalisation step. The trap is specifically
   when chaining across a non-canonical registered table.

Coming from Bedtools?
---------------------

``DISJOIN`` has no direct bedtools equivalent -- no single bedtools command
splits intervals at breakpoints. See the :doc:`bedtools-migration` guide for
the GIQL operators that do map onto bedtools commands.
