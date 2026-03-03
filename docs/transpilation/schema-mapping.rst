Schema Mapping
==============

This guide explains how to configure GIQL to work with your genomic data by
defining table configurations that map logical genomic columns to physical columns.

GIQL needs to know how your genomic data is structured in order to translate
genomic operators into SQL. This is done through ``Table`` objects, which
map a logical "genomic column" (used in your queries) to the physical columns
in your files, data frames, or database tables.

In GIQL queries, you use a logical genomic column name like ``interval``:

.. code-block:: sql

   SELECT * FROM variants WHERE interval INTERSECTS 'chr1:1000-2000'

Behind the scenes, GIQL expands this to actual column comparisons:

.. code-block:: sql

   SELECT * FROM variants
   WHERE "chrom" = 'chr1' AND "start" < 2000 AND "end" > 1000

The ``Table`` configuration tells GIQL which physical columns (``chrom``,
``start``, ``end``) correspond to the logical ``interval`` column.

Configuring Tables
------------------

Basic Configuration
~~~~~~~~~~~~~~~~~~~

For tables that use the default column names (``chrom``, ``start``, ``end``,
``strand``), pass the table name as a string:

.. code-block:: python

   from giql import transpile

   sql = transpile(
       """
       SELECT * FROM variants
       WHERE interval INTERSECTS 'chr1:1000-2000'
       """,
       tables=["variants"],
   )

Default Columns
~~~~~~~~~~~~~~~

GIQL uses these default column names:

- **chrom**: The chromosome/contig identifier (e.g., 'chr1', 'chrX')
- **start**: The start position of the genomic interval (0-based, inclusive)
- **end**: The end position of the genomic interval (0-based, exclusive)
- **strand**: Strand orientation ('+', '-', or '.'), optional

The default genomic pseudo-column name is ``interval``.

.. _custom-column-names:

Custom Column Names
~~~~~~~~~~~~~~~~~~~

If your table uses different column names, create a ``Table`` object with
the mapping:

.. code-block:: python

   from giql import Table, transpile

   sql = transpile(
       """
       SELECT * FROM my_table
       WHERE interval INTERSECTS 'chr1:1000-2000'
       """,
       tables=[
           Table(
               "my_table",
               chrom_col="chrom",          # Your chromosome column
               start_col="chromStart",     # Your start column (UCSC-style)
               end_col="chromEnd",         # Your end column
           )
       ],
   )

Configuring Multiple Tables
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Pass all tables that participate in genomic queries:

.. code-block:: python

   from giql import transpile

   # Tables with default column names
   sql = transpile(
       """
       SELECT v.*, g.name AS gene_name
       FROM variants v
       JOIN genes g ON v.interval INTERSECTS g.interval
       """,
       tables=["variants", "genes"],
   )

Different Schemas Per Table
~~~~~~~~~~~~~~~~~~~~~~~~~~~

Tables can have different column names and even different genomic column names.
Mix strings (for default columns) with ``Table`` objects (for custom columns):

.. code-block:: python

   from giql import Table, transpile

   sql = transpile(
       """
       SELECT v.ID, g.gene_name
       FROM variants v
       JOIN genes g ON v.var_interval INTERSECTS g.gene_interval
       """,
       tables=[
           # VCF-style columns
           Table(
               "variants",
               genomic_col="var_interval",
               chrom_col="CHROM",
               start_col="POS",
               end_col="END",
               strand_col=None,
           ),
           # BED-style columns (defaults)
           Table(
               "genes",
               genomic_col="gene_interval",
           ),
       ],
   )

Coordinate Systems
------------------

Understanding BED Coordinates
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

GIQL uses the BED coordinate convention by default:

- **0-based start**: The first base of a chromosome is position 0
- **Half-open intervals**: Start is inclusive, end is exclusive
- **Interval [start, end)**: Contains positions from start to end-1

Example: An interval ``chr1:100-200`` covers bases 100 through 199 (100 bases total).

Working with 1-Based Coordinates
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

If your data uses 1-based coordinates (like VCF or GFF), configure the
``Table`` accordingly:

.. code-block:: python

   from giql import Table, transpile

   sql = transpile(
       query,
       tables=[
           Table(
               "variants",
               coordinate_system="1based",
               interval_type="closed",
           )
       ],
   )

Working with Point Features
~~~~~~~~~~~~~~~~~~~~~~~~~~~

For point features (like SNPs), create an interval of length 1:

.. code-block:: python

   # For a SNP at position 1000 (1-based)
   # 0-based interval: [999, 1000)
   start = 999
   end = 1000
