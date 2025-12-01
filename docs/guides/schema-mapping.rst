Schema Mapping Guide
====================

This guide explains how to configure GIQL to work with your genomic data by
registering table schemas and mapping logical genomic columns.

.. contents::
   :local:
   :depth: 2

Understanding Schema Mapping
----------------------------

GIQL needs to know how your genomic data is structured in order to translate
genomic operators into SQL. This is done through schema registration, which
maps a logical "genomic column" (used in your queries) to the physical columns
in your database tables.

The Core Concept
~~~~~~~~~~~~~~~~

In GIQL queries, you use a logical genomic column name like ``interval``:

.. code-block:: sql

   SELECT * FROM variants WHERE interval INTERSECTS 'chr1:1000-2000'

Behind the scenes, GIQL expands this to actual column comparisons:

.. code-block:: sql

   SELECT * FROM variants
   WHERE chromosome = 'chr1' AND start_pos < 2000 AND end_pos > 1000

Schema registration tells GIQL which physical columns (``chromosome``,
``start_pos``, ``end_pos``) correspond to the logical ``interval`` column.

Registering Table Schemas
-------------------------

Basic Registration
~~~~~~~~~~~~~~~~~~

Register a table schema using ``register_table_schema()``:

.. code-block:: python

   from giql import GIQLEngine

   with GIQLEngine(target_dialect="duckdb") as engine:
       # Load data
       engine.load_csv("variants", "variants.csv")

       # Register schema
       engine.register_table_schema(
           "variants",                    # Table name
           {
               "id": "INTEGER",
               "chromosome": "VARCHAR",
               "start_pos": "BIGINT",
               "end_pos": "BIGINT",
               "name": "VARCHAR",
               "quality": "FLOAT",
           },
           genomic_column="interval",     # Logical column name for queries
       )

       # Now you can use 'interval' in queries
       cursor = engine.execute("""
           SELECT * FROM variants
           WHERE interval INTERSECTS 'chr1:1000-2000'
       """)

Required Columns
~~~~~~~~~~~~~~~~

For schema registration, your table must have columns that map to:

- **chromosome**: The chromosome/contig identifier (e.g., 'chr1', 'chrX')
- **start_pos**: The start position of the genomic interval (0-based, inclusive)
- **end_pos**: The end position of the genomic interval (0-based, exclusive)

GIQL looks for these column names by default. If your columns have different
names, see :ref:`custom-column-names`.

Optional Strand Column
~~~~~~~~~~~~~~~~~~~~~~

If your data includes strand information, include it in the schema:

.. code-block:: python

   engine.register_table_schema(
       "features",
       {
           "chromosome": "VARCHAR",
           "start_pos": "BIGINT",
           "end_pos": "BIGINT",
           "strand": "VARCHAR",       # '+', '-', or '.'
           "name": "VARCHAR",
       },
       genomic_column="interval",
   )

The strand column enables strand-specific operations in operators like
CLUSTER and NEAREST.

.. _custom-column-names:

Custom Column Names
~~~~~~~~~~~~~~~~~~~

If your table uses different column names for genomic coordinates, specify
the mapping explicitly:

.. code-block:: python

   engine.register_table_schema(
       "my_table",
       {
           "chrom": "VARCHAR",        # Your chromosome column
           "chromStart": "BIGINT",    # Your start column (UCSC-style)
           "chromEnd": "BIGINT",      # Your end column
           "name": "VARCHAR",
       },
       genomic_column="interval",
       chromosome_column="chrom",      # Map to your column name
       start_column="chromStart",      # Map to your column name
       end_column="chromEnd",          # Map to your column name
   )

Multiple Tables
---------------

Register Multiple Tables
~~~~~~~~~~~~~~~~~~~~~~~~

Register all tables that will participate in genomic queries:

.. code-block:: python

   with GIQLEngine(target_dialect="duckdb") as engine:
       # Load data files
       engine.load_csv("variants", "variants.bed")
       engine.load_csv("genes", "genes.bed")
       engine.load_csv("regulatory", "regulatory.bed")

       # Define common schema
       bed_schema = {
           "chromosome": "VARCHAR",
           "start_pos": "BIGINT",
           "end_pos": "BIGINT",
           "name": "VARCHAR",
           "score": "FLOAT",
           "strand": "VARCHAR",
       }

       # Register each table
       for table in ["variants", "genes", "regulatory"]:
           engine.register_table_schema(
               table,
               bed_schema,
               genomic_column="interval",
           )

       # Now you can join tables using genomic operators
       cursor = engine.execute("""
           SELECT v.*, g.name AS gene_name
           FROM variants v
           JOIN genes g ON v.interval INTERSECTS g.interval
       """)

Different Schemas Per Table
~~~~~~~~~~~~~~~~~~~~~~~~~~~

Tables can have different schemas and even different genomic column names:

.. code-block:: python

   # Variants table with VCF-style columns
   engine.register_table_schema(
       "variants",
       {
           "CHROM": "VARCHAR",
           "POS": "BIGINT",
           "END": "BIGINT",
           "ID": "VARCHAR",
           "QUAL": "FLOAT",
       },
       genomic_column="var_interval",
       chromosome_column="CHROM",
       start_column="POS",
       end_column="END",
   )

   # Genes table with BED-style columns
   engine.register_table_schema(
       "genes",
       {
           "chromosome": "VARCHAR",
           "start_pos": "BIGINT",
           "end_pos": "BIGINT",
           "gene_name": "VARCHAR",
           "strand": "VARCHAR",
       },
       genomic_column="gene_interval",
   )

   # Query using different genomic column names
   cursor = engine.execute("""
       SELECT v.ID, g.gene_name
       FROM variants v
       JOIN genes g ON v.var_interval INTERSECTS g.gene_interval
   """)

Coordinate Systems
------------------

Understanding BED Coordinates
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

GIQL uses the BED coordinate convention:

- **0-based start**: The first base of a chromosome is position 0
- **Half-open intervals**: Start is inclusive, end is exclusive
- **Interval [start, end)**: Contains positions from start to end-1

Example: An interval ``chr1:100-200`` covers bases 100 through 199 (100 bases total).

Converting from 1-Based Coordinates
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

If your data uses 1-based coordinates (like VCF or GFF), convert when loading:

.. code-block:: python

   import pandas as pd

   # Load 1-based data
   df = pd.read_csv("variants.vcf", sep="\t")

   # Convert to 0-based
   df['start_pos'] = df['POS'] - 1  # Convert 1-based to 0-based
   df['end_pos'] = df['POS']        # For SNPs, end = start + 1

   # Load into engine
   engine.conn.execute("CREATE TABLE variants AS SELECT * FROM df")

   # Register schema
   engine.register_table_schema(
       "variants",
       {
           "chromosome": "VARCHAR",
           "start_pos": "BIGINT",
           "end_pos": "BIGINT",
           # ... other columns
       },
       genomic_column="interval",
   )

Working with Point Features
~~~~~~~~~~~~~~~~~~~~~~~~~~~

For point features (like SNPs), create an interval of length 1:

.. code-block:: python

   # For a SNP at position 1000 (1-based)
   # 0-based interval: [999, 1000)
   start_pos = 999
   end_pos = 1000

Data Types
----------

Recommended Column Types
~~~~~~~~~~~~~~~~~~~~~~~~

For optimal performance, use appropriate data types:

.. list-table::
   :header-rows: 1
   :widths: 25 25 50

   * - Column
     - Recommended Type
     - Notes
   * - chromosome
     - VARCHAR
     - String type for chromosome names
   * - start_pos
     - BIGINT
     - 64-bit integer for large genomes
   * - end_pos
     - BIGINT
     - 64-bit integer for large genomes
   * - strand
     - VARCHAR(1) or CHAR(1)
     - Single character: '+', '-', '.'
   * - score
     - FLOAT or DOUBLE
     - Numeric scores
   * - name
     - VARCHAR
     - Feature identifiers

Type Compatibility
~~~~~~~~~~~~~~~~~~

GIQL schemas use SQL type names. Common mappings:

.. list-table::
   :header-rows: 1
   :widths: 30 35 35

   * - GIQL Schema Type
     - DuckDB Type
     - SQLite Type
   * - INTEGER
     - INTEGER
     - INTEGER
   * - BIGINT
     - BIGINT
     - INTEGER
   * - VARCHAR
     - VARCHAR
     - TEXT
   * - FLOAT
     - FLOAT
     - REAL
   * - DOUBLE
     - DOUBLE
     - REAL

Loading Data
------------

From CSV Files
~~~~~~~~~~~~~~

Load CSV files directly:

.. code-block:: python

   engine.load_csv("features", "features.csv")

   # With custom options
   engine.load_csv(
       "features",
       "features.tsv",
       delimiter="\t",
       header=True,
   )

From Pandas DataFrames
~~~~~~~~~~~~~~~~~~~~~~

Load data from pandas:

.. code-block:: python

   import pandas as pd

   df = pd.read_csv("features.bed", sep="\t", header=None,
                    names=["chromosome", "start_pos", "end_pos", "name"])

   # Register the DataFrame as a table
   engine.conn.execute("CREATE TABLE features AS SELECT * FROM df")

   # Then register the schema
   engine.register_table_schema(
       "features",
       {
           "chromosome": "VARCHAR",
           "start_pos": "BIGINT",
           "end_pos": "BIGINT",
           "name": "VARCHAR",
       },
       genomic_column="interval",
   )

From Existing Database Tables
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

If tables already exist in your database, just register their schemas:

.. code-block:: python

   # Connect to existing database
   with GIQLEngine(target_dialect="duckdb", db_path="my_database.duckdb") as engine:
       # Register schemas for existing tables
       engine.register_table_schema(
           "existing_table",
           {
               "chromosome": "VARCHAR",
               "start_pos": "BIGINT",
               "end_pos": "BIGINT",
               "name": "VARCHAR",
           },
           genomic_column="interval",
       )

       # Query existing data
       cursor = engine.execute("""
           SELECT * FROM existing_table
           WHERE interval INTERSECTS 'chr1:1000-2000'
       """)

Troubleshooting
---------------

Common Issues
~~~~~~~~~~~~~

**"Unknown column" errors:**

- Ensure the table schema is registered before querying
- Check that the genomic column name in your query matches the registered name
- Verify column names in the schema match actual table columns

**Incorrect results:**

- Verify your coordinate system (0-based vs 1-based)
- Check that start_pos < end_pos for all intervals
- Ensure chromosome names match between tables (e.g., 'chr1' vs '1')

**Performance issues:**

- See the :doc:`performance` guide for optimization tips
- Consider adding indexes on genomic columns

Verifying Schema Registration
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Check that schemas are registered correctly:

.. code-block:: python

   # After registration, test with a simple query
   sql = engine.transpile("""
       SELECT * FROM variants
       WHERE interval INTERSECTS 'chr1:1000-2000'
   """)
   print(sql)
   # Should show expanded SQL with chromosome, start_pos, end_pos comparisons
