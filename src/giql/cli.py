"""Command-line interface for GIQL.

This module provides a CLI that mirrors bedtools intersect functionality
using GIQL's genomic query capabilities.
"""

import sys
from pathlib import Path

import click
import duckdb
from oxbow import from_bam
from oxbow import from_bed
from oxbow import from_gff
from oxbow import from_gtf
from oxbow import from_vcf

from giql import GIQLEngine


@click.group()
@click.version_option()
def cli():
    """GIQL - Genomic Interval Query Language.

    SQL-based toolkit for genomic range queries.
    """
    pass


def _detect_file_format(file_path: Path) -> str:
    """Detect genomic file format from file extension.

    :param file_path:
        Path to the file
    :return:
        Format identifier: 'bed', 'bam', 'vcf', 'gff', 'gtf'
    :raises click.ClickException:
        If format cannot be determined
    """
    # Handle compressed files
    suffixes = file_path.suffixes
    if suffixes[-1] == ".gz":
        # Remove .gz and check the actual format
        ext = suffixes[-2] if len(suffixes) >= 2 else ""
    else:
        ext = file_path.suffix

    ext = ext.lower()

    format_map = {
        ".bed": "bed",
        ".bam": "bam",
        ".vcf": "vcf",
        ".gff": "gff",
        ".gff3": "gff",
        ".gtf": "gtf",
    }

    if ext in format_map:
        return format_map[ext]

    raise click.ClickException(
        f"Unsupported file format: {ext}. Supported formats: BED, BAM, VCF, GFF, GTF"
    )


def _load_genomic_file(
    conn: duckdb.DuckDBPyConnection, file_path: Path, table_name: str
) -> dict[str, str]:
    """Load genomic file using appropriate oxbow function.

    :param conn:
        DuckDB connection
    :param file_path:
        Path to genomic file
    :param table_name:
        Name for the table to create
    :return:
        Dictionary mapping column names to types
    :raises click.ClickException:
        If file cannot be loaded
    """
    fmt = _detect_file_format(file_path)
    compression = "gzip" if file_path.suffix == ".gz" else None

    try:
        match fmt:
            case "bed":
                df = from_bed(str(file_path), compression=compression).to_duckdb(conn)
            case "bam":
                df = from_bam(str(file_path)).to_duckdb(conn)
            case "vcf":
                df = from_vcf(str(file_path), compression=compression).to_duckdb(conn)
            case "gff":
                df = from_gff(str(file_path), compression=compression).to_duckdb(conn)
            case "gtf":
                df = from_gtf(str(file_path), compression=compression).to_duckdb(conn)
            case _:
                raise click.ClickException(f"Unsupported format: {fmt}")

        conn.execute(f"CREATE TABLE {table_name} AS SELECT * FROM df")

        # Get column information
        col_info = conn.execute(f"DESCRIBE {table_name}").fetchall()
        return {col[0]: col[1] for col in col_info}

    except Exception as e:
        raise click.ClickException(f"Failed to load {file_path}: {e}")


def _expand_rest_columns(df):
    """Expand 'rest' columns from BED files into separate columns.

    BED files store extra fields beyond chrom/start/end in a 'rest' column
    as a tab-delimited string. This function expands those into separate columns
    to match bedtools output format.

    :param df:
        DataFrame with potential 'rest' columns
    :return:
        DataFrame with rest columns expanded
    """
    import pandas as pd

    # pandas.read_sql can return duplicate column names when joining
    # Find all 'rest' column positions
    rest_indices = [i for i, col in enumerate(df.columns) if col == "rest"]

    if not rest_indices:
        return df

    # Build new dataframe with expanded columns
    # We need to handle duplicate column names, so we can't use a dict
    new_data = {}
    new_col_names = []

    for i, col in enumerate(df.columns):
        if col == "rest" and i in rest_indices:
            # Expand this rest column
            col_data = df.iloc[:, i]
            expanded = col_data.fillna("").astype(str).str.split("\t", expand=True)

            # Add all expanded columns with unique names
            for j in range(expanded.shape[1]):
                col_name = f"field_{j + 4}"
                # Make unique if duplicate
                base_name = col_name
                counter = 0
                while col_name in new_col_names:
                    counter += 1
                    col_name = f"{base_name}_{counter}"
                new_col_names.append(col_name)
                new_data[col_name] = expanded[j]
        else:
            # Keep non-rest columns as-is
            # Make unique names for duplicates
            col_name = col
            base_name = col_name
            counter = 0
            while col_name in new_col_names:
                counter += 1
                col_name = f"{base_name}_{counter}"
            new_col_names.append(col_name)
            new_data[col_name] = df.iloc[:, i]

    # Rebuild dataframe with explicit column order
    result = pd.DataFrame(new_data, columns=new_col_names)
    return result


def _detect_genomic_columns(columns: dict[str, str]) -> dict[str, str | None]:
    """Detect genomic coordinate columns from available columns.

    :param columns:
        Dictionary of column name -> type
    :return:
        Dictionary with keys: chrom_col, start_col, end_col, strand_col
    """
    col_names = {c.lower(): c for c in columns.keys()}

    # Chromosome column patterns (in priority order)
    chrom_col = None
    for pattern in ["chrom", "seqid", "chr", "chromosome", "contig", "seqname"]:
        if pattern in col_names:
            chrom_col = col_names[pattern]
            break

    # Start column patterns
    start_col = None
    for pattern in [
        "start",
        "chromstart",
        "pos",
        "begin",
        "txstart",
        "cdsstart",
        "thickstart",
    ]:
        if pattern in col_names:
            start_col = col_names[pattern]
            break

    # End column patterns
    end_col = None
    for pattern in [
        "end",
        "chromend",
        "stop",
        "txend",
        "cdsend",
        "thickend",
    ]:
        if pattern in col_names:
            end_col = col_names[pattern]
            break

    # Strand column patterns
    strand_col = None
    for pattern in ["strand", "str", "orientation"]:
        if pattern in col_names:
            strand_col = col_names[pattern]
            break

    return {
        "chrom_col": chrom_col,
        "start_col": start_col,
        "end_col": end_col,
        "strand_col": strand_col,
    }


@cli.command()
@click.option(
    "-a",
    "--file-a",
    required=True,
    type=click.Path(exists=True),
    help="BAM/BED/GFF/VCF file 'A'. Each feature in A is compared to B.",
)
@click.option(
    "-b",
    "--file-b",
    required=True,
    multiple=True,
    type=click.Path(exists=True),
    help="One or more BAM/BED/GFF/VCF files for comparison.",
)
@click.option(
    "-wa",
    "--write-a",
    is_flag=True,
    help="Write the original entry in A for each overlap.",
)
@click.option(
    "-wb",
    "--write-b",
    is_flag=True,
    help="Write the original entry in B for each overlap.",
)
@click.option(
    "-loj",
    "--left-outer-join",
    is_flag=True,
    help="Perform left outer join. Report all A features with NULL B when no overlap.",
)
@click.option(
    "-wo",
    "--write-overlap",
    is_flag=True,
    help="Write the number of overlapping base pairs between features.",
)
@click.option(
    "-wao",
    "--write-all-overlap",
    is_flag=True,
    help="Like -wo but includes A features with zero overlap.",
)
@click.option(
    "-u",
    "--unique",
    is_flag=True,
    help="Report each A feature only once if any overlap exists in B.",
)
@click.option(
    "-c",
    "--count",
    is_flag=True,
    help="For each entry in A, report the number of overlaps in B.",
)
@click.option(
    "-v",
    "--invert",
    is_flag=True,
    help="Only report entries in A that have no overlap in B.",
)
@click.option(
    "-f",
    "--fraction-a",
    type=float,
    help="Minimum overlap as fraction of A.",
)
@click.option(
    "-F",
    "--fraction-b",
    type=float,
    help="Minimum overlap as fraction of B.",
)
@click.option(
    "-r",
    "--reciprocal",
    is_flag=True,
    help="Require reciprocal overlap fraction for both A and B.",
)
@click.option(
    "-e",
    "--either",
    is_flag=True,
    help="Require that -f OR -F be satisfied (not both).",
)
@click.option(
    "-s",
    "--same-strand",
    is_flag=True,
    help="Require same strand for overlaps.",
)
@click.option(
    "-S",
    "--opposite-strand",
    is_flag=True,
    help="Require opposite strand for overlaps.",
)
@click.option(
    "--header",
    is_flag=True,
    help="Print the header from A before results.",
)
@click.option(
    "--names",
    multiple=True,
    help="Aliases for B files (instead of file numbers).",
)
@click.option(
    "-sorted",
    "--sorted-input",
    is_flag=True,
    help="For compatibility with bedtools (currently ignored).",
)
def intersect(
    file_a,
    file_b,
    write_a,
    write_b,
    left_outer_join,
    write_overlap,
    write_all_overlap,
    unique,
    count,
    invert,
    fraction_a,
    fraction_b,
    reciprocal,
    either,
    same_strand,
    opposite_strand,
    header,
    names,
    sorted_input,
):
    """Find overlaps between genomic features.

    Similar to bedtools intersect, this command finds overlapping intervals
    between files A and B using GIQL's spatial operators.

    Supports BED, BAM, VCF, GFF, and GTF formats (gzip compressed or uncompressed).
    """
    # Validate conflicting options
    if same_strand and opposite_strand:
        raise click.UsageError("Cannot use -s and -S together")

    output_modes = [
        write_a,
        write_b,
        left_outer_join,
        write_overlap,
        write_all_overlap,
        unique,
        count,
        invert,
    ]
    if sum(output_modes) > 1:
        raise click.UsageError("Can only specify one output mode")

    # Create DuckDB connection
    conn = duckdb.connect()

    # Initialize engine with existing connection
    engine = GIQLEngine(target_dialect="duckdb", connection=conn)

    try:
        # Load file A
        file_a_path = Path(file_a)
        table_a = "file_a"
        columns_a = _load_genomic_file(conn, file_a_path, table_a)

        # Detect genomic columns
        genomic_cols_a = _detect_genomic_columns(columns_a)

        if not all(
            [
                genomic_cols_a["chrom_col"],
                genomic_cols_a["start_col"],
                genomic_cols_a["end_col"],
            ]
        ):
            raise click.ClickException(
                f"Could not detect genomic columns in {file_a}. "
                f"Found columns: {list(columns_a.keys())}"
            )

        # Register schema for file A
        engine.register_table_schema(
            table_a,
            columns_a,
            genomic_column="position",
            chrom_col=genomic_cols_a["chrom_col"],
            start_col=genomic_cols_a["start_col"],
            end_col=genomic_cols_a["end_col"],
            strand_col=genomic_cols_a["strand_col"],
        )

        # Process file(s) B
        results = []
        for idx, b_file in enumerate(file_b):
            b_path = Path(b_file)
            table_b = f"file_b_{idx}"

            # Load file B
            columns_b = _load_genomic_file(conn, b_path, table_b)

            # Detect genomic columns in B
            genomic_cols_b = _detect_genomic_columns(columns_b)

            if not all(
                [
                    genomic_cols_b["chrom_col"],
                    genomic_cols_b["start_col"],
                    genomic_cols_b["end_col"],
                ]
            ):
                raise click.ClickException(
                    f"Could not detect genomic columns in {b_file}"
                )

            # Register schema for file B
            engine.register_table_schema(
                table_b,
                columns_b,
                genomic_column="region",
                chrom_col=genomic_cols_b["chrom_col"],
                start_col=genomic_cols_b["start_col"],
                end_col=genomic_cols_b["end_col"],
                strand_col=genomic_cols_b["strand_col"],
            )

            # Build query based on options
            query = _build_intersect_query(
                table_a=table_a,
                table_b=table_b,
                chrom_a=genomic_cols_a["chrom_col"],
                start_a=genomic_cols_a["start_col"],
                end_a=genomic_cols_a["end_col"],
                strand_a=genomic_cols_a["strand_col"],
                chrom_b=genomic_cols_b["chrom_col"],
                start_b=genomic_cols_b["start_col"],
                end_b=genomic_cols_b["end_col"],
                strand_b=genomic_cols_b["strand_col"],
                write_a=write_a,
                write_b=write_b,
                left_outer_join=left_outer_join,
                write_overlap=write_overlap,
                write_all_overlap=write_all_overlap,
                unique=unique,
                count=count,
                invert=invert,
                same_strand=same_strand,
                opposite_strand=opposite_strand,
                fraction_a=fraction_a,
                fraction_b=fraction_b,
                reciprocal=reciprocal,
                either=either,
            )

            result = engine.query(query)

            # Expand "rest" columns if present (BED format compatibility)
            result = _expand_rest_columns(result)

            # If names provided, add a column indicating which B file
            if names and idx < len(names):
                result["b_file"] = names[idx]
            elif len(file_b) > 1:
                result["b_file"] = b_path.name

            results.append(result)

        # Combine results if multiple B files
        if results:
            if len(results) > 1:
                import pandas as pd

                final_result = pd.concat(results, ignore_index=True)
            else:
                final_result = results[0]

            # Output results
            if header:
                # Print header
                print("\t".join(final_result.columns))

            # Print results as TSV
            final_result.to_csv(sys.stdout, sep="\t", index=False, header=False)

    finally:
        engine.close()


def _build_intersect_query(
    table_a: str,
    table_b: str,
    chrom_a: str,
    start_a: str,
    end_a: str,
    strand_a: str | None,
    chrom_b: str,
    start_b: str,
    end_b: str,
    strand_b: str | None,
    write_a: bool = False,
    write_b: bool = False,
    left_outer_join: bool = False,
    write_overlap: bool = False,
    write_all_overlap: bool = False,
    unique: bool = False,
    count: bool = False,
    invert: bool = False,
    same_strand: bool = False,
    opposite_strand: bool = False,
    fraction_a: float | None = None,
    fraction_b: float | None = None,
    reciprocal: bool = False,
    either: bool = False,
) -> str:
    """Build GIQL query based on intersect options."""

    # Build strand filter if needed
    strand_filter = ""
    if same_strand and strand_a and strand_b:
        strand_filter = f' AND a."{strand_a}" = b."{strand_b}"'
    elif opposite_strand and strand_a and strand_b:
        strand_filter = f' AND a."{strand_a}" != b."{strand_b}"'

    # Build fraction filter if needed
    fraction_filter = ""
    if fraction_a or fraction_b:
        filters = []

        if fraction_a:
            # Overlap must be at least fraction_a of A's length
            overlap_expr = (
                f'LEAST(a."{end_a}", b."{end_b}") - '
                f'GREATEST(a."{start_a}", b."{start_b}")'
            )
            a_length = f'(a."{end_a}" - a."{start_a}")'
            filters.append(f"({overlap_expr}::FLOAT / {a_length} >= {fraction_a})")

        if fraction_b:
            # Overlap must be at least fraction_b of B's length
            overlap_expr = (
                f'LEAST(a."{end_a}", b."{end_b}") - '
                f'GREATEST(a."{start_a}", b."{start_b}")'
            )
            b_length = f'(b."{end_b}" - b."{start_b}")'
            filters.append(f"({overlap_expr}::FLOAT / {b_length} >= {fraction_b})")

        # Combine filters based on reciprocal/either flags
        if reciprocal and len(filters) == 2:
            # Both must be satisfied (AND)
            fraction_filter = f" AND ({filters[0]} AND {filters[1]})"
        elif either and len(filters) == 2:
            # Either must be satisfied (OR)
            fraction_filter = f" AND ({filters[0]} OR {filters[1]})"
        elif filters:
            # Just one filter or default behavior
            fraction_filter = f" AND {' AND '.join(filters)}"

    if invert:
        # Only features in A with no overlap in B
        where_clause = f"a.position INTERSECTS b.region{strand_filter}{fraction_filter}"
        return f"""
            SELECT a.*
            FROM {table_a} a
            WHERE NOT EXISTS (
                SELECT 1 FROM {table_b} b
                WHERE {where_clause}
            )
        """

    if count:
        # Count overlaps
        # Get all columns from table A for GROUP BY
        on_clause = f"a.position INTERSECTS b.region{strand_filter}{fraction_filter}"
        return f"""
            SELECT a.*, COUNT(b.\"{chrom_b}\") as overlap_count
            FROM {table_a} a
            LEFT JOIN {table_b} b ON {on_clause}
            GROUP BY ALL
        """

    if unique:
        # Report each A feature only once if overlaps exist
        on_clause = f"a.position INTERSECTS b.region{strand_filter}{fraction_filter}"
        return f"""
            SELECT DISTINCT a.*
            FROM {table_a} a
            JOIN {table_b} b ON {on_clause}
        """

    if left_outer_join or write_all_overlap:
        # Left outer join
        join_type = "LEFT JOIN"
    else:
        join_type = "JOIN"

    # Build select clause
    if write_a and not write_b:
        select_clause = "a.*"
    elif write_b and not write_a:
        select_clause = "b.*"
    else:
        # Default: write both A and B
        select_clause = "a.*, b.*"

    # Add overlap calculation if requested
    if write_overlap or write_all_overlap:
        # Calculate overlap size: min(end_a, end_b) - max(start_a, start_b)
        overlap_expr = f"""
            CASE
                WHEN b.\"{chrom_b}\" IS NULL THEN 0
                ELSE GREATEST(0,
                    LEAST(a.\"{end_a}\", b.\"{end_b}\") -
                    GREATEST(a.\"{start_a}\", b.\"{start_b}\")
                )
            END as overlap_bp
        """
        select_clause = f"{select_clause}, {overlap_expr}"

    # Build ON clause
    on_clause = f"a.position INTERSECTS b.region{strand_filter}{fraction_filter}"

    # Build base query
    query = f"""
        SELECT {select_clause}
        FROM {table_a} a
        {join_type} {table_b} b ON {on_clause}
    """

    return query


if __name__ == "__main__":
    cli()
