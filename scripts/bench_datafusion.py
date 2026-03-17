"""Benchmark DataFusion vs DuckDB on GIQL queries over a Parquet file.

Usage:
    python scripts/bench_datafusion.py <parquet_file>

The Parquet file should contain genomic interval data with columns:
chrom, start, end, strand.
"""

import argparse
import time

QUERIES = {
    "INTERSECTS filter": (
        "SELECT * FROM data WHERE interval INTERSECTS 'chr1:1000000-2000000'"
    ),
    "MERGE": "SELECT * FROM MERGE(data)",
    "CLUSTER": "SELECT *, CLUSTER(interval) AS cid FROM data",
    "DISTANCE self-join": (
        "SELECT DISTANCE(a.interval, b.interval) AS dist"
        " FROM data AS a CROSS JOIN data AS b"
        " WHERE a.chrom = 'chr1' AND b.chrom = 'chr1'"
        " LIMIT 1000"
    ),
    "Standalone NEAREST": (
        "SELECT * FROM data"
        " ORDER BY NEAREST(interval, 'chr1:1500000-1500001')"
        " LIMIT 5"
    ),
}


def bench_datafusion(parquet_path: str) -> dict[str, float]:
    from giql.engines.datafusion import DataFusionEngine

    engine = DataFusionEngine()
    engine.register_parquet("data", parquet_path)
    timings = {}
    for name, query in QUERIES.items():
        t0 = time.perf_counter()
        try:
            engine.query(query)
        except Exception as e:
            timings[name] = float("nan")
            print(f"  [DataFusion] {name}: ERROR - {e}")
            continue
        elapsed = time.perf_counter() - t0
        timings[name] = elapsed
    return timings


def bench_duckdb(parquet_path: str) -> dict[str, float]:
    import duckdb

    from giql import Table
    from giql import transpile

    conn = duckdb.connect()
    conn.execute(
        f"CREATE TABLE data AS SELECT * FROM read_parquet('{parquet_path}')"
    )
    tables = [Table("data")]
    timings = {}
    for name, query in QUERIES.items():
        try:
            sql = transpile(query, tables=tables)
        except Exception as e:
            timings[name] = float("nan")
            print(f"  [DuckDB] {name}: TRANSPILE ERROR - {e}")
            continue
        t0 = time.perf_counter()
        try:
            conn.execute(sql).fetchall()
        except Exception as e:
            timings[name] = float("nan")
            print(f"  [DuckDB] {name}: ERROR - {e}")
            continue
        elapsed = time.perf_counter() - t0
        timings[name] = elapsed
    conn.close()
    return timings


def main():
    parser = argparse.ArgumentParser(
        description="Benchmark DataFusion vs DuckDB on GIQL queries"
    )
    parser.add_argument("parquet_file", help="Path to a genomic Parquet file")
    args = parser.parse_args()

    print(f"Benchmarking with: {args.parquet_file}\n")

    print("--- DataFusion ---")
    df_timings = bench_datafusion(args.parquet_file)
    for name, elapsed in df_timings.items():
        print(f"  {name}: {elapsed:.4f}s")

    print("\n--- DuckDB ---")
    duck_timings = bench_duckdb(args.parquet_file)
    for name, elapsed in duck_timings.items():
        print(f"  {name}: {elapsed:.4f}s")

    print("\n--- Comparison (DataFusion / DuckDB) ---")
    for name in QUERIES:
        df_t = df_timings.get(name, float("nan"))
        dk_t = duck_timings.get(name, float("nan"))
        if df_t != df_t or dk_t != dk_t:  # NaN check
            print(f"  {name}: skipped (error)")
        else:
            ratio = df_t / dk_t if dk_t > 0 else float("inf")
            print(f"  {name}: {ratio:.2f}x")


if __name__ == "__main__":
    main()
