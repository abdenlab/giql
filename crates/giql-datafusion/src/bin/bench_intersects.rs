//! Benchmark binary for the INTERSECTS join optimizer.
//!
//! Usage:
//!   bench_intersects <left.parquet> <right.parquet> [--reps N] [--op join|pairs]
//!
//! Outputs one JSON line per timed run:
//!   {"op":"intersect_join","engine":"giql-optimized","rep":0,"time_s":0.123,"n_rows":4567}

use std::path::PathBuf;
use std::time::Instant;

use datafusion::execution::SessionStateBuilder;
use datafusion::prelude::*;
use giql_datafusion::{register_optimizer, IntersectsOptimizerConfig};

const INTERSECT_JOIN_SQL: &str = "\
    SELECT DISTINCT a.chrom, a.start, a.\"end\" \
    FROM a JOIN b \
    ON a.chrom = b.chrom \
    AND a.start < b.\"end\" \
    AND a.\"end\" > b.start";

const INTERSECT_PAIRS_SQL: &str = "\
    SELECT a.chrom, a.start, a.\"end\", \
           b.chrom AS chrom_b, b.start AS start_b, b.\"end\" AS end_b \
    FROM a JOIN b \
    ON a.chrom = b.chrom \
    AND a.start < b.\"end\" \
    AND a.\"end\" > b.start";

#[tokio::main]
async fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args: Vec<String> = std::env::args().collect();

    if args.len() < 3 {
        eprintln!(
            "Usage: {} <left.parquet> <right.parquet> \
             [--reps N] [--op join|pairs]",
            args[0]
        );
        std::process::exit(1);
    }

    let left_path = PathBuf::from(&args[1]);
    let right_path = PathBuf::from(&args[2]);

    let mut reps = 3;
    let mut op = "join".to_string();
    let mut no_optimizer = false;
    let mut force_binned = false;
    let mut sql_binned: Option<usize> = None;

    let mut i = 3;
    while i < args.len() {
        match args[i].as_str() {
            "--reps" => {
                i += 1;
                reps = args[i].parse()?;
            }
            "--op" => {
                i += 1;
                op = args[i].clone();
            }
            "--no-optimizer" => {
                no_optimizer = true;
            }
            "--force-binned" => {
                force_binned = true;
            }
            "--sql-binned" => {
                i += 1;
                sql_binned = Some(args[i].parse()?);
            }
            _ => {
                eprintln!("Unknown arg: {}", args[i]);
                std::process::exit(1);
            }
        }
        i += 1;
    }

    let sql: String = if let Some(bs) = sql_binned {
        // Run the pure SQL binned approach through Rust DF
        format!(
            "WITH __giql_left AS (\
             SELECT *, UNNEST(range(CAST(\"start\" / {bs} AS BIGINT), \
             CAST((\"end\" - 1) / {bs} + 1 AS BIGINT))) AS __giql_bin \
             FROM a), \
             __giql_right AS (\
             SELECT *, UNNEST(range(CAST(\"start\" / {bs} AS BIGINT), \
             CAST((\"end\" - 1) / {bs} + 1 AS BIGINT))) AS __giql_bin \
             FROM b) \
             SELECT DISTINCT \
             l.\"chrom\", l.\"start\", l.\"end\", \
             r.\"chrom\" AS chrom_r, r.\"start\" AS start_r, r.\"end\" AS end_r \
             FROM __giql_left AS l \
             JOIN __giql_right AS r \
             ON l.\"chrom\" = r.\"chrom\" AND l.__giql_bin = r.__giql_bin \
             WHERE l.\"start\" < r.\"end\" AND l.\"end\" > r.\"start\""
        )
    } else {
        match op.as_str() {
            "join" => INTERSECT_JOIN_SQL.to_string(),
            "pairs" => INTERSECT_PAIRS_SQL.to_string(),
            _ => {
                eprintln!("Unknown op: {op}. Use 'join' or 'pairs'.");
                std::process::exit(1);
            }
        }
    };

    let op_name = format!("intersect_{op}");

    let ctx = if no_optimizer {
        SessionContext::new()
    } else {
        let config = if force_binned {
            // Set thresholds so high that sweep-line is never chosen
            IntersectsOptimizerConfig {
                p99_median_threshold: f64::MAX,
                cv_threshold: f64::MAX,
                max_sample_row_groups: 3,
            }
        } else {
            IntersectsOptimizerConfig::default()
        };
        let state = SessionStateBuilder::new()
            .with_default_features()
            .build();
        let state = register_optimizer(state, config);
        SessionContext::from(state)
    };

    ctx.register_parquet("a", left_path.to_str().unwrap(), Default::default())
        .await?;
    ctx.register_parquet("b", right_path.to_str().unwrap(), Default::default())
        .await?;

    // Warmup
    let _ = ctx.sql(&sql).await?.collect().await?;

    // Timed reps
    for rep in 0..reps {
        let t0 = Instant::now();
        let batches = ctx.sql(&sql).await?.collect().await?;
        let elapsed = t0.elapsed().as_secs_f64();
        let n_rows: usize = batches.iter().map(|b| b.num_rows()).sum();

        println!(
            "{{\"op\":\"{op_name}\",\"engine\":\"giql-optimized\",\
             \"rep\":{rep},\"time_s\":{elapsed:.6},\"n_rows\":{n_rows}}}"
        );
    }

    Ok(())
}
