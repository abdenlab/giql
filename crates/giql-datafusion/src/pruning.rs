use crate::stats::IntervalStats;

/// Domain bounds for one side of the join, derived from row group
/// metadata.
#[derive(Debug, Clone)]
pub struct DomainBounds {
    /// Global minimum start position across all row groups.
    pub min_start: i64,
    /// Global maximum end position across all row groups.
    pub max_end: i64,
}

impl From<&IntervalStats> for DomainBounds {
    fn from(stats: &IntervalStats) -> Self {
        Self {
            min_start: stats.domain_min,
            max_end: stats.domain_max,
        }
    }
}

/// Generate pruning predicates from domain bounds.
///
/// These predicates are always true for the join result set but help
/// the database engine skip row groups that are entirely outside the
/// other side's range. They should be injected as additional filter
/// predicates on the `ParquetExec` nodes before the join.
///
/// Returns predicate expressions as `(column_name, op, literal_value)`
/// tuples that can be converted to DataFusion `PhysicalExpr` nodes.
pub fn generate_pruning_predicates(
    left_bounds: &DomainBounds,
    right_bounds: &DomainBounds,
) -> Vec<PruningPredicate> {
    vec![
        // Left rows that start after right's max end cannot overlap
        PruningPredicate {
            side: JoinSide::Left,
            column: PruningColumn::Start,
            op: PruningOp::Lt,
            value: right_bounds.max_end,
        },
        // Left rows that end before right's min start cannot overlap
        PruningPredicate {
            side: JoinSide::Left,
            column: PruningColumn::End,
            op: PruningOp::Gt,
            value: right_bounds.min_start,
        },
        // Right rows that start after left's max end cannot overlap
        PruningPredicate {
            side: JoinSide::Right,
            column: PruningColumn::Start,
            op: PruningOp::Lt,
            value: left_bounds.max_end,
        },
        // Right rows that end before left's min start cannot overlap
        PruningPredicate {
            side: JoinSide::Right,
            column: PruningColumn::End,
            op: PruningOp::Gt,
            value: left_bounds.min_start,
        },
    ]
}

/// A pruning predicate to inject on a ParquetExec node.
#[derive(Debug, Clone)]
pub struct PruningPredicate {
    /// Which side of the join this predicate applies to.
    pub side: JoinSide,
    /// Which interval column to filter.
    pub column: PruningColumn,
    /// Comparison operator.
    pub op: PruningOp,
    /// Literal value to compare against.
    pub value: i64,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum JoinSide {
    Left,
    Right,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PruningColumn {
    Start,
    End,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PruningOp {
    Lt,
    Gt,
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_generate_pruning_predicates() {
        let left = DomainBounds {
            min_start: 1000,
            max_end: 50000,
        };
        let right = DomainBounds {
            min_start: 2000,
            max_end: 60000,
        };

        let preds = generate_pruning_predicates(&left, &right);
        assert_eq!(preds.len(), 4);

        // Left start < 60000 (right max_end)
        assert_eq!(preds[0].side, JoinSide::Left);
        assert_eq!(preds[0].column, PruningColumn::Start);
        assert_eq!(preds[0].op, PruningOp::Lt);
        assert_eq!(preds[0].value, 60000);

        // Left end > 2000 (right min_start)
        assert_eq!(preds[1].side, JoinSide::Left);
        assert_eq!(preds[1].column, PruningColumn::End);
        assert_eq!(preds[1].op, PruningOp::Gt);
        assert_eq!(preds[1].value, 2000);

        // Right start < 50000 (left max_end)
        assert_eq!(preds[2].side, JoinSide::Right);
        assert_eq!(preds[2].column, PruningColumn::Start);
        assert_eq!(preds[2].op, PruningOp::Lt);
        assert_eq!(preds[2].value, 50000);

        // Right end > 1000 (left min_start)
        assert_eq!(preds[3].side, JoinSide::Right);
        assert_eq!(preds[3].column, PruningColumn::End);
        assert_eq!(preds[3].op, PruningOp::Gt);
        assert_eq!(preds[3].value, 1000);
    }
}
