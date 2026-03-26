"""Query transformers for GIQL operations.

This module contains transformers that rewrite queries containing GIQL-specific
operations (like CLUSTER and MERGE) into equivalent SQL with CTEs.
"""

from sqlglot import exp

from giql.constants import DEFAULT_CHROM_COL
from giql.constants import DEFAULT_END_COL
from giql.constants import DEFAULT_START_COL
from giql.constants import DEFAULT_STRAND_COL
from giql.expressions import GIQLCluster
from giql.expressions import GIQLCoverage
from giql.expressions import GIQLMerge
from giql.table import Tables

# Mapping from COVERAGE stat parameter to SQL aggregate function
COVERAGE_STAT_MAP = {
    "count": "COUNT",
    "mean": "AVG",
    "sum": "SUM",
    "min": "MIN",
    "max": "MAX",
}


class ClusterTransformer:
    """Transforms queries containing CLUSTER into CTE-based queries.

    CLUSTER cannot be a simple window function because it requires nested
    window functions (LAG inside SUM). Instead, we transform:

        SELECT *, CLUSTER(interval) AS cluster_id FROM features

    Into:

        WITH lag_calc AS (
            SELECT *, LAG(end_pos) OVER (...) AS prev_end FROM features
        )
        SELECT *, SUM(CASE WHEN prev_end >= start_pos ...) AS cluster_id
        FROM lag_calc
    """

    def __init__(self, tables: Tables):
        """Initialize transformer.

        :param tables:
            Table configurations for column mapping
        """
        self.tables = tables

    def _get_table_name(self, query: exp.Select) -> str | None:
        """Extract table name from query's FROM clause.

        :param query:
            Query to extract table name from
        :return:
            Table name if FROM contains a simple table, None otherwise
        """
        from_clause = query.args.get("from_")
        if not from_clause:
            return None

        if isinstance(from_clause.this, exp.Table):
            return from_clause.this.name

        return None

    def _get_genomic_columns(self, query: exp.Select) -> tuple[str, str, str, str]:
        """Get genomic column names from table config or defaults.

        :param query:
            Query to extract table and column info from
        :return:
            Tuple of (chrom_col, start_col, end_col, strand_col)
        """
        table_name = self._get_table_name(query)

        # Default column names
        chrom_col = DEFAULT_CHROM_COL
        start_col = DEFAULT_START_COL
        end_col = DEFAULT_END_COL
        strand_col = DEFAULT_STRAND_COL

        if table_name:
            table = self.tables.get(table_name)
            if table:
                chrom_col = table.chrom_col
                start_col = table.start_col
                end_col = table.end_col
                if table.strand_col:
                    strand_col = table.strand_col

        return chrom_col, start_col, end_col, strand_col

    def transform(self, query: exp.Expression) -> exp.Expression:
        """Transform query if it contains CLUSTER expressions.

        :param query:
            Parsed query AST
        :return:
            Transformed query AST
        """
        if not isinstance(query, exp.Select):
            return query

        # First, recursively transform any CTEs that might contain CLUSTER
        if query.args.get("with_"):
            cte = query.args["with_"]
            for cte_expr in cte.expressions:
                if isinstance(cte_expr, exp.CTE):
                    # Transform the CTE's subquery
                    cte_expr.set("this", self.transform(cte_expr.this))

        # Recursively transform subqueries in FROM clause
        if query.args.get("from_"):
            from_clause = query.args["from_"]
            self._transform_subqueries_in_node(from_clause)

        # Recursively transform subqueries in JOIN clauses
        if query.args.get("joins"):
            for join in query.args["joins"]:
                self._transform_subqueries_in_node(join)

        # Recursively transform subqueries in WHERE clause
        if query.args.get("where"):
            self._transform_subqueries_in_node(query.args["where"])

        # Find all CLUSTER expressions in the SELECT clause
        cluster_exprs = self._find_cluster_expressions(query)

        if not cluster_exprs:
            return query

        # Transform query for each CLUSTER expression
        for cluster_expr in cluster_exprs:
            query = self._transform_for_cluster(query, cluster_expr)

        return query

    def _transform_subqueries_in_node(self, node: exp.Expression):
        """Recursively transform subqueries within an expression node.

        :param node:
            Expression node to search for subqueries
        """
        # Find and transform any Subquery nodes
        for subquery in node.find_all(exp.Subquery):
            if isinstance(subquery.this, exp.Select):
                transformed = self.transform(subquery.this)
                subquery.set("this", transformed)

    def _find_cluster_expressions(self, query: exp.Select) -> list[GIQLCluster]:
        """Find all CLUSTER expressions in query.

        :param query:
            Query to search
        :return:
            List of CLUSTER expressions
        """
        cluster_exprs = []

        for expression in query.expressions:
            # Check if this is a CLUSTER expression or an alias containing one
            if isinstance(expression, GIQLCluster):
                cluster_exprs.append(expression)
            elif isinstance(expression, exp.Alias):
                if isinstance(expression.this, GIQLCluster):
                    cluster_exprs.append(expression.this)

        return cluster_exprs

    def _transform_for_cluster(
        self, query: exp.Select, cluster_expr: GIQLCluster
    ) -> exp.Select:
        """Transform query to compute CLUSTER using CTEs.

        :param query:
            Original query
        :param cluster_expr:
            CLUSTER expression to transform
        :return:
            Transformed query with CTEs
        """
        # Extract CLUSTER parameters
        distance_expr = cluster_expr.args.get("distance")

        # Handle distance parameter - could be int literal or None
        if distance_expr:
            if isinstance(distance_expr, exp.Literal):
                distance = int(distance_expr.this)
            else:
                # Try to extract as string and convert
                try:
                    distance = int(str(distance_expr.this))
                except (ValueError, AttributeError):
                    distance = 0
        else:
            distance = 0

        stranded_expr = cluster_expr.args.get("stranded")
        if stranded_expr:
            # Handle different types of boolean expressions
            if isinstance(stranded_expr, exp.Boolean):
                stranded = stranded_expr.this
            elif isinstance(stranded_expr, exp.Literal):
                stranded = str(stranded_expr.this).upper() == "TRUE"
            else:
                # Try to extract the value as a string
                stranded = str(stranded_expr).upper() in ("TRUE", "1", "YES")
        else:
            stranded = False

        # Get column names from table config or use defaults
        chrom_col, start_col, end_col, strand_col = self._get_genomic_columns(query)

        # Build partition clause
        partition_cols = [exp.column(chrom_col, quoted=True)]
        if stranded:
            partition_cols.append(exp.column(strand_col, quoted=True))

        # Build ORDER BY for window
        order_by = [exp.Ordered(this=exp.column(start_col, quoted=True))]

        # Create LAG window spec
        lag_window = exp.Window(
            this=exp.Anonymous(
                this="LAG", expressions=[exp.column(end_col, quoted=True)]
            ),
            partition_by=partition_cols,
            order=exp.Order(expressions=order_by),
        )

        # Add distance offset if specified
        if distance > 0:
            lag_with_distance = exp.Add(
                this=lag_window, expression=exp.Literal.number(distance)
            )
        else:
            lag_with_distance = lag_window

        # Create CASE expression for is_new_cluster
        case_expr = exp.Case(
            ifs=[
                exp.If(
                    this=exp.GTE(
                        this=lag_with_distance,
                        expression=exp.column(start_col, quoted=True),
                    ),
                    true=exp.Literal.number(0),
                )
            ],
            default=exp.Literal.number(1),
        )

        # Build CTE SELECT expressions (all original except CLUSTER, plus is_new_cluster)
        cte_expressions = []
        for expression in query.expressions:
            # Skip CLUSTER expressions
            if isinstance(expression, GIQLCluster):
                continue
            elif isinstance(expression, exp.Alias) and isinstance(
                expression.this, GIQLCluster
            ):
                continue
            else:
                cte_expressions.append(expression)

        # Ensure required columns for window functions are included
        required_cols = {chrom_col, start_col, end_col}
        if stranded:
            required_cols.add(strand_col)

        # Check if required columns are already in the select list
        selected_cols = set()
        for expr in cte_expressions:
            if isinstance(expr, exp.Column):
                selected_cols.add(expr.name)
            elif isinstance(expr, exp.Alias):
                # Don't count aliases as the source column
                pass
            elif isinstance(expr, exp.Star):
                # SELECT * includes all columns
                selected_cols = required_cols  # Assume all are covered
                break

        # Add missing required columns
        for col in required_cols - selected_cols:
            cte_expressions.append(exp.column(col, quoted=True))

        # Add is_new_cluster calculation
        cte_expressions.append(exp.alias_(case_expr, "is_new_cluster", quoted=False))

        # Build CTE query
        cte_select = exp.Select()
        cte_select.select(*cte_expressions, copy=False)

        # Copy FROM, WHERE, GROUP BY, HAVING from original (but not ORDER BY)
        # Use copy() to avoid sharing references between queries
        if query.args.get("from_"):
            from_clause = query.args["from_"].copy()
            cte_select.set("from_", from_clause)
        if query.args.get("where"):
            cte_select.set("where", query.args["where"].copy())
        if query.args.get("group"):
            cte_select.set("group", query.args["group"].copy())
        if query.args.get("having"):
            cte_select.set("having", query.args["having"].copy())

        # Create outer query with SUM over is_new_cluster
        sum_window = exp.Window(
            this=exp.Sum(this=exp.column("is_new_cluster")),
            partition_by=partition_cols,
            order=exp.Order(expressions=order_by),
        )

        # Build outer SELECT expressions (replace CLUSTER with SUM)
        new_expressions = []
        for expression in query.expressions:
            if isinstance(expression, GIQLCluster):
                new_expressions.append(sum_window)
            elif isinstance(expression, exp.Alias) and isinstance(
                expression.this, GIQLCluster
            ):
                # Keep the alias but replace the expression
                new_expressions.append(
                    exp.alias_(sum_window, expression.alias, quoted=False)
                )
            else:
                new_expressions.append(expression)

        # Build new query
        new_query = exp.Select()
        new_query.select(*new_expressions, copy=False)

        # Wrap CTE in subquery and set as FROM clause
        subquery = exp.Subquery(
            this=cte_select,
            alias=exp.TableAlias(this=exp.Identifier(this="lag_calc")),
        )
        new_query.from_(subquery, copy=False)

        # Copy ORDER BY from original to outer query
        if query.args.get("order"):
            new_query.order_by(*query.args["order"].expressions, copy=False)

        return new_query


class MergeTransformer:
    """Transforms queries containing MERGE into GROUP BY queries.

    MERGE combines overlapping intervals using CLUSTER + aggregation:

        SELECT MERGE(interval) FROM features

    Into:

        WITH clustered AS (
            SELECT *, CLUSTER(interval) AS __giql_cluster_id FROM features
        )
        SELECT
            chromosome,
            MIN(start_pos) AS start_pos,
            MAX(end_pos) AS end_pos
        FROM clustered
        GROUP BY chromosome, __giql_cluster_id
        ORDER BY chromosome, start_pos
    """

    def __init__(self, tables: Tables):
        """Initialize transformer.

        :param tables:
            Table configurations for column mapping
        """
        self.tables = tables
        self.cluster_transformer = ClusterTransformer(tables)

    def transform(self, query: exp.Expression) -> exp.Expression:
        """Transform query if it contains MERGE expressions.

        :param query:
            Parsed query AST
        :return:
            Transformed query AST
        """
        if not isinstance(query, exp.Select):
            return query

        # First, recursively transform any CTEs that might contain MERGE
        if query.args.get("with_"):
            cte = query.args["with_"]
            for cte_expr in cte.expressions:
                if isinstance(cte_expr, exp.CTE):
                    # Transform the CTE's subquery
                    cte_expr.set("this", self.transform(cte_expr.this))

        # Recursively transform subqueries in FROM clause
        if query.args.get("from_"):
            from_clause = query.args["from_"]
            self._transform_subqueries_in_node(from_clause)

        # Recursively transform subqueries in JOIN clauses
        if query.args.get("joins"):
            for join in query.args["joins"]:
                self._transform_subqueries_in_node(join)

        # Recursively transform subqueries in WHERE clause
        if query.args.get("where"):
            self._transform_subqueries_in_node(query.args["where"])

        # Find all MERGE expressions in the SELECT clause
        merge_exprs = self._find_merge_expressions(query)

        if not merge_exprs:
            return query

        # For now, support only one MERGE expression
        if len(merge_exprs) > 1:
            raise ValueError("Multiple MERGE expressions not yet supported")

        merge_expr = merge_exprs[0]
        return self._transform_for_merge(query, merge_expr)

    def _transform_subqueries_in_node(self, node: exp.Expression):
        """Recursively transform subqueries within an expression node.

        :param node:
            Expression node to search for subqueries
        """
        # Find and transform any Subquery nodes
        for subquery in node.find_all(exp.Subquery):
            if isinstance(subquery.this, exp.Select):
                transformed = self.transform(subquery.this)
                subquery.set("this", transformed)

    def _find_merge_expressions(self, query: exp.Select) -> list[GIQLMerge]:
        """Find all MERGE expressions in query.

        :param query:
            Query to search
        :return:
            List of MERGE expressions
        """
        merge_exprs = []

        for expression in query.expressions:
            if isinstance(expression, GIQLMerge):
                merge_exprs.append(expression)
            elif isinstance(expression, exp.Alias):
                if isinstance(expression.this, GIQLMerge):
                    merge_exprs.append(expression.this)

        return merge_exprs

    def _transform_for_merge(
        self, query: exp.Select, merge_expr: GIQLMerge
    ) -> exp.Select:
        """Transform query to compute MERGE using CLUSTER + GROUP BY.

        :param query:
            Original query
        :param merge_expr:
            MERGE expression to transform
        :return:
            Transformed query with clustering and aggregation
        """
        # Extract MERGE parameters (same as CLUSTER)
        distance_expr = merge_expr.args.get("distance")
        stranded_expr = merge_expr.args.get("stranded")

        # Get column names from table config or use defaults
        (
            chrom_col,
            start_col,
            end_col,
            strand_col,
        ) = self.cluster_transformer._get_genomic_columns(query)

        # Build CLUSTER expression with same parameters
        cluster_kwargs = {"this": merge_expr.this}
        if distance_expr:
            cluster_kwargs["distance"] = distance_expr
        if stranded_expr:
            cluster_kwargs["stranded"] = stranded_expr

        cluster_expr = GIQLCluster(**cluster_kwargs)

        # Create intermediate query with CLUSTER
        # Start with original query's FROM/WHERE/etc
        cluster_query = exp.Select()
        cluster_query.select(exp.Star(), copy=False)
        cluster_query.select(
            exp.alias_(cluster_expr, "__giql_cluster_id", quoted=False),
            append=True,
            copy=False,
        )

        # Copy FROM, WHERE from original
        # Use copy() to avoid sharing references between queries
        if query.args.get("from_"):
            cluster_query.set("from_", query.args["from_"].copy())
        if query.args.get("where"):
            cluster_query.set("where", query.args["where"].copy())

        # Apply CLUSTER transformation to get the CTE-based query
        cluster_query = self.cluster_transformer.transform(cluster_query)

        # Build GROUP BY columns
        group_by_cols = [exp.column(chrom_col)]

        # Handle stranded parameter
        if stranded_expr:
            if isinstance(stranded_expr, exp.Boolean):
                stranded = stranded_expr.this
            elif isinstance(stranded_expr, exp.Literal):
                stranded = str(stranded_expr.this).upper() == "TRUE"
            else:
                stranded = str(stranded_expr).upper() in ("TRUE", "1", "YES")
        else:
            stranded = False

        if stranded:
            group_by_cols.append(exp.column(strand_col, quoted=True))

        group_by_cols.append(exp.column("__giql_cluster_id"))

        # Build SELECT expressions for merged intervals
        select_exprs = []

        # Add group-by columns (non-aggregated)
        select_exprs.append(exp.column(chrom_col, quoted=True))
        if stranded:
            select_exprs.append(exp.column(strand_col, quoted=True))

        # Add merged interval bounds
        select_exprs.append(
            exp.alias_(
                exp.Min(this=exp.column(start_col, quoted=True)), start_col, quoted=False
            )
        )
        select_exprs.append(
            exp.alias_(
                exp.Max(this=exp.column(end_col, quoted=True)), end_col, quoted=False
            )
        )

        # Process other columns from original SELECT
        for expression in query.expressions:
            # Skip the MERGE expression itself
            if isinstance(expression, GIQLMerge):
                continue
            elif isinstance(expression, exp.Alias) and isinstance(
                expression.this, GIQLMerge
            ):
                continue
            # Include other columns (they should be aggregates or in GROUP BY)
            else:
                select_exprs.append(expression)

        # Build final query
        final_query = exp.Select()
        final_query.select(*select_exprs, copy=False)

        # FROM the clustered subquery
        subquery = exp.Subquery(
            this=cluster_query,
            alias=exp.TableAlias(this=exp.Identifier(this="clustered")),
        )
        final_query.from_(subquery, copy=False)

        # Add GROUP BY
        final_query.group_by(*group_by_cols, copy=False)

        # Add ORDER BY (chromosome, start)
        final_query.order_by(
            exp.Ordered(this=exp.column(chrom_col, quoted=True)), copy=False
        )
        final_query.order_by(
            exp.Ordered(this=exp.column(start_col, quoted=True)), append=True, copy=False
        )

        # Preserve any existing CTEs from the original query
        if query.args.get("with_"):
            final_query.set("with_", query.args["with_"].copy())

        return final_query


class CoverageTransformer:
    """Transforms queries containing COVERAGE into binned coverage queries.

    COVERAGE tiles the genome into fixed-width bins and aggregates overlapping
    intervals per bin:

        SELECT COVERAGE(interval, 1000) FROM features

    Into:

        WITH __giql_bins AS (
            SELECT chrom, bin_start AS start, bin_start + 1000 AS "end"
            FROM (
                SELECT DISTINCT chrom, MAX("end") AS __max_end
                FROM features GROUP BY chrom
            ) AS __giql_chroms,
            LATERAL generate_series(0, __max_end, 1000) AS t(bin_start)
        )
        SELECT bins.chrom, bins.start, bins."end", COUNT(source.*)
        FROM __giql_bins AS bins
        LEFT JOIN features AS source
          ON source.start < bins."end"
          AND source."end" > bins.start
          AND source.chrom = bins.chrom
        GROUP BY bins.chrom, bins.start, bins."end"
        ORDER BY bins.chrom, bins.start
    """

    def __init__(self, tables: Tables):
        """Initialize transformer.

        :param tables:
            Table configurations for column mapping
        """
        self.tables = tables

    def _get_table_name(self, query: exp.Select) -> str | None:
        """Extract table name from query's FROM clause.

        :param query:
            Query to extract table name from
        :return:
            Table name if FROM contains a simple table, None otherwise
        """
        from_clause = query.args.get("from_")
        if not from_clause:
            return None
        if isinstance(from_clause.this, exp.Table):
            return from_clause.this.name
        return None

    def _get_table_alias(self, query: exp.Select) -> str | None:
        """Extract table alias from query's FROM clause.

        :param query:
            Query to extract alias from
        :return:
            Table alias if present, None otherwise
        """
        from_clause = query.args.get("from_")
        if not from_clause:
            return None
        if isinstance(from_clause.this, exp.Table):
            return from_clause.this.alias
        return None

    def _get_genomic_columns(self, query: exp.Select) -> tuple[str, str, str]:
        """Get genomic column names from table config or defaults.

        :param query:
            Query to extract table and column info from
        :return:
            Tuple of (chrom_col, start_col, end_col)
        """
        table_name = self._get_table_name(query)

        chrom_col = DEFAULT_CHROM_COL
        start_col = DEFAULT_START_COL
        end_col = DEFAULT_END_COL

        if table_name:
            table = self.tables.get(table_name)
            if table:
                chrom_col = table.chrom_col
                start_col = table.start_col
                end_col = table.end_col

        return chrom_col, start_col, end_col

    def transform(self, query: exp.Expression) -> exp.Expression:
        """Transform query if it contains COVERAGE expressions.

        :param query:
            Parsed query AST
        :return:
            Transformed query AST
        """
        if not isinstance(query, exp.Select):
            return query

        # Recursively transform CTEs
        if query.args.get("with_"):
            cte = query.args["with_"]
            for cte_expr in cte.expressions:
                if isinstance(cte_expr, exp.CTE):
                    cte_expr.set("this", self.transform(cte_expr.this))

        # Recursively transform subqueries in FROM/JOIN/WHERE
        for key in ("from_", "where"):
            if query.args.get(key):
                self._transform_subqueries_in_node(query.args[key])
        if query.args.get("joins"):
            for join in query.args["joins"]:
                self._transform_subqueries_in_node(join)

        # Find COVERAGE expressions in SELECT
        coverage_exprs = self._find_coverage_expressions(query)
        if not coverage_exprs:
            return query

        if len(coverage_exprs) > 1:
            raise ValueError("Multiple COVERAGE expressions not yet supported")

        return self._transform_for_coverage(query, coverage_exprs[0])

    def _transform_subqueries_in_node(self, node: exp.Expression):
        """Recursively transform subqueries within an expression node.

        :param node:
            Expression node to search for subqueries
        """
        for subquery in node.find_all(exp.Subquery):
            if isinstance(subquery.this, exp.Select):
                transformed = self.transform(subquery.this)
                subquery.set("this", transformed)

    def _find_coverage_expressions(self, query: exp.Select) -> list[GIQLCoverage]:
        """Find all COVERAGE expressions in query.

        :param query:
            Query to search
        :return:
            List of COVERAGE expressions
        """
        coverage_exprs = []
        for expression in query.expressions:
            if isinstance(expression, GIQLCoverage):
                coverage_exprs.append(expression)
            elif isinstance(expression, exp.Alias):
                if isinstance(expression.this, GIQLCoverage):
                    coverage_exprs.append(expression.this)
        return coverage_exprs

    def _transform_for_coverage(
        self, query: exp.Select, coverage_expr: GIQLCoverage
    ) -> exp.Select:
        """Transform query to compute COVERAGE using bins CTE + JOIN + GROUP BY.

        :param query:
            Original query
        :param coverage_expr:
            COVERAGE expression to transform
        :return:
            Transformed query
        """
        # Extract parameters
        resolution_expr = coverage_expr.args.get("resolution")
        if isinstance(resolution_expr, exp.Literal):
            resolution = int(resolution_expr.this)
        else:
            try:
                resolution = int(str(resolution_expr.this))
            except (ValueError, AttributeError):
                raise ValueError("COVERAGE resolution must be an integer literal")

        stat_expr = coverage_expr.args.get("stat")
        if stat_expr:
            if isinstance(stat_expr, exp.Literal):
                stat = stat_expr.this.strip("'\"").lower()
            else:
                stat = str(stat_expr).strip("'\"").lower()
        else:
            stat = "count"

        if stat not in COVERAGE_STAT_MAP:
            raise ValueError(
                f"Unknown COVERAGE stat '{stat}'. "
                f"Must be one of: {', '.join(COVERAGE_STAT_MAP)}"
            )

        sql_agg = COVERAGE_STAT_MAP[stat]

        # Extract target parameter
        target_expr = coverage_expr.args.get("target")
        if target_expr:
            if isinstance(target_expr, exp.Literal):
                target_col = target_expr.this.strip("'\"")
            else:
                target_col = str(target_expr).strip("'\"")
        else:
            target_col = None

        # Get column names and table info
        chrom_col, start_col, end_col = self._get_genomic_columns(query)
        table_name = self._get_table_name(query)
        table_alias = self._get_table_alias(query)
        source_ref = table_alias or table_name or "source"

        # Build __giql_chroms subquery:
        #   SELECT DISTINCT chrom, MAX("end") AS __max_end FROM <table> GROUP BY chrom
        chroms_select = exp.Select()
        chroms_select.select(
            exp.column(chrom_col, quoted=True),
            copy=False,
        )
        chroms_select.select(
            exp.alias_(
                exp.Max(this=exp.column(end_col, quoted=True)),
                "__max_end",
                quoted=False,
            ),
            append=True,
            copy=False,
        )

        if table_name:
            chroms_select.from_(exp.to_table(table_name), copy=False)

        # Apply WHERE from original query to the chroms subquery too,
        # qualifying unqualified column references with the table name
        if query.args.get("where"):
            chroms_where = query.args["where"].copy()
            if table_name:
                for col in chroms_where.find_all(exp.Column):
                    if not col.table:
                        col.set("table", exp.Identifier(this=table_name))
            chroms_select.set("where", chroms_where)

        chroms_select.group_by(exp.column(chrom_col, quoted=True), copy=False)

        chroms_subquery = exp.Subquery(
            this=chroms_select,
            alias=exp.TableAlias(this=exp.Identifier(this="__giql_chroms")),
        )

        # Build bins CTE using raw SQL for generate_series + LATERAL
        # since SQLGlot doesn't natively support generate_series
        bins_select = exp.Select()
        bins_select.select(
            exp.column(chrom_col, table="__giql_chroms", quoted=True),
            copy=False,
        )
        bins_select.select(
            exp.alias_(
                exp.column("bin_start"),
                start_col,
                quoted=True,
            ),
            append=True,
            copy=False,
        )
        bins_select.select(
            exp.alias_(
                exp.Add(
                    this=exp.column("bin_start"),
                    expression=exp.Literal.number(resolution),
                ),
                end_col,
                quoted=True,
            ),
            append=True,
            copy=False,
        )

        # FROM __giql_chroms subquery
        bins_select.from_(chroms_subquery, copy=False)

        # CROSS JOIN LATERAL generate_series(0, __max_end, resolution) AS t(bin_start)
        lateral_join = exp.Join(
            this=exp.Lateral(
                this=exp.Anonymous(
                    this="generate_series",
                    expressions=[
                        exp.Literal.number(0),
                        exp.column("__max_end"),
                        exp.Literal.number(resolution),
                    ],
                ),
                alias=exp.TableAlias(
                    this=exp.Identifier(this="t"),
                    columns=[exp.Identifier(this="bin_start")],
                ),
            ),
            kind="CROSS",
        )
        bins_select.append("joins", lateral_join)

        # Wrap bins_select as a CTE named __giql_bins
        bins_cte = exp.CTE(
            this=bins_select,
            alias=exp.TableAlias(this=exp.Identifier(this="__giql_bins")),
        )
        with_clause = exp.With(expressions=[bins_cte])

        # Build the aggregate expression
        if stat == "count":
            if target_col:
                agg_expr = exp.Anonymous(
                    this="COUNT",
                    expressions=[
                        exp.column(target_col, table=source_ref, quoted=True),
                    ],
                )
            else:
                agg_expr = exp.Anonymous(
                    this="COUNT",
                    expressions=[
                        exp.Column(
                            this=exp.Star(),
                            table=exp.Identifier(this=source_ref),
                        )
                    ],
                )
        else:
            if target_col:
                agg_expr = exp.Anonymous(
                    this=sql_agg,
                    expressions=[
                        exp.column(target_col, table=source_ref, quoted=True),
                    ],
                )
            else:
                agg_expr = exp.Anonymous(
                    this=sql_agg,
                    expressions=[
                        exp.Sub(
                            this=exp.column(end_col, table=source_ref, quoted=True),
                            expression=exp.column(
                                start_col, table=source_ref, quoted=True
                            ),
                        )
                    ],
                )

        # Build main SELECT
        final_query = exp.Select()

        # Add bin coordinate columns
        final_query.select(
            exp.column(chrom_col, table="bins", quoted=True),
            copy=False,
        )
        final_query.select(
            exp.column(start_col, table="bins", quoted=True),
            append=True,
            copy=False,
        )
        final_query.select(
            exp.column(end_col, table="bins", quoted=True),
            append=True,
            copy=False,
        )

        # Replace COVERAGE(...) in select list with aggregate, and add other columns
        for expression in query.expressions:
            if isinstance(expression, GIQLCoverage):
                final_query.select(
                    exp.alias_(agg_expr, "value", quoted=False),
                    append=True,
                    copy=False,
                )
            elif isinstance(expression, exp.Alias) and isinstance(
                expression.this, GIQLCoverage
            ):
                final_query.select(
                    exp.alias_(agg_expr, expression.alias, quoted=False),
                    append=True,
                    copy=False,
                )
            else:
                final_query.select(expression, append=True, copy=False)

        # FROM __giql_bins AS bins
        final_query.from_(
            exp.Table(
                this=exp.Identifier(this="__giql_bins"),
                alias=exp.TableAlias(this=exp.Identifier(this="bins")),
            ),
            copy=False,
        )

        # LEFT JOIN source ON overlap conditions
        source_table = exp.to_table(table_name) if table_name else exp.to_table("source")
        source_table.set("alias", exp.TableAlias(this=exp.Identifier(this=source_ref)))

        join_condition = exp.And(
            this=exp.And(
                this=exp.LT(
                    this=exp.column(start_col, table=source_ref, quoted=True),
                    expression=exp.column(end_col, table="bins", quoted=True),
                ),
                expression=exp.GT(
                    this=exp.column(end_col, table=source_ref, quoted=True),
                    expression=exp.column(start_col, table="bins", quoted=True),
                ),
            ),
            expression=exp.EQ(
                this=exp.column(chrom_col, table=source_ref, quoted=True),
                expression=exp.column(chrom_col, table="bins", quoted=True),
            ),
        )

        # Merge original WHERE into the JOIN ON condition so that
        # LEFT JOIN still produces zero-coverage bins (WHERE would filter
        # them out because source columns are NULL for non-matching bins)
        if query.args.get("where"):
            where_condition = query.args["where"].this.copy()
            # Qualify unqualified column references with source_ref
            for col in where_condition.find_all(exp.Column):
                if not col.table:
                    col.set("table", exp.Identifier(this=source_ref))
            join_condition = exp.And(
                this=join_condition,
                expression=where_condition,
            )

        left_join = exp.Join(
            this=source_table,
            on=join_condition,
            kind="LEFT",
        )
        final_query.append("joins", left_join)

        # GROUP BY bins.chrom, bins.start, bins.end
        final_query.group_by(
            exp.column(chrom_col, table="bins", quoted=True),
            copy=False,
        )
        final_query.group_by(
            exp.column(start_col, table="bins", quoted=True),
            append=True,
            copy=False,
        )
        final_query.group_by(
            exp.column(end_col, table="bins", quoted=True),
            append=True,
            copy=False,
        )

        # ORDER BY bins.chrom, bins.start
        final_query.order_by(
            exp.Ordered(this=exp.column(chrom_col, table="bins", quoted=True)),
            copy=False,
        )
        final_query.order_by(
            exp.Ordered(this=exp.column(start_col, table="bins", quoted=True)),
            append=True,
            copy=False,
        )

        # Attach the WITH clause
        final_query.set("with_", with_clause)

        return final_query
