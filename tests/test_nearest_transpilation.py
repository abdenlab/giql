"""Transpilation tests for NEAREST operator SQL generation.

Tests verify that NEAREST() is correctly transpiled to SQL
(LATERAL joins for correlated queries, ORDER BY + LIMIT for standalone).
"""

import pytest
from sqlglot import parse_one

from giql import Table
from giql.dialect import GIQLDialect
from giql.generators import BaseGIQLGenerator
from giql.table import Tables


@pytest.fixture
def tables_with_peaks_and_genes():
    """Tables container with peaks and genes tables."""
    tables = Tables()
    tables.register("peaks", Table("peaks"))
    tables.register("genes", Table("genes"))
    return tables


class TestNearestTranspilation:
    """Tests for NEAREST transpilation to SQL."""

    def test_nearest_basic_k3(self, tables_with_peaks_and_genes):
        """
        GIVEN a GIQL query with NEAREST(genes, k := 3)
        WHEN transpiling to SQL
        THEN should generate LATERAL join with DISTANCE and LIMIT 3
        """
        sql = """
        SELECT *
        FROM peaks
        CROSS JOIN LATERAL NEAREST(genes, reference := peaks.interval, k := 3)
        """

        ast = parse_one(sql, dialect=GIQLDialect)
        generator = BaseGIQLGenerator(tables=tables_with_peaks_and_genes)
        output = generator.generate(ast)

        # Expectations:
        # - LATERAL subquery
        # - DISTANCE(...) AS distance in SELECT
        # - WHERE peaks.chromosome = genes.chromosome (pre-filter)
        # - ORDER BY distance
        # - LIMIT 3
        assert "LATERAL" in output.upper()
        assert "CASE" in output or "DISTANCE" in output  # Distance calculation
        assert " AS distance" in output or " as distance" in output.lower()
        assert "LIMIT 3" in output
        assert "ORDER BY" in output

    def test_nearest_with_max_distance(self, tables_with_peaks_and_genes):
        """
        GIVEN a GIQL query with NEAREST(genes, k := 5, max_distance := 100000)
        WHEN transpiling to SQL
        THEN should generate LATERAL join with distance filter
        """
        sql = """
        SELECT *
        FROM peaks
        CROSS JOIN LATERAL NEAREST(genes, reference := peaks.interval, k := 5, max_distance := 100000)
        """

        ast = parse_one(sql, dialect=GIQLDialect)
        generator = BaseGIQLGenerator(tables=tables_with_peaks_and_genes)
        output = generator.generate(ast)

        # Expectations:
        # - LATERAL subquery
        # - Distance filter: <= 100000
        # - LIMIT 5
        assert "LATERAL" in output.upper()
        assert "100000" in output
        assert "LIMIT 5" in output

    def test_nearest_standalone_literal(self, tables_with_peaks_and_genes):
        """
        GIVEN a GIQL query with literal reference NEAREST(genes, reference := 'chr1:1000-2000', k := 3)
        WHEN transpiling to SQL
        THEN should generate standalone query without LATERAL
        """
        sql = """
        SELECT *
        FROM NEAREST(genes, reference := 'chr1:1000-2000', k := 3)
        """

        ast = parse_one(sql, dialect=GIQLDialect)
        generator = BaseGIQLGenerator(tables=tables_with_peaks_and_genes)
        output = generator.generate(ast)

        # Expectations:
        # - No LATERAL (standalone mode)
        # - Distance calculation with literal 'chr1', 1000, 2000
        # - ORDER BY distance
        # - LIMIT 3
        assert "LATERAL" not in output.upper()
        assert "chr1" in output.lower()
        assert "LIMIT 3" in output

    def test_nearest_with_stranded(self, tables_with_peaks_and_genes):
        """
        GIVEN a GIQL query with NEAREST(genes, k := 3, stranded := true)
        WHEN transpiling to SQL
        THEN should generate SQL with strand filtering
        """
        sql = """
        SELECT *
        FROM peaks
        CROSS JOIN LATERAL NEAREST(genes, reference := peaks.interval, k := 3, stranded := true)
        """

        ast = parse_one(sql, dialect=GIQLDialect)
        generator = BaseGIQLGenerator(tables=tables_with_peaks_and_genes)
        output = generator.generate(ast)

        # Expectations:
        # - LATERAL subquery
        # - Strand filtering in WHERE clause
        # - LIMIT 3
        assert "LATERAL" in output.upper()
        assert "strand" in output.lower()
        assert "LIMIT 3" in output

    def test_nearest_with_signed(self, tables_with_peaks_and_genes):
        """
        GIVEN a GIQL query with NEAREST(genes, k := 3, signed := true)
        WHEN transpiling to SQL
        THEN should generate SQL with signed distance column
            (negative for upstream, positive for downstream)
        """
        sql = """
        SELECT *
        FROM peaks
        CROSS JOIN LATERAL NEAREST(genes, reference := peaks.interval, k := 3, signed := true)
        """

        ast = parse_one(sql, dialect=GIQLDialect)
        generator = BaseGIQLGenerator(tables=tables_with_peaks_and_genes)
        output = generator.generate(ast)

        # Expectations:
        # - LATERAL subquery
        # - Signed distance calculation (includes negation for upstream)
        # - LIMIT 3
        assert "LATERAL" in output.upper()
        assert "LIMIT 3" in output
        # Check for signed distance: the ELSE branch should have a negation
        # for upstream features (B before A)
        assert "ELSE -(" in output, (
            f"Expected signed distance with negation for upstream, got:\n{output}"
        )

    def test_giqlnearest_sql_should_raise_legacy_error_when_target_is_unknown(
        self, tables_with_peaks_and_genes
    ):
        """Test that NEAREST keeps the legacy error wording for an unknown target.

        Given:
            A NEAREST call whose target is neither a registered table nor a CTE.
        When:
            Generating SQL via the base generator.
        Then:
            It should raise ``ValueError`` matching the legacy ``Target table
            'missing' not found`` wording and not borrow DISJOIN-specific
            phrasing — a regression guard for the DISJOIN/NEAREST branching
            introduced for issue #105.
        """
        # Arrange
        sql = (
            "SELECT * FROM peaks "
            "CROSS JOIN LATERAL NEAREST(missing, reference := peaks.interval, k := 1)"
        )
        ast = parse_one(sql, dialect=GIQLDialect)
        generator = BaseGIQLGenerator(tables=tables_with_peaks_and_genes)

        # Act & assert
        with pytest.raises(ValueError) as excinfo:
            generator.generate(ast)
        message = str(excinfo.value)
        assert "Target table 'missing' not found in tables" in message
        assert "DISJOIN target" not in message
        assert "CTE defined in this query" not in message

    def test_giqlnearest_sql_should_reject_cte_target_with_cte_aware_message(
        self, tables_with_peaks_and_genes
    ):
        """Test that NEAREST rejects a CTE target with a CTE-aware message.

        Given:
            A WITH clause defines a CTE named ``__cte`` and a NEAREST call
            names it as its target; only ``peaks`` is a registered table.
        When:
            Generating SQL via the base generator.
        Then:
            It should raise ``ValueError`` whose message names the CTE
            match and directs the user to register the relation as a
            table — the legacy ``not found in tables`` wording masked the
            in-scope CTE case.
        """
        # Arrange
        sql = (
            "WITH __cte AS (SELECT * FROM peaks) "
            "SELECT * FROM peaks CROSS JOIN LATERAL "
            "NEAREST(__cte, reference := peaks.interval, k := 1)"
        )
        ast = parse_one(sql, dialect=GIQLDialect)
        generator = BaseGIQLGenerator(tables=tables_with_peaks_and_genes)

        # Act & assert
        with pytest.raises(ValueError) as excinfo:
            generator.generate(ast)
        message = str(excinfo.value)
        assert "NEAREST target '__cte'" in message
        assert "matches an enclosing CTE" in message
        assert "register the relation as a table" in message

    def test_giqlnearest_sql_should_raise_ambiguity_when_target_is_both_registered_and_cte(
        self, tables_with_peaks_and_genes
    ):
        """Test that NEAREST rejects a target name that is both registered and a CTE.

        Given:
            A query where the NEAREST target name ``peaks`` matches both a
            registered table and an enclosing CTE.
        When:
            Generating SQL via the base generator.
        Then:
            It should raise ``ValueError`` whose message names the
            ambiguity and asks the user to rename one of the two
            bindings, preventing the silent SQL-scoping shadow where the
            CTE rows win at execution time while column expressions are
            built against the registered table's schema.
        """
        # Arrange
        sql = (
            'WITH peaks AS (SELECT \'chr1\' AS chrom, 500 AS "start", 600 AS "end") '
            "SELECT * FROM NEAREST(peaks, reference := 'chr1:50-60', k := 1)"
        )
        ast = parse_one(sql, dialect=GIQLDialect)
        generator = BaseGIQLGenerator(tables=tables_with_peaks_and_genes)

        # Act & assert
        with pytest.raises(ValueError) as excinfo:
            generator.generate(ast)
        message = str(excinfo.value)
        assert "NEAREST target 'peaks'" in message
        assert "resolves to both" in message
        assert "rename one of them" in message

    def test_giqlnearest_sql_should_raise_without_DISJOIN_prefix_when_target_is_qualified(
        self, tables_with_peaks_and_genes
    ):
        """Test that a qualified NEAREST target raises a NEAREST-prefixed error.

        Given:
            A NEAREST query whose target is db-qualified (``db.genes``).
        When:
            Generating SQL via the base generator.
        Then:
            It should raise ``ValueError`` whose message begins with
            ``Target`` (capitalised role label) and does NOT contain
            ``DISJOIN`` — a regression guard against the shared
            ``_extract_bare_name`` helper hardcoding the wrong operator
            name on the NEAREST branch.
        """
        # Arrange
        sql = "SELECT * FROM NEAREST(db.genes, reference := 'chr1:1000-2000', k := 1)"
        ast = parse_one(sql, dialect=GIQLDialect)
        generator = BaseGIQLGenerator(tables=tables_with_peaks_and_genes)

        # Act & assert
        with pytest.raises(ValueError) as excinfo:
            generator.generate(ast)
        message = str(excinfo.value)
        assert "is qualified" in message
        assert message.startswith("Target")
        assert "DISJOIN" not in message
