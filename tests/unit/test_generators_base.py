"""Tests for BaseGIQLGenerator.

Test specification: specs/test_generators_base.md
Test IDs: BG-001 through BG-020
"""

import pytest
from sqlglot import parse_one

from giql.dialect import GIQLDialect
from giql.generators import BaseGIQLGenerator
from giql.table import Table
from giql.table import Tables


@pytest.fixture
def tables_two():
    """Tables with two tables for column-to-column tests."""
    tables = Tables()
    tables.register("features_a", Table("features_a"))
    tables.register("features_b", Table("features_b"))
    return tables


@pytest.fixture
def tables_peaks_and_genes():
    """Tables with peaks and genes for NEAREST/DISTANCE tests."""
    tables = Tables()
    tables.register("peaks", Table("peaks"))
    tables.register("genes", Table("genes"))
    return tables


def _normalize(sql: str) -> str:
    """Collapse whitespace for easier assertion."""
    return " ".join(sql.split())


class TestBaseGIQLGenerator:
    """Tests for BaseGIQLGenerator class (BG-001 to BG-020)."""

    # ------------------------------------------------------------------
    # Instantiation
    # ------------------------------------------------------------------

    def test_bg_001_no_args_defaults(self):
        """
        GIVEN no arguments
        WHEN BaseGIQLGenerator is instantiated
        THEN instance has empty Tables and SUPPORTS_LATERAL is True.
        """
        generator = BaseGIQLGenerator()

        assert generator.tables is not None
        assert generator.SUPPORTS_LATERAL is True
        # Empty tables: looking up any name returns None
        assert generator.tables.get("anything") is None

    def test_bg_002_with_tables(self):
        """
        GIVEN a Tables instance with a registered table
        WHEN BaseGIQLGenerator is instantiated with tables=
        THEN the instance uses the provided tables for column resolution.
        """
        tables = Tables()
        tables.register("peaks", Table("peaks"))
        generator = BaseGIQLGenerator(tables=tables)

        assert generator.tables is tables
        assert "peaks" in generator.tables

    # ------------------------------------------------------------------
    # Spatial predicates
    # ------------------------------------------------------------------

    def test_bg_003_intersects_literal(self):
        """
        GIVEN an Intersects AST node with a literal range 'chr1:1000-2000'
        WHEN generate is called
        THEN output contains chrom = 'chr1' AND start < 2000 AND end > 1000.
        """
        tables = Tables()
        tables.register("peaks", Table("peaks"))
        generator = BaseGIQLGenerator(tables=tables)

        ast = parse_one(
            "SELECT * FROM peaks WHERE interval INTERSECTS 'chr1:1000-2000'",
            dialect=GIQLDialect,
        )
        sql = generator.generate(ast)

        assert "\"chrom\" = 'chr1'" in sql
        assert '"start" < 2000' in sql
        assert '"end" > 1000' in sql

    def test_bg_004_intersects_column_to_column(self, tables_two):
        """
        GIVEN an Intersects AST node with column-to-column (a.interval INTERSECTS b.interval)
        WHEN generate is called
        THEN output contains chrom equality and overlap conditions using both table prefixes.
        """
        generator = BaseGIQLGenerator(tables=tables_two)

        ast = parse_one(
            "SELECT * FROM features_a AS a CROSS JOIN features_b AS b "
            "WHERE a.interval INTERSECTS b.interval",
            dialect=GIQLDialect,
        )
        sql = generator.generate(ast)

        assert 'a."chrom" = b."chrom"' in sql
        assert 'a."start" < b."end"' in sql
        assert 'a."end" > b."start"' in sql

    def test_bg_005_contains_point(self):
        """
        GIVEN a Contains AST node with a point range 'chr1:1500'
        WHEN generate is called
        THEN output contains point containment predicate.
        """
        generator = BaseGIQLGenerator()

        ast = parse_one(
            "SELECT * FROM peaks WHERE interval CONTAINS 'chr1:1500'",
            dialect=GIQLDialect,
        )
        sql = generator.generate(ast)

        assert "\"chrom\" = 'chr1'" in sql
        assert '"start" <= 1500' in sql
        assert '"end" > 1500' in sql

    def test_bg_006_contains_range(self):
        """
        GIVEN a Contains AST node with a range 'chr1:1000-2000'
        WHEN generate is called
        THEN output contains range containment predicate.
        """
        generator = BaseGIQLGenerator()

        ast = parse_one(
            "SELECT * FROM peaks WHERE interval CONTAINS 'chr1:1000-2000'",
            dialect=GIQLDialect,
        )
        sql = generator.generate(ast)

        assert "\"chrom\" = 'chr1'" in sql
        assert '"start" <= 1000' in sql
        assert '"end" >= 2000' in sql

    def test_bg_007_within_range(self):
        """
        GIVEN a Within AST node with a range 'chr1:1000-5000'
        WHEN generate is called
        THEN output contains within predicate.
        """
        generator = BaseGIQLGenerator()

        ast = parse_one(
            "SELECT * FROM peaks WHERE interval WITHIN 'chr1:1000-5000'",
            dialect=GIQLDialect,
        )
        sql = generator.generate(ast)

        assert "\"chrom\" = 'chr1'" in sql
        assert '"start" >= 1000' in sql
        assert '"end" <= 5000' in sql

    # ------------------------------------------------------------------
    # Spatial set predicates
    # ------------------------------------------------------------------

    def test_bg_008_intersects_any(self):
        """
        GIVEN a SpatialSetPredicate with INTERSECTS ANY and two ranges
        WHEN generate is called
        THEN output contains two conditions joined by OR.
        """
        generator = BaseGIQLGenerator()

        ast = parse_one(
            "SELECT * FROM peaks "
            "WHERE interval INTERSECTS ANY('chr1:1000-2000', 'chr1:5000-6000')",
            dialect=GIQLDialect,
        )
        sql = generator.generate(ast)

        assert " OR " in sql
        assert '"end" > 1000' in sql
        assert '"end" > 5000' in sql

    def test_bg_009_intersects_all(self):
        """
        GIVEN a SpatialSetPredicate with INTERSECTS ALL and two ranges
        WHEN generate is called
        THEN output contains two conditions joined by AND.
        """
        generator = BaseGIQLGenerator()

        ast = parse_one(
            "SELECT * FROM peaks "
            "WHERE interval INTERSECTS ALL('chr1:1000-2000', 'chr1:1500-1800')",
            dialect=GIQLDialect,
        )
        sql = generator.generate(ast)

        # The outer WHERE already has AND, but the set predicate wraps
        # its conditions in parens joined by AND.
        norm = _normalize(sql)
        # Both range predicates should appear
        assert '"start" < 2000' in sql
        assert '"start" < 1800' in sql
        # They are joined by AND (inside the set predicate parentheses)
        # Check the pattern: one condition AND another condition
        idx_first = norm.index('"start" < 2000')
        idx_second = norm.index('"start" < 1800')
        between = norm[idx_first:idx_second]
        assert "AND" in between

    # ------------------------------------------------------------------
    # DISTANCE
    # ------------------------------------------------------------------

    def test_bg_010_distance_basic(self, tables_two):
        """
        GIVEN a GIQLDistance node with two column references
        WHEN generate is called
        THEN output contains CASE WHEN with chromosome check, overlap check, and distance calculations.
        """
        generator = BaseGIQLGenerator(tables=tables_two)

        ast = parse_one(
            "SELECT DISTANCE(a.interval, b.interval) AS dist "
            "FROM features_a a CROSS JOIN features_b b",
            dialect=GIQLDialect,
        )
        sql = generator.generate(ast)

        assert 'a."chrom" != b."chrom" THEN NULL' in sql
        assert "THEN 0" in sql
        assert 'b."start" - a."end"' in sql
        assert 'a."start" - b."end"' in sql
        assert sql.startswith("SELECT CASE WHEN")

    def test_bg_011_distance_stranded(self, tables_two):
        """
        GIVEN a GIQLDistance node with stranded := true
        WHEN generate is called
        THEN output contains strand NULL checks and strand flip logic.
        """
        generator = BaseGIQLGenerator(tables=tables_two)

        ast = parse_one(
            "SELECT DISTANCE(a.interval, b.interval, stranded := true) AS dist "
            "FROM features_a a CROSS JOIN features_b b",
            dialect=GIQLDialect,
        )
        sql = generator.generate(ast)

        assert 'a."strand" IS NULL' in sql
        assert 'b."strand" IS NULL' in sql
        assert "a.\"strand\" = '.'" in sql
        assert "a.\"strand\" = '?'" in sql
        assert "a.\"strand\" = '-'" in sql

    def test_bg_012_distance_signed(self, tables_two):
        """
        GIVEN a GIQLDistance node with signed := true
        WHEN generate is called
        THEN output contains signed distance (negative for upstream).
        """
        generator = BaseGIQLGenerator(tables=tables_two)

        ast = parse_one(
            "SELECT DISTANCE(a.interval, b.interval, signed := true) AS dist "
            "FROM features_a a CROSS JOIN features_b b",
            dialect=GIQLDialect,
        )
        sql = generator.generate(ast)

        # Signed: ELSE branch has negative sign
        assert "-(" in sql
        # Unsigned ELSE would be (a."start" - b."end") without negation
        # Signed ELSE is -(a."start" - b."end")
        assert '-(a."start" - b."end")' in sql

    def test_bg_013_distance_stranded_and_signed(self, tables_two):
        """
        GIVEN a GIQLDistance node with stranded := true and signed := true
        WHEN generate is called
        THEN output contains both strand flip and signed distance.
        """
        generator = BaseGIQLGenerator(tables=tables_two)

        ast = parse_one(
            "SELECT DISTANCE(a.interval, b.interval, stranded := true, signed := true) AS dist "
            "FROM features_a a CROSS JOIN features_b b",
            dialect=GIQLDialect,
        )
        sql = generator.generate(ast)

        # Should have strand NULL checks
        assert 'a."strand" IS NULL' in sql
        # Should have strand flip
        assert "a.\"strand\" = '-'" in sql
        # Stranded+signed: the ELSE for '-' strand flips sign differently
        # from stranded-only
        # In stranded+signed: ELSE WHEN strand='-' THEN (a.start - b.end)
        # In stranded-only:   ELSE WHEN strand='-' THEN -(a.start - b.end)
        assert '(a."start" - b."end")' in sql
        assert '-(a."start" - b."end")' in sql

    def test_bg_014_distance_closed_intervals(self):
        """
        GIVEN tables with interval_type="closed" for one table
        WHEN generate is called for a DISTANCE expression
        THEN output contains '+ 1' gap adjustment.
        """
        tables = Tables()
        tables.register("bed_a", Table("bed_a", interval_type="closed"))
        tables.register("bed_b", Table("bed_b", interval_type="closed"))
        generator = BaseGIQLGenerator(tables=tables)

        ast = parse_one(
            "SELECT DISTANCE(a.interval, b.interval) AS dist "
            "FROM bed_a a CROSS JOIN bed_b b",
            dialect=GIQLDialect,
        )
        sql = generator.generate(ast)

        assert "+ 1)" in sql

    # ------------------------------------------------------------------
    # NEAREST
    # ------------------------------------------------------------------

    def test_bg_015_nearest_standalone(self, tables_peaks_and_genes):
        """
        GIVEN a GIQLNearest node with explicit reference (standalone mode)
        WHEN generate is called
        THEN output is a subquery with WHERE, ORDER BY ABS(distance), LIMIT.
        """
        generator = BaseGIQLGenerator(tables=tables_peaks_and_genes)

        ast = parse_one(
            "SELECT * FROM NEAREST(genes, reference := 'chr1:1000-2000')",
            dialect=GIQLDialect,
        )
        sql = generator.generate(ast)
        norm = _normalize(sql)

        assert "WHERE" in norm
        assert "ORDER BY ABS(" in norm
        assert "LIMIT 1" in norm
        assert "'chr1' = genes.\"chrom\"" in sql
        assert "AS distance" in sql

    def test_bg_016_nearest_k5(self, tables_peaks_and_genes):
        """
        GIVEN a GIQLNearest node with k := 5
        WHEN generate is called
        THEN output has LIMIT 5.
        """
        generator = BaseGIQLGenerator(tables=tables_peaks_and_genes)

        ast = parse_one(
            "SELECT * FROM NEAREST(genes, reference := 'chr1:1000-2000', k := 5)",
            dialect=GIQLDialect,
        )
        sql = generator.generate(ast)

        assert "LIMIT 5" in sql

    def test_bg_017_nearest_max_distance(self, tables_peaks_and_genes):
        """
        GIVEN a GIQLNearest node with max_distance := 100000
        WHEN generate is called
        THEN the distance threshold appears in the WHERE clause.
        """
        generator = BaseGIQLGenerator(tables=tables_peaks_and_genes)

        ast = parse_one(
            "SELECT * FROM NEAREST(genes, reference := 'chr1:1000-2000', max_distance := 100000)",
            dialect=GIQLDialect,
        )
        sql = generator.generate(ast)
        norm = _normalize(sql)

        assert "100000" in norm
        assert "<= 100000" in norm

    def test_bg_018_nearest_correlated_lateral(self, tables_peaks_and_genes):
        """
        GIVEN a GIQLNearest node in correlated mode (no standalone reference, in LATERAL context)
        WHEN generate is called
        THEN output is a LATERAL-compatible subquery referencing the outer table columns.
        """
        generator = BaseGIQLGenerator(tables=tables_peaks_and_genes)

        ast = parse_one(
            "SELECT * FROM peaks "
            "CROSS JOIN LATERAL NEAREST(genes, reference := peaks.interval, k := 3)",
            dialect=GIQLDialect,
        )
        sql = generator.generate(ast)
        norm = _normalize(sql)

        assert "LATERAL" in norm
        assert 'peaks."chrom"' in sql
        assert 'genes."chrom"' in sql
        assert "LIMIT 3" in sql

    def test_bg_019_nearest_stranded(self, tables_peaks_and_genes):
        """
        GIVEN a GIQLNearest node with stranded := true
        WHEN generate is called
        THEN output includes strand matching in WHERE clause.
        """
        generator = BaseGIQLGenerator(tables=tables_peaks_and_genes)

        ast = parse_one(
            "SELECT * FROM peaks "
            "CROSS JOIN LATERAL NEAREST(genes, reference := peaks.interval, k := 3, stranded := true)",
            dialect=GIQLDialect,
        )
        sql = generator.generate(ast)

        assert 'peaks."strand"' in sql
        assert 'genes."strand"' in sql
        # Strand matching in WHERE
        assert 'peaks."strand" = genes."strand"' in sql

    # ------------------------------------------------------------------
    # SELECT override
    # ------------------------------------------------------------------

    def test_bg_020_select_alias_mapping(self):
        """
        GIVEN a SELECT with aliased FROM and JOIN tables
        WHEN generate is called
        THEN alias-to-table mapping is built correctly, verified through correct column resolution in a spatial op.
        """
        tables = Tables()
        tables.register("features_a", Table("features_a"))
        tables.register("features_b", Table("features_b"))
        generator = BaseGIQLGenerator(tables=tables)

        ast = parse_one(
            "SELECT * FROM features_a AS a "
            "JOIN features_b AS b ON a.id = b.id "
            "WHERE a.interval INTERSECTS b.interval",
            dialect=GIQLDialect,
        )
        sql = generator.generate(ast)

        # The aliases 'a' and 'b' should resolve to the registered tables
        # and produce correctly qualified column references
        assert 'a."chrom" = b."chrom"' in sql
        assert 'a."start" < b."end"' in sql
        assert 'a."end" > b."start"' in sql
