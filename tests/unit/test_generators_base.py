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

    def test___init___should_use_defaults_when_no_args(self):
        """Test __init__ uses default state when no arguments are supplied.

        Given:
            No arguments
        When:
            BaseGIQLGenerator is instantiated
        Then:
            It should have empty Tables and SUPPORTS_LATERAL set to True
        """
        # Arrange / Act
        generator = BaseGIQLGenerator()

        # Assert
        assert generator.tables is not None
        assert generator.SUPPORTS_LATERAL is True
        # Empty tables: looking up any name returns None
        assert generator.tables.get("anything") is None

    def test___init___should_use_provided_tables_when_given(self):
        """Test __init__ adopts a caller-supplied Tables instance.

        Given:
            A Tables instance with a registered table
        When:
            BaseGIQLGenerator is instantiated with tables=
        Then:
            It should use the provided tables for column resolution
        """
        # Arrange
        tables = Tables()
        tables.register("peaks", Table("peaks"))

        # Act
        generator = BaseGIQLGenerator(tables=tables)

        # Assert
        assert generator.tables is tables
        assert "peaks" in generator.tables

    # ------------------------------------------------------------------
    # Spatial predicates
    # ------------------------------------------------------------------

    def test_generate_should_emit_overlap_conditions_when_intersects_literal(self):
        """Test generate emits overlap SQL for an INTERSECTS literal range.

        Given:
            An Intersects AST node with a literal range 'chr1:1000-2000'
        When:
            generate is called
        Then:
            It should contain chrom = 'chr1' AND start < 2000 AND end > 1000
        """
        # Arrange
        tables = Tables()
        tables.register("peaks", Table("peaks"))
        generator = BaseGIQLGenerator(tables=tables)
        ast = parse_one(
            "SELECT * FROM peaks WHERE interval INTERSECTS 'chr1:1000-2000'",
            dialect=GIQLDialect,
        )

        # Act
        sql = generator.generate(ast)

        # Assert
        assert "\"chrom\" = 'chr1'" in sql
        assert '"start" < 2000' in sql
        assert '"end" > 1000' in sql

    def test_generate_should_emit_qualified_overlap_when_intersects_column_to_column(self, tables_two):
        """Test generate emits table-qualified overlap for column-to-column INTERSECTS.

        Given:
            An Intersects AST node with column-to-column (a.interval INTERSECTS b.interval)
        When:
            generate is called
        Then:
            It should contain chrom equality and overlap conditions using both table prefixes
        """
        # Arrange
        generator = BaseGIQLGenerator(tables=tables_two)
        ast = parse_one(
            "SELECT * FROM features_a AS a CROSS JOIN features_b AS b "
            "WHERE a.interval INTERSECTS b.interval",
            dialect=GIQLDialect,
        )

        # Act
        sql = generator.generate(ast)

        # Assert
        assert 'a."chrom" = b."chrom"' in sql
        assert 'a."start" < b."end"' in sql
        assert 'a."end" > b."start"' in sql

    def test_generate_should_emit_point_containment_when_contains_point(self):
        """Test generate emits point containment SQL when CONTAINS targets a point.

        Given:
            A Contains AST node with a point range 'chr1:1500'
        When:
            generate is called
        Then:
            It should contain point containment predicate
        """
        # Arrange
        generator = BaseGIQLGenerator()
        ast = parse_one(
            "SELECT * FROM peaks WHERE interval CONTAINS 'chr1:1500'",
            dialect=GIQLDialect,
        )

        # Act
        sql = generator.generate(ast)

        # Assert
        assert "\"chrom\" = 'chr1'" in sql
        assert '"start" <= 1500' in sql
        assert '"end" > 1500' in sql

    def test_generate_should_emit_range_containment_when_contains_range(self):
        """Test generate emits range containment SQL when CONTAINS targets a range.

        Given:
            A Contains AST node with a range 'chr1:1000-2000'
        When:
            generate is called
        Then:
            It should contain range containment predicate
        """
        # Arrange
        generator = BaseGIQLGenerator()
        ast = parse_one(
            "SELECT * FROM peaks WHERE interval CONTAINS 'chr1:1000-2000'",
            dialect=GIQLDialect,
        )

        # Act
        sql = generator.generate(ast)

        # Assert
        assert "\"chrom\" = 'chr1'" in sql
        assert '"start" <= 1000' in sql
        assert '"end" >= 2000' in sql

    def test_generate_should_emit_within_predicate_when_within_range(self):
        """Test generate emits within-range SQL when the predicate is WITHIN.

        Given:
            A Within AST node with a range 'chr1:1000-5000'
        When:
            generate is called
        Then:
            It should contain within predicate
        """
        # Arrange
        generator = BaseGIQLGenerator()
        ast = parse_one(
            "SELECT * FROM peaks WHERE interval WITHIN 'chr1:1000-5000'",
            dialect=GIQLDialect,
        )

        # Act
        sql = generator.generate(ast)

        # Assert
        assert "\"chrom\" = 'chr1'" in sql
        assert '"start" >= 1000' in sql
        assert '"end" <= 5000' in sql

    # ------------------------------------------------------------------
    # Spatial set predicates
    # ------------------------------------------------------------------

    def test_generate_should_join_with_or_when_intersects_any(self):
        """Test generate joins predicates with OR for INTERSECTS ANY.

        Given:
            A SpatialSetPredicate with INTERSECTS ANY and two ranges
        When:
            generate is called
        Then:
            It should contain two conditions joined by OR
        """
        # Arrange
        generator = BaseGIQLGenerator()
        ast = parse_one(
            "SELECT * FROM peaks "
            "WHERE interval INTERSECTS ANY('chr1:1000-2000', 'chr1:5000-6000')",
            dialect=GIQLDialect,
        )

        # Act
        sql = generator.generate(ast)

        # Assert
        assert " OR " in sql
        assert '"end" > 1000' in sql
        assert '"end" > 5000' in sql

    def test_generate_should_join_with_and_when_intersects_all(self):
        """Test generate joins predicates with AND for INTERSECTS ALL.

        Given:
            A SpatialSetPredicate with INTERSECTS ALL and two ranges
        When:
            generate is called
        Then:
            It should contain two conditions joined by AND
        """
        # Arrange
        generator = BaseGIQLGenerator()
        ast = parse_one(
            "SELECT * FROM peaks "
            "WHERE interval INTERSECTS ALL('chr1:1000-2000', 'chr1:1500-1800')",
            dialect=GIQLDialect,
        )

        # Act
        sql = generator.generate(ast)

        # Assert
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

    def test_generate_should_emit_case_when_distance_basic(self, tables_two):
        """Test generate emits a CASE WHEN expression for basic DISTANCE.

        Given:
            A GIQLDistance node with two column references
        When:
            generate is called
        Then:
            It should contain CASE WHEN with chromosome check, overlap check, and distance calculations
        """
        # Arrange
        generator = BaseGIQLGenerator(tables=tables_two)
        ast = parse_one(
            "SELECT DISTANCE(a.interval, b.interval) AS dist "
            "FROM features_a a CROSS JOIN features_b b",
            dialect=GIQLDialect,
        )

        # Act
        sql = generator.generate(ast)

        # Assert
        assert 'a."chrom" != b."chrom" THEN NULL' in sql
        assert "THEN 0" in sql
        assert 'b."start" - a."end"' in sql
        assert 'a."start" - b."end"' in sql
        assert sql.startswith("SELECT CASE WHEN")

    def test_generate_should_emit_strand_logic_when_distance_stranded(self, tables_two):
        """Test generate emits strand NULL checks and flip logic when DISTANCE is stranded.

        Given:
            A GIQLDistance node with stranded := true
        When:
            generate is called
        Then:
            It should contain strand NULL checks and strand flip logic
        """
        # Arrange
        generator = BaseGIQLGenerator(tables=tables_two)
        ast = parse_one(
            "SELECT DISTANCE(a.interval, b.interval, stranded := true) AS dist "
            "FROM features_a a CROSS JOIN features_b b",
            dialect=GIQLDialect,
        )

        # Act
        sql = generator.generate(ast)

        # Assert
        assert 'a."strand" IS NULL' in sql
        assert 'b."strand" IS NULL' in sql
        assert "a.\"strand\" = '.'" in sql
        assert "a.\"strand\" = '?'" in sql
        assert "a.\"strand\" = '-'" in sql

    def test_generate_should_emit_signed_distance_when_distance_signed(self, tables_two):
        """Test generate emits a negated upstream branch when DISTANCE is signed.

        Given:
            A GIQLDistance node with signed := true
        When:
            generate is called
        Then:
            It should contain signed distance (negative for upstream)
        """
        # Arrange
        generator = BaseGIQLGenerator(tables=tables_two)
        ast = parse_one(
            "SELECT DISTANCE(a.interval, b.interval, signed := true) AS dist "
            "FROM features_a a CROSS JOIN features_b b",
            dialect=GIQLDialect,
        )

        # Act
        sql = generator.generate(ast)

        # Assert
        # Signed: ELSE branch has negative sign
        assert "-(" in sql
        # Unsigned ELSE would be (a."start" - b."end") without negation
        # Signed ELSE is -(a."start" - b."end")
        assert '-(a."start" - b."end")' in sql

    def test_generate_should_combine_strand_and_sign_when_distance_stranded_and_signed(self, tables_two):
        """Test generate combines strand flipping and signed output when both flags are set.

        Given:
            A GIQLDistance node with stranded := true and signed := true
        When:
            generate is called
        Then:
            It should contain both strand flip and signed distance
        """
        # Arrange
        generator = BaseGIQLGenerator(tables=tables_two)
        ast = parse_one(
            "SELECT DISTANCE(a.interval, b.interval, stranded := true, signed := true) AS dist "
            "FROM features_a a CROSS JOIN features_b b",
            dialect=GIQLDialect,
        )

        # Act
        sql = generator.generate(ast)

        # Assert
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

    def test_generate_should_add_gap_adjustment_when_distance_uses_closed_intervals(self):
        """Test generate adds a +1 gap adjustment for closed-interval DISTANCE.

        Given:
            Tables with interval_type="closed" for one table
        When:
            generate is called for a DISTANCE expression
        Then:
            It should contain '+ 1' gap adjustment
        """
        # Arrange
        tables = Tables()
        tables.register("bed_a", Table("bed_a", interval_type="closed"))
        tables.register("bed_b", Table("bed_b", interval_type="closed"))
        generator = BaseGIQLGenerator(tables=tables)
        ast = parse_one(
            "SELECT DISTANCE(a.interval, b.interval) AS dist "
            "FROM bed_a a CROSS JOIN bed_b b",
            dialect=GIQLDialect,
        )

        # Act
        sql = generator.generate(ast)

        # Assert
        assert "+ 1)" in sql

    # ------------------------------------------------------------------
    # NEAREST
    # ------------------------------------------------------------------

    def test_generate_should_emit_order_by_and_limit_when_nearest_standalone(self, tables_peaks_and_genes):
        """Test generate emits an ORDER BY / LIMIT subquery for standalone NEAREST.

        Given:
            A GIQLNearest node with explicit reference (standalone mode)
        When:
            generate is called
        Then:
            It should produce a subquery with WHERE, ORDER BY ABS(distance), and LIMIT
        """
        # Arrange
        generator = BaseGIQLGenerator(tables=tables_peaks_and_genes)
        ast = parse_one(
            "SELECT * FROM NEAREST(genes, reference := 'chr1:1000-2000')",
            dialect=GIQLDialect,
        )

        # Act
        sql = generator.generate(ast)
        norm = _normalize(sql)

        # Assert
        assert "WHERE" in norm
        assert "ORDER BY ABS(" in norm
        assert "LIMIT 1" in norm
        assert "'chr1' = genes.\"chrom\"" in sql
        assert "AS distance" in sql

    def test_generate_should_limit_five_when_nearest_k_is_five(self, tables_peaks_and_genes):
        """Test generate applies LIMIT 5 when NEAREST is given k := 5.

        Given:
            A GIQLNearest node with k := 5
        When:
            generate is called
        Then:
            It should produce LIMIT 5
        """
        # Arrange
        generator = BaseGIQLGenerator(tables=tables_peaks_and_genes)
        ast = parse_one(
            "SELECT * FROM NEAREST(genes, reference := 'chr1:1000-2000', k := 5)",
            dialect=GIQLDialect,
        )

        # Act
        sql = generator.generate(ast)

        # Assert
        assert "LIMIT 5" in sql

    def test_generate_should_embed_threshold_when_nearest_max_distance(self, tables_peaks_and_genes):
        """Test generate embeds the max_distance threshold in the WHERE clause.

        Given:
            A GIQLNearest node with max_distance := 100000
        When:
            generate is called
        Then:
            It should place the distance threshold in the WHERE clause
        """
        # Arrange
        generator = BaseGIQLGenerator(tables=tables_peaks_and_genes)
        ast = parse_one(
            "SELECT * FROM NEAREST(genes, reference := 'chr1:1000-2000', max_distance := 100000)",
            dialect=GIQLDialect,
        )

        # Act
        sql = generator.generate(ast)
        norm = _normalize(sql)

        # Assert
        assert "100000" in norm
        assert "<= 100000" in norm

    def test_generate_should_reference_outer_columns_when_nearest_correlated_lateral(self, tables_peaks_and_genes):
        """Test generate emits a LATERAL-compatible subquery referencing outer columns.

        Given:
            A GIQLNearest node in correlated mode (no standalone reference, in LATERAL context)
        When:
            generate is called
        Then:
            It should produce a LATERAL-compatible subquery referencing the outer table columns
        """
        # Arrange
        generator = BaseGIQLGenerator(tables=tables_peaks_and_genes)
        ast = parse_one(
            "SELECT * FROM peaks "
            "CROSS JOIN LATERAL NEAREST(genes, reference := peaks.interval, k := 3)",
            dialect=GIQLDialect,
        )

        # Act
        sql = generator.generate(ast)
        norm = _normalize(sql)

        # Assert
        assert "LATERAL" in norm
        assert 'peaks."chrom"' in sql
        assert 'genes."chrom"' in sql
        assert "LIMIT 3" in sql

    def test_generate_should_match_strand_when_nearest_stranded(self, tables_peaks_and_genes):
        """Test generate includes strand matching in WHERE when NEAREST is stranded.

        Given:
            A GIQLNearest node with stranded := true
        When:
            generate is called
        Then:
            It should include strand matching in the WHERE clause
        """
        # Arrange
        generator = BaseGIQLGenerator(tables=tables_peaks_and_genes)
        ast = parse_one(
            "SELECT * FROM peaks "
            "CROSS JOIN LATERAL NEAREST(genes, reference := peaks.interval, k := 3, stranded := true)",
            dialect=GIQLDialect,
        )

        # Act
        sql = generator.generate(ast)

        # Assert
        assert 'peaks."strand"' in sql
        assert 'genes."strand"' in sql
        # Strand matching in WHERE
        assert 'peaks."strand" = genes."strand"' in sql

    # ------------------------------------------------------------------
    # SELECT override
    # ------------------------------------------------------------------

    def test_generate_should_resolve_aliases_when_select_has_alias_mapping(self):
        """Test generate resolves FROM/JOIN aliases to registered tables.

        Given:
            A SELECT with aliased FROM and JOIN tables
        When:
            generate is called
        Then:
            It should build alias-to-table mapping correctly, verified through correct column resolution in a spatial op
        """
        # Arrange
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

        # Act
        sql = generator.generate(ast)

        # Assert
        # The aliases 'a' and 'b' should resolve to the registered tables
        # and produce correctly qualified column references
        assert 'a."chrom" = b."chrom"' in sql
        assert 'a."start" < b."end"' in sql
        assert 'a."end" > b."start"' in sql
