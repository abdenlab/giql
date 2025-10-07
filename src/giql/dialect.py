"""
Custom SQL dialect with genomic extensions.
"""

from typing import Final

from sqlglot.dialects import Dialect
from sqlglot.parser import Parser
from sqlglot.tokens import Tokenizer
from sqlglot.tokens import TokenType

from giql.expressions import Contains
from giql.expressions import Intersects
from giql.expressions import SpatialSetPredicate
from giql.expressions import Within

# Token type constants
INTERSECTS: Final = "INTERSECTS"
CONTAINS: Final = "CONTAINS"
WITHIN: Final = "WITHIN"

# Register custom token types
setattr(TokenType, INTERSECTS, INTERSECTS)
setattr(TokenType, CONTAINS, CONTAINS)
setattr(TokenType, WITHIN, WITHIN)


class GIQLDialect(Dialect):
    """Generic SQL dialect with genomic spatial operators."""

    class Tokenizer(Tokenizer):
        """Tokenizer with genomic keywords."""

        KEYWORDS = {
            **Tokenizer.KEYWORDS,
            INTERSECTS: getattr(TokenType, INTERSECTS),
            CONTAINS: getattr(TokenType, CONTAINS),
            WITHIN: getattr(TokenType, WITHIN),
        }

    class Parser(Parser):
        """Parser with genomic predicate support."""

        def _parse_comparison(self):
            """Override to handle spatial operators."""
            return self._parse_spatial() or super()._parse_comparison()

        def _parse_spatial(self):
            """
            Parse spatial predicates.

            Handles:
                - column INTERSECTS 'chr1:1000-2000'
                - column INTERSECTS ANY('chr1:1000-2000', 'chr1:5000-6000')
            """
            start_index = self._index
            this = self._parse_term()

            if self._match(getattr(TokenType, INTERSECTS)):
                return self._parse_spatial_predicate(this, INTERSECTS, Intersects)
            elif self._match(getattr(TokenType, CONTAINS)):
                return self._parse_spatial_predicate(this, CONTAINS, Contains)
            elif self._match(getattr(TokenType, WITHIN)):
                return self._parse_spatial_predicate(this, WITHIN, Within)

            # No spatial operator found - retreat and return None to allow fallback
            self._retreat(start_index)
            return None

        def _parse_spatial_predicate(self, left, operator, expr_class):
            """Parse right side of spatial predicate."""
            # Check for ANY/ALL quantifier
            if self._match_set((TokenType.ANY, TokenType.ALL, TokenType.SOME)):
                assert self._prev is not None, "Expected token after successful match"
                quantifier = self._prev.text.upper()
                if quantifier == "SOME":
                    quantifier = "ANY"

                # Parse range list
                self._match_l_paren()
                ranges = self._parse_csv(self._parse_expression)
                self._match_r_paren()

                return self.expression(
                    SpatialSetPredicate,
                    this=left,
                    operator=operator,
                    quantifier=quantifier,
                    ranges=ranges,
                )
            else:
                # Simple spatial predicate
                right = self._parse_term()
                return self.expression(expr_class, this=left, expression=right)
