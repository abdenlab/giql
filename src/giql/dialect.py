"""Custom SQL dialect with genomic extensions.

This module defines the GIQL dialect, which extends standard SQL with
spatial operators for genomic interval queries.
"""

from typing import Final

from sqlglot.dialects import Dialect
from sqlglot.parser import Parser
from sqlglot.tokens import Tokenizer
from sqlglot.tokens import TokenType

from giql.expressions import Contains
from giql.expressions import GIQLCluster
from giql.expressions import GIQLMerge
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
        """Tokenizer with genomic keywords.

        Extends the base tokenizer to recognize GIQL spatial operators
        (INTERSECTS, CONTAINS, WITHIN).
        """

        KEYWORDS = {
            **Tokenizer.KEYWORDS,
            INTERSECTS: getattr(TokenType, INTERSECTS),
            CONTAINS: getattr(TokenType, CONTAINS),
            WITHIN: getattr(TokenType, WITHIN),
        }

    class Parser(Parser):
        """Parser with genomic predicate support."""

        FUNCTIONS = {
            **Parser.FUNCTIONS,
            "CLUSTER": GIQLCluster.from_arg_list,
            "MERGE": GIQLMerge.from_arg_list,
        }

        def _parse_comparison(self):
            """Override to handle spatial operators.

            :return:
                Parsed spatial expression or falls back to parent's comparison parsing
            """
            return self._parse_spatial() or super()._parse_comparison()

        def _parse_spatial(self):
            """Parse spatial predicates.

            Handles:
                - column INTERSECTS 'chr1:1000-2000'
                - column INTERSECTS ANY('chr1:1000-2000', 'chr1:5000-6000')
                - column CONTAINS 'chr1:1500'
                - column WITHIN 'chr1:1000-5000'

            :return:
                Parsed spatial expression or None if no spatial operator found
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
            """Parse right side of spatial predicate.

            :param left:
                Left side expression (column reference)
            :param operator:
                Spatial operator token (INTERSECTS, CONTAINS, WITHIN)
            :param expr_class:
                Expression class to instantiate (Intersects, Contains, Within)
            :return:
                Parsed spatial predicate expression
            """
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
