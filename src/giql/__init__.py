"""GIQL - Genomic Interval Query Language.

A SQL dialect for genomic range queries.
"""

from giql.expander import REGISTRY
from giql.expander import ExpanderRegistry
from giql.expander import ExpansionContext
from giql.expander import OperatorExpander
from giql.expander import StatementFinalizer
from giql.expander import register
from giql.table import Table
from giql.targets import Capabilities
from giql.targets import DataFusionTarget
from giql.targets import DuckDBTarget
from giql.targets import GenericTarget
from giql.targets import Target
from giql.transpile import transpile

__version__ = "0.1.0"


__all__ = [
    "Table",
    "transpile",
    # Extension hook: register custom targets and override operator expanders.
    "register",
    "REGISTRY",
    "ExpanderRegistry",
    "ExpansionContext",
    "OperatorExpander",
    "StatementFinalizer",
    "Target",
    "Capabilities",
    "GenericTarget",
    "DuckDBTarget",
    "DataFusionTarget",
]
