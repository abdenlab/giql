# Test Interface Contracts

This directory contains interface contracts (protocols and abstract base classes) for the bedtools integration test suite.

## Purpose

These contracts define the expected behavior of test infrastructure components without specifying implementation details. They serve as:

1. **Design documentation**: Clear specification of component responsibilities
2. **Type checking**: Enable static type checking with mypy/pyright
3. **Testability**: Allow mocking and testing of test infrastructure itself
4. **Flexibility**: Multiple implementations can satisfy same contract

## Files

### `test_interfaces.py`

Defines all protocols and contracts for the test suite:

- **`IntervalGeneratorProtocol`**: Interface for generating simulated genomic datasets
- **`BedtoolsWrapperProtocol`**: Interface for executing bedtools commands
- **`ResultComparatorProtocol`**: Interface for comparing GIQL and bedtools results
- **`DatabaseFixtureProtocol`**: Interface for database test fixtures
- **`TestScenarioProtocol`**: Interface for test scenario definitions
- **Error classes**: Custom exceptions for test execution failures

## Usage

### Implementing a Protocol

```python
from contracts.test_interfaces import IntervalGeneratorProtocol

class IntervalGenerator:
    """Concrete implementation of IntervalGeneratorProtocol."""

    def __init__(self, seed: int = 42):
        self.rng = random.Random(seed)

    def generate_overlapping(
        self,
        chromosome: str,
        count: int
    ) -> List[Tuple[str, int, int, str, int, str]]:
        # Implementation details...
        pass
```

### Type Hints

```python
def run_test(
    generator: IntervalGeneratorProtocol,
    bedtools: BedtoolsWrapperProtocol,
    comparator: ResultComparatorProtocol
) -> ComparisonResult:
    # Function accepts any object implementing these protocols
    pass
```

### Using ComparisonResult

```python
def test_intersect_basic():
    """Test basic intersection operation."""
    # Arrange
    # ...

    # Act
    result = compare_results(giql_output, bedtools_output)

    # Assert
    assert result, result.failure_message()
```

## Design Principles

1. **Protocol over inheritance**: Use `Protocol` for structural subtyping, not ABC
2. **Clear contracts**: Each protocol has single, well-defined responsibility
3. **Error handling**: Custom exceptions for different failure modes
4. **Testability**: Contracts enable easy mocking for unit tests

## Contract Guarantees

Each protocol implementation MUST:

- Validate inputs and raise appropriate exceptions for invalid data
- Return consistent types as specified in protocol
- Be deterministic when given same inputs (especially with seeded generation)
- Clean up resources (temp files, connections) in error cases

## Evolution Guidelines

When modifying contracts:

1. **Backward compatibility**: Don't break existing implementations
2. **Optional parameters**: Add new parameters as optional with defaults
3. **New methods**: Add as optional (not required by protocol)
4. **Breaking changes**: Require major version bump and migration guide
