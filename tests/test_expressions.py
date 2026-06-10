"""Unit tests for the GIQL expression-node slot descriptors.

Tests pin the declarative ``SlotSpec`` contract on the operator expression
classes: which slots each operator declares and which of them are table-shaped
reference slots that the ``ResolveOperatorRefs`` pass resolves.
"""

from giql.expressions import Contains
from giql.expressions import GIQLDisjoin
from giql.expressions import GIQLDistance
from giql.expressions import GIQLNearest
from giql.expressions import Intersects
from giql.expressions import SlotSpec
from giql.expressions import SpatialSetPredicate
from giql.expressions import Within

OPERATOR_CLASSES = (
    GIQLDisjoin,
    GIQLNearest,
    GIQLDistance,
    Intersects,
    Contains,
    Within,
    SpatialSetPredicate,
)


class TestSlotSpec:
    """Tests for the SlotSpec slot descriptor."""

    def test_is_ref_slot_with_table_trichotomy_shapes(self):
        """Test is_ref_slot for a slot accepting only table shapes.

        Given:
            A SlotSpec whose accepted shapes are the registered-table / CTE /
            subquery trichotomy.
        When:
            Reading its is_ref_slot property.
        Then:
            It should be True.
        """
        # Arrange
        spec = SlotSpec("reference", frozenset({"registered_table", "cte", "subquery"}))

        # Act & assert
        assert spec.is_ref_slot is True

    def test_is_ref_slot_with_column_shape(self):
        """Test is_ref_slot for a slot accepting a non-table shape.

        Given:
            A SlotSpec whose accepted shapes include a column reference.
        When:
            Reading its is_ref_slot property.
        Then:
            It should be False.
        """
        # Arrange
        spec = SlotSpec("expression", frozenset({"literal_range", "column"}))

        # Act & assert
        assert spec.is_ref_slot is False

    def test_is_ref_slot_with_empty_accepts(self):
        """Test is_ref_slot for a slot accepting no shapes.

        Given:
            A SlotSpec with an empty accepted-shapes set.
        When:
            Reading its is_ref_slot property.
        Then:
            It should be False.
        """
        # Arrange
        spec = SlotSpec("this", frozenset())

        # Act & assert
        assert spec.is_ref_slot is False

    def test_giql_slots_declared_on_every_operator_class(self):
        """Test the declarative slot contract across the operator classes.

        Given:
            The seven GIQL operator expression classes.
        When:
            Reading each class's GIQL_SLOTS declaration.
        Then:
            It should declare a tuple of SlotSpecs on every operator, with
            exactly DISJOIN's this/reference and NEAREST's this as the
            table-shaped reference slots.
        """
        # Arrange
        expected_ref_slots = {
            (GIQLDisjoin, "this"),
            (GIQLDisjoin, "reference"),
            (GIQLNearest, "this"),
        }

        # Act
        declared_ref_slots = {
            (cls, spec.arg)
            for cls in OPERATOR_CLASSES
            for spec in cls.GIQL_SLOTS
            if spec.is_ref_slot
        }

        # Assert
        for cls in OPERATOR_CLASSES:
            assert isinstance(cls.GIQL_SLOTS, tuple)
            assert all(isinstance(spec, SlotSpec) for spec in cls.GIQL_SLOTS)
        assert declared_ref_slots == expected_ref_slots
