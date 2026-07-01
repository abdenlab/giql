"""Tests for the top-level ``giql`` public API surface.

#146 promoted the extension hook to a supported public API and removed the
custom generator. These tests pin the package facade: the extension-hook
symbols are importable from the ``giql`` root, are the same objects as their
submodule definitions, and the deleted ``giql.generators`` package is gone.
"""

import importlib

import pytest

import giql

#: The extension-hook symbols #146 added to ``giql.__all__`` (beyond the
#: pre-existing ``DEFAULT_BIN_SIZE`` / ``Table`` / ``transpile``).
_EXTENSION_HOOK_EXPORTS = [
    "register",
    "REGISTRY",
    "ExpanderRegistry",
    "ExpansionContext",
    "OperatorExpander",
    "Target",
    "Capabilities",
    "GenericTarget",
    "DuckDBTarget",
    "DataFusionTarget",
]

#: Each re-export paired with the submodule it is defined in, so the facade can
#: be checked for object identity (not a shadow copy).
_REEXPORT_ORIGINS = [
    ("register", "giql.expander"),
    ("REGISTRY", "giql.expander"),
    ("ExpanderRegistry", "giql.expander"),
    ("ExpansionContext", "giql.expander"),
    ("OperatorExpander", "giql.expander"),
    ("Target", "giql.targets"),
    ("Capabilities", "giql.targets"),
    ("GenericTarget", "giql.targets"),
    ("DuckDBTarget", "giql.targets"),
    ("DataFusionTarget", "giql.targets"),
]


class TestPublicApiSurface:
    """The extension hook is a first-class part of the ``giql`` package API."""

    @pytest.mark.parametrize("name", _EXTENSION_HOOK_EXPORTS)
    def test_extension_hook_symbol_is_importable_from_package_root(self, name):
        """Test that each extension-hook symbol is exported from ``giql``.

        Given:
            An extension-hook name #146 added to ``giql.__all__``.
        When:
            Accessing it as an attribute of the top-level ``giql`` package and
            checking ``giql.__all__``.
        Then:
            It should resolve to a real attribute and be listed in ``__all__``.
        """
        # Act & assert
        assert hasattr(giql, name)
        assert name in giql.__all__

    @pytest.mark.parametrize("name, origin_module", _REEXPORT_ORIGINS)
    def test_reexport_is_identical_to_submodule_origin(self, name, origin_module):
        """Test that a package-root export is the submodule's own object.

        Given:
            A name re-exported from ``giql`` and the submodule that defines it.
        When:
            Comparing ``giql.<name>`` to ``<submodule>.<name>``.
        Then:
            They should be the same object (a re-export, not a shadow copy).
        """
        # Arrange
        submodule = importlib.import_module(origin_module)

        # Act & assert
        assert getattr(giql, name) is getattr(submodule, name)

    def test_generators_package_is_removed(self):
        """Test that the deleted ``giql.generators`` package is gone.

        Given:
            The ``giql.generators`` package that #146 removed.
        When:
            Attempting to import the package.
        Then:
            It should raise ModuleNotFoundError, pinning the deletion so the
            generator cannot silently reappear.
        """
        # Act & assert
        with pytest.raises(ModuleNotFoundError):
            importlib.import_module("giql.generators")
