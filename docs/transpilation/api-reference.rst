API Reference
=============

Transpilation
-------------

.. autofunction:: giql.transpile.transpile

.. autoclass:: giql.table.Table

Extension hook
--------------

The registry-based extension point for custom targets and operator expanders
(see :doc:`extending`). Every symbol below is re-exported from the top-level
``giql`` package.

.. autofunction:: giql.expander.register

.. autodata:: giql.expander.REGISTRY

.. autoclass:: giql.expander.ExpanderRegistry
   :members: register, register_target, target, resolve, has_override,
             unregister, clear, snapshot, restore

.. autoclass:: giql.expander.RegistrySnapshot
   :members:

.. autoclass:: giql.expander.ExpansionContext
   :members:

.. autoclass:: giql.expander.OperatorExpander
   :members:

Targets and capabilities
------------------------

.. autoclass:: giql.targets.Target
   :members:

.. autoclass:: giql.targets.Capabilities
   :members:

.. autoclass:: giql.targets.GenericTarget

.. autoclass:: giql.targets.DuckDBTarget

.. autoclass:: giql.targets.DataFusionTarget
