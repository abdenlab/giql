"""Built-in operator expanders for epic #137.

Importing this package registers every built-in expander as a side effect:
each submodule decorates its expander(s) with ``@register(...)`` at import
time, and this package imports all of them. The import is wired once (in
:mod:`giql.transpile`) so the process-wide ``REGISTRY`` is populated before the
first transpile.

New operator modules are picked up automatically: drop a ``<operator>.py`` into
this package and it is imported here without editing this file.

Modules whose name starts with ``_`` are skipped (private helpers, not
expanders). Submodules import in :func:`pkgutil.iter_modules` order, which sets
last-write-wins resolution-order precedence for overlapping registrations; an
import error here aborts the whole package import by design (a broken built-in
expander must not be silently skipped).
"""

from __future__ import annotations

import importlib
import pkgutil

for _module_info in pkgutil.iter_modules(__path__):
    if _module_info.name.startswith("_"):
        continue
    importlib.import_module(f"{__name__}.{_module_info.name}")
