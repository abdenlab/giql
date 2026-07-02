"""Shared operator-parameter coercions for the expanders.

The private ``_``-prefix keeps this module out of the ``expanders`` package
auto-discovery (see :mod:`giql.expanders`), so it is a plain helper module rather
than a registered expander.
"""

from __future__ import annotations

from sqlglot import exp


def coerce_bool_param(param: exp.Expression | None) -> bool:
    """Coerce an optional operator boolean argument to a Python ``bool``.

    The ``stranded`` / ``signed`` keyword arguments of DISTANCE and NEAREST are
    read the same way: an absent argument is ``False``, an
    :class:`sqlglot.exp.Boolean` node yields ``bool(node.this)``, and anything
    else is truthy when its rendered text is ``TRUE`` / ``1`` / ``YES``. The
    :func:`bool` coercion guarantees a real ``bool`` regardless of what the
    parser stored on the node.

    Parameters
    ----------
    param : exp.Expression | None
        The ``stranded`` or ``signed`` argument node, or ``None`` if absent.

    Returns
    -------
    bool
        The coerced Python boolean (``False`` when the argument is absent).
    """
    if not param:
        return False
    if isinstance(param, exp.Boolean):
        return bool(param.this)
    return str(param).upper() in ("TRUE", "1", "YES")
