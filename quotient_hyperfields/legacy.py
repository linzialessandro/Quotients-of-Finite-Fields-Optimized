"""
Quarantined legacy helpers for prime-field quotients GF(p)/G only.

Prefer :class:`quotient_hyperfields.hyperfield.QuotientHyperfield` for all new code.
These remain only for historical notebooks that may still import them.
"""

from __future__ import annotations

# Re-export nothing critical; kept as a marker module.
LEGACY_NOTE = (
    "Legacy GF(p)-only hyperaddition helpers lived in utils.py. "
    "Use QuotientHyperfield instead."
)
