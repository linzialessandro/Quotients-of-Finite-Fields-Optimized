"""
Research library for Krasner quotient hyperfields of finite fields.

Objects are hyperfields of the form F_q / G, where G is the unique multiplicative
subgroup of F_q^x of a given order d (equivalently, of a given index r = (q-1)/d).
The hyperfield has order r + 1.

Primary references
------------------
- Baker–Jin, arXiv:1912.11496 — isomorphism of large finite-field quotients
- Kędzierski–Linzi–Stojałowska — characteristic and C-characteristic
"""

from quotient_hyperfields.hyperfield import QuotientHyperfield
from quotient_hyperfields.isomorphism import (
    IsoResult,
    are_isomorphic,
    are_isomorphic_baker_jin,
    are_isomorphic_general,
    baker_jin_bound,
    baker_jin_bound_safe,
    baker_jin_bound_sharp,
    is_baker_jin_large,
    weil_excess,
)
from quotient_hyperfields.census import (
    classify_hyperfields,
    generate_triples,
)

# experiments is intentionally not imported here so that
# `python -m quotient_hyperfields.experiments` / the qh-open-questions
# entry point can run cleanly without a double-import warning.

__all__ = [
    "QuotientHyperfield",
    "IsoResult",
    "are_isomorphic",
    "are_isomorphic_baker_jin",
    "are_isomorphic_general",
    "baker_jin_bound",
    "baker_jin_bound_safe",
    "baker_jin_bound_sharp",
    "is_baker_jin_large",
    "weil_excess",
    "generate_triples",
    "classify_hyperfields",
]

__version__ = "0.1.0"
