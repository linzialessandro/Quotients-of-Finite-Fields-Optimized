"""
Bitmask encodings of finite hyperoperation table cells.

Inspired by Talotti's binary-string representation of finite magmata
(subsets of an n-element set as integers in ``[0, 2^n - 1]``), specialised
here to the additive tables of Krasner quotient hyperfields.

Labels live in ``{-1} ∪ {0, ..., r-1}`` (``ZERO = -1`` for the additive
identity).  We map them to bit positions:

  - bit 0  ↔  ZERO (additive 0)
  - bit i+1 ↔ label i  (i = 0..r-1)

For ``r + 1 ≤ 63`` the mask fits in a machine word; larger ``r`` still works
via Python's arbitrary-precision integers.

Reference
---------
E. Talotti, *Hyperstructures in computer science*, 2025
(local PDF: ``papers/hyperstructure_comp_science.pdf``).
"""

from __future__ import annotations

from typing import Iterable

from quotient_hyperfields.hyperfield import ZERO


def label_bit(label: int) -> int:
    """Bit index for a hyperfield label (ZERO → 0, label i → i+1)."""
    if label == ZERO:
        return 0
    if label < 0:
        raise ValueError(f"invalid label {label}")
    return label + 1


def labels_to_mask(labels: Iterable[int]) -> int:
    """Encode a set of labels as a bitmask integer."""
    mask = 0
    for x in labels:
        mask |= 1 << label_bit(x)
    return mask


def mask_to_labels(mask: int) -> frozenset[int]:
    """Decode a bitmask into a frozenset of labels."""
    if mask < 0:
        raise ValueError("mask must be non-negative")
    out: set[int] = set()
    bit = 0
    while mask:
        if mask & 1:
            out.add(ZERO if bit == 0 else bit - 1)
        mask >>= 1
        bit += 1
    return frozenset(out)


def mask_contains(mask: int, label: int) -> bool:
    """True iff ``label`` is set in ``mask``."""
    return bool(mask & (1 << label_bit(label)))
