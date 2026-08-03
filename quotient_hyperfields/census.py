"""Enumeration of (p, k, d) triples and classification into iso classes."""

from __future__ import annotations

from quotient_hyperfields.hyperfield import QuotientHyperfield
from quotient_hyperfields.isomorphism import are_isomorphic
from quotient_hyperfields.primes import is_prime


def generate_triples(n: int, max_p: int, max_k: int) -> list[tuple[int, int, int]]:
    """
    Valid triples (p, k, d) with hyperfield order n = (p^k - 1)/d + 1,
    p prime <= max_p, 1 <= k <= max_k.
    """
    if n <= 1:
        return []

    r = n - 1  # index
    triples: list[tuple[int, int, int]] = []
    for p in range(2, max_p + 1):
        if not is_prime(p):
            continue
        for k in range(1, max_k + 1):
            q = p**k
            if (q - 1) % r == 0:
                d = (q - 1) // r
                triples.append((p, k, d))
    return triples


def classify_hyperfields(
    n: int,
    max_p: int,
    max_k: int,
    *,
    method: str = "auto",
) -> list[list[tuple[int, int, int]]]:
    """
    Partition triples of hyperfield order n into isomorphism classes.
    """
    triples = generate_triples(n, max_p, max_k)
    if not triples:
        return []

    # Build once
    objects = [QuotientHyperfield.from_params(p, k, d) for p, k, d in triples]
    assigned = [False] * len(triples)
    classes: list[list[tuple[int, int, int]]] = []

    for i, hi in enumerate(objects):
        if assigned[i]:
            continue
        current = [triples[i]]
        assigned[i] = True
        for j in range(i + 1, len(objects)):
            if assigned[j]:
                continue
            if are_isomorphic(hi, objects[j], method=method):  # type: ignore[arg-type]
                current.append(triples[j])
                assigned[j] = True
        classes.append(current)

    return classes


def analyze_characteristics(
    isomorphism_classes: list[list[tuple[int, int, int]]],
) -> tuple[list[float | None], list[float | None]]:
    """Characteristic and C-characteristic of one representative per class."""
    characteristics: list[float | None] = []
    c_characteristics: list[float | None] = []

    for iso_class in isomorphism_classes:
        p, k, d = iso_class[0]
        try:
            h = QuotientHyperfield.from_params(p, k, d)
            characteristics.append(h.char)
            c_characteristics.append(h.c_char)
        except Exception:
            characteristics.append(None)
            c_characteristics.append(None)

    return characteristics, c_characteristics
