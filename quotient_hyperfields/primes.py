"""Primality and prime-power factorisation helpers."""

from __future__ import annotations

import math
from functools import lru_cache


def is_prime(n: int) -> bool:
    """Trial-division primality test (adequate for laptop-scale parameters)."""
    if n < 2:
        return False
    if n == 2:
        return True
    if n % 2 == 0:
        return False
    limit = int(math.isqrt(n))
    for d in range(3, limit + 1, 2):
        if n % d == 0:
            return False
    return True


def factor_prime_power(q: int) -> tuple[int, int]:
    """
    Factor q = p^k with p prime and k >= 1.

    Raises
    ------
    ValueError
        If q < 2 or q is not a prime power.
    """
    if q < 2:
        raise ValueError(f"Expected prime power >= 2, got {q}")
    if is_prime(q):
        return q, 1

    # Smallest prime factor
    p = None
    if q % 2 == 0:
        p = 2
    else:
        limit = int(math.isqrt(q))
        for d in range(3, limit + 1, 2):
            if q % d == 0:
                p = d
                break
    if p is None:
        # q is prime (isqrt edge cases)
        return q, 1

    k = 0
    n = q
    while n % p == 0:
        n //= p
        k += 1
    if n != 1:
        raise ValueError(f"{q} is not a prime power")
    if p**k != q:
        raise ValueError(f"{q} is not a prime power")
    return p, k


def gcd(a: int, b: int) -> int:
    while b:
        a, b = b, a % b
    return abs(a)


def prime_factors(n: int) -> list[int]:
    """Distinct prime factors of n (n >= 1)."""
    if n < 1:
        raise ValueError(f"n must be positive, got {n}")
    if n == 1:
        return []
    factors: list[int] = []
    d = 2
    while d * d <= n:
        if n % d == 0:
            factors.append(d)
            while n % d == 0:
                n //= d
        d += 1
    if n > 1:
        factors.append(n)
    return factors


@lru_cache(maxsize=None)
def primitive_root_mod_p(p: int) -> int:
    """
    Smallest primitive root modulo a prime p.

    Used by the prime-field backend of :class:`QuotientHyperfield` to avoid
    constructing a ``galois`` field for every prime order.
    """
    if not is_prime(p):
        raise ValueError(f"p must be prime, got {p}")
    if p == 2:
        return 1
    if p == 3:
        return 2
    factors = prime_factors(p - 1)
    for g in range(2, p):
        if all(pow(g, (p - 1) // q, p) != 1 for q in factors):
            return g
    raise RuntimeError(f"no primitive root found modulo {p}")
