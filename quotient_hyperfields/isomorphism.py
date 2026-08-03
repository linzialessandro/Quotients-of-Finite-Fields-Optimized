"""
Isomorphism of finite-field quotient hyperfields.

Layers
------
1. Baker–Jin (Theorem 1.1, arXiv:1912.11496) — O(1) for q >= N_r
   Default N_r is the sharp Remark 1.2 Weil bound (not the coarser r^4).
2. General — Aut(C_r) search matching 1 ⊞ x tables (gold standard / sporadics)
3. Policy wrapper — Baker–Jin when applicable, else general
"""

from __future__ import annotations

from dataclasses import dataclass
from enum import Enum
from functools import lru_cache
from math import gcd, sqrt
from typing import Literal

from quotient_hyperfields.hyperfield import ZERO, QuotientHyperfield
from quotient_hyperfields.primes import factor_prime_power


class IsoMethod(str, Enum):
    BAKER_JIN = "baker_jin"
    GENERAL = "general"
    AUTO = "auto"


@dataclass(frozen=True)
class IsoResult:
    """Result of an isomorphism query."""

    isomorphic: bool | None
    method: str
    reason: str
    sporadic: bool = False

    def __bool__(self) -> bool:
        return bool(self.isomorphic)


def weil_excess(N: float, r: int) -> float:
    """
    Lower-bound surplus from the Baker–Jin / Davenport–Hasse estimate:

        (N + 1) - (r - 1)(r - 2) √N - 3r

    The large-q regime of Theorem 1.1 applies when this quantity is positive
    (Remark 1.2).
    """
    if N < 0:
        raise ValueError(f"N must be non-negative, got {N}")
    if r < 2:
        raise ValueError(f"r must be >= 2 for the Weil excess, got {r}")
    return (N + 1) - (r - 1) * (r - 2) * sqrt(N) - 3 * r


def is_baker_jin_large(q: int, r: int) -> bool:
    """
    True iff field order q is in the large regime for index r (Remark 1.2).

    Equivalent to q >= baker_jin_bound(r) for integer q.
    """
    if r < 2:
        return True
    if q < 2:
        return False
    return weil_excess(q, r) > 0


@lru_cache(maxsize=None)
def baker_jin_bound_safe(r: int) -> int:
    """Coarse bound N_r = r^4 from the statement of Theorem 1.1."""
    if r < 2:
        return 0
    return r**4


@lru_cache(maxsize=None)
def baker_jin_bound_sharp(r: int) -> int:
    """
    Sharp Baker–Jin bound N_r from Remark 1.2.

    Smallest positive integer N such that

        (N + 1) - (r - 1)(r - 2) √N - 3r > 0.

    Paper checks: N_2 = 6, N_3 = 17. Always N_r <= r^4.
    """
    if r < 2:
        return 0
    if r == 2:
        return 6

    # f(N) = N + 1 - c√N - 3r with c = (r-1)(r-2).
    # f → ∞, single valley, so the first N with f(N) > 0 is a valid threshold.
    # Search in [1, r^4] by binary search (r^4 is always safe).
    hi = baker_jin_bound_safe(r)
    if weil_excess(hi, r) <= 0:
        # Should never happen for r >= 2; fall back to hi.
        return hi

    lo = 1
    while lo < hi:
        mid = (lo + hi) // 2
        if weil_excess(mid, r) > 0:
            hi = mid
        else:
            lo = mid + 1
    return lo


def baker_jin_bound(r: int) -> int:
    """
    Default N_r used by Baker–Jin iso checks: sharp Remark 1.2 bound.

    Use :func:`baker_jin_bound_safe` for the coarser r^4 threshold.
    """
    return baker_jin_bound_sharp(r)


def are_isomorphic_baker_jin(
    h1: QuotientHyperfield | tuple,
    h2: QuotientHyperfield | tuple,
    *,
    n_r: int | None = None,
) -> IsoResult:
    """
    Apply Baker–Jin Theorem 1.1 with Remark 1.2 bound by default.

    Returns
    -------
    IsoResult
        isomorphic is True/False when the theorem applies; None when the pair
        is sporadic (some q < N_r) and a general check is required.

    Parameters
    ----------
    n_r :
        Optional override for the large-q threshold. Default is
        :func:`baker_jin_bound` (Remark 1.2). Pass ``baker_jin_bound_safe(r)``
        to force the coarser r^4 cutoff.
    """
    a = _as_hyperfield(h1)
    b = _as_hyperfield(h2)

    if a.r != b.r:
        return IsoResult(
            False,
            "baker_jin",
            f"different indices r1={a.r}, r2={b.r}",
        )

    r = a.r
    if r < 2:
        # r = 1: Krasner K vs F_2 — both are quotients; compare orders / fields directly
        # F_q / F_q^x has 2 elements. For q=2 get F_2; for q>2 get Krasner K.
        # They are isomorphic iff both are K or both are F_2.
        both_f2 = a.q == 2 and b.q == 2
        both_k = a.q > 2 and b.q > 2
        return IsoResult(
            both_f2 or both_k,
            "baker_jin",
            "r=1 special case (Krasner K vs F_2)",
        )

    bound = n_r if n_r is not None else baker_jin_bound(r)
    # Prefer the analytic Remark 1.2 test when using the default bound:
    # q is large iff weil_excess(q, r) > 0. When n_r is overridden, use q >= n_r.
    if n_r is None:
        sporadic = not is_baker_jin_large(a.q, r) or not is_baker_jin_large(b.q, r)
        bound_desc = f"N_r={bound} (Remark 1.2)"
    else:
        sporadic = a.q < bound or b.q < bound
        bound_desc = f"N_r={bound}"

    if sporadic:
        return IsoResult(
            None,
            "baker_jin",
            f"sporadic: need general check ({bound_desc}, q1={a.q}, q2={b.q})",
            sporadic=True,
        )

    if r % 2 == 1:
        return IsoResult(
            True,
            "baker_jin",
            f"odd r={r}, both q >= {bound}",
        )

    # even r: class determined by q mod 2r
    if (a.q % (2 * r)) == (b.q % (2 * r)):
        return IsoResult(
            True,
            "baker_jin",
            f"even r={r}, same residue mod {2 * r}",
        )
    return IsoResult(
        False,
        "baker_jin",
        f"even r={r}, different residues mod {2 * r} "
        f"({a.q % (2 * r)} vs {b.q % (2 * r)})",
    )


def are_isomorphic_general(
    h1: QuotientHyperfield | tuple,
    h2: QuotientHyperfield | tuple,
) -> IsoResult:
    """
    Gold-standard isomorphism test for finite-field quotients.

    Both multiplicative groups are cyclic of order r. A hyperfield isomorphism is
    determined by a group automorphism phi([g]) = [g]^k (gcd(k, r) = 1) that
    satisfies phi(1 ⊞ x) = 1 ⊞ phi(x) for all x (cf. Ameri et al., Prop. 3.23;
    Baker–Jin structure of large quotients).
    """
    a = _as_hyperfield(h1)
    b = _as_hyperfield(h2)

    if a.r != b.r:
        return IsoResult(False, "general", f"different r ({a.r} vs {b.r})")

    r = a.r
    if r == 0:
        return IsoResult(True, "general", "trivial r=0")

    # Precompute 1 ⊞ x tables for both (labels in {-1} ∪ Z/rZ)
    domain = [ZERO] + list(range(r))
    table_a = {x: a.one_plus(x) for x in domain}
    table_b = {x: b.one_plus(x) for x in domain}

    # Try all automorphisms of C_r
    for k in range(r):
        if gcd(k, r) != 1 and r > 1:
            continue
        if r == 1:
            k = 0  # only identity on the trivial group; label 0 is [1]

        def phi(x: int, kk: int = k) -> int:
            if x == ZERO:
                return ZERO
            return (kk * x) % r if r > 0 else 0

        ok = True
        for x in domain:
            left = frozenset(phi(y) for y in table_a[x])
            right = table_b[phi(x)]
            if left != right:
                ok = False
                break
        if ok:
            return IsoResult(
                True,
                "general",
                f"matched under Aut multiplier k={k}",
            )
        if r == 1:
            break

    return IsoResult(False, "general", "no Aut(C_r) matches 1+x tables")


def are_isomorphic(
    h1: QuotientHyperfield | tuple,
    h2: QuotientHyperfield | tuple,
    method: Literal["auto", "baker_jin", "general"] = "auto",
    *,
    n_r: int | None = None,
) -> IsoResult:
    """
    Policy wrapper.

    - baker_jin: Theorem 1.1 only (isomorphic may be None if sporadic)
    - general: always the Aut-based table check
    - auto: Baker–Jin when decisive; general on sporadic or r < 2 edge handling
    """
    if method == "baker_jin":
        return are_isomorphic_baker_jin(h1, h2, n_r=n_r)
    if method == "general":
        return are_isomorphic_general(h1, h2)
    if method != "auto":
        raise ValueError(f"Unknown method {method!r}")

    bj = are_isomorphic_baker_jin(h1, h2, n_r=n_r)
    if bj.sporadic or bj.isomorphic is None:
        gen = are_isomorphic_general(h1, h2)
        return IsoResult(
            gen.isomorphic,
            "auto/general",
            f"{bj.reason}; {gen.reason}",
            sporadic=True,
        )
    return IsoResult(
        bj.isomorphic,
        "auto/baker_jin",
        bj.reason,
        sporadic=False,
    )


def _as_hyperfield(
    obj: QuotientHyperfield | tuple,
) -> QuotientHyperfield:
    if isinstance(obj, QuotientHyperfield):
        return obj
    if isinstance(obj, tuple):
        if len(obj) == 2:
            q, d = obj
            return QuotientHyperfield.from_q_d(int(q), int(d))
        if len(obj) == 3:
            p, k, d = obj
            return QuotientHyperfield.from_params(int(p), int(k), int(d))
    raise TypeError(
        "Expected QuotientHyperfield or (q, d) or (p, k, d), "
        f"got {type(obj).__name__}"
    )


# Back-compat helpers used by census / CLI ---------------------------------

def are_isomorphic_optimized(q1: int, d1: int, q2: int, d2: int) -> bool:
    """Legacy boolean API: auto policy, True/False only."""
    result = are_isomorphic((q1, d1), (q2, d2), method="auto")
    return bool(result.isomorphic)
