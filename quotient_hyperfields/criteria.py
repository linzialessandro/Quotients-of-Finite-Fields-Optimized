"""
Paper criteria that finite-field quotient hyperfields must satisfy.

References
----------
- Massouros–Massouros, AIMS Mathematics 10 (2025), 21287–21421,
  Propositions on sum cardinalities in a quotient F/Q
  (cf. DOI 10.3934/math.2025951).
- Kędzierski–Linzi–Stojałowska, Mathematics 2023, Proposition 7,
  Lemma 3, Corollary 1 (characteristic bounds from |G|).

These checks are used as sanity tests on :class:`QuotientHyperfield` and as
experiment harnesses in :mod:`papers_experiments`.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field

from quotient_hyperfields.hyperfield import ZERO, QuotientHyperfield
from quotient_hyperfields.primes import prime_factors


@dataclass
class CheckResult:
    """Outcome of a single paper criterion check."""

    name: str
    passed: bool
    detail: str
    failures: list[str] = field(default_factory=list)

    def __bool__(self) -> bool:
        return self.passed


# ---------------------------------------------------------------------------
# Massouros–Massouros (order-7 paper): sum cardinality
# ---------------------------------------------------------------------------


def massouros_prop1_sum_bound(h: QuotientHyperfield) -> CheckResult:
    """
    Massouros Proposition 1.

    In a quotient hyperfield F/Q, the number of summands in x+y cannot exceed
    |Q|. Here Q = G_d has order ``h.d``, and x+y is the hyperaddition of cosets
    (as a set of hyperfield elements, including 0).
    """
    d = h.d
    labels = [ZERO] + list(range(h.r))
    failures: list[str] = []
    for i in labels:
        for j in labels:
            if j < i:
                continue  # commutative
            s = h.hyperadd(i, j)
            if len(s) > d:
                failures.append(
                    f"|[{i}]+[{j}]|={len(s)} > d={d}"
                )
    return CheckResult(
        name="Massouros Prop.1 (|x+y| ≤ |G|)",
        passed=len(failures) == 0,
        detail=f"checked all pairs; d=|G|={d}, order={h.order}",
        failures=failures[:20],
    )


def _difference_set(h: QuotientHyperfield, x: int) -> frozenset[int]:
    """Cosets appearing in x − x = x + (−x)."""
    # −[g^i] = [g^i · (−1)] ; −1 = alpha^{(q-1)/2} has discrete log (q-1)/2 mod r
    # Label of −1: (q-1)/2 mod r = (d*r)/2 mod r. If r odd, (q-1)/2 ≡ 0 mod r
    # when −1 ∈ G (r odd ⇒ −1 in G for the unique index-r subgroup? Actually
    # −1 = alpha^{(q-1)/2}, coset log is ((q-1)/2) mod r = (d r / 2) mod r.
    # If d even, (d/2)*r ≡ 0 mod r so −1 ∈ G, label of −1 is 0? Wait coset of -1.
    # coset label of alpha^k is k mod r. k=(q-1)/2 = d*r/2.
    # d*r/2 mod r = 0 if d even, else (r/2) if d odd and r even.
    if h.q % 2 == 0:
        # char 2 field: −1 = 1, so −x = x
        return h.hyperadd(x, x) if x != ZERO else frozenset({ZERO})
    # odd characteristic
    half = (h.q - 1) // 2
    neg_one_label = half % h.r  # label of [−1]
    if x == ZERO:
        return frozenset({ZERO})
    neg_x = (x + neg_one_label) % h.r
    return h.hyperadd(x, neg_x)


def massouros_prop2_when_applicable(h: QuotientHyperfield) -> CheckResult:
    """
    Massouros Proposition 2 (conditional).

    If for every nonzero coset x the differences x−x share only 0 pairwise in the
    strong form of the hypothesis, then non-opposite unequal pairs have
    |x+y| = |G|.

    We check a practical version used in their order-7 arguments:
    - Always re-verify Prop.1.
    - For pairs of distinct nonzero non-opposite labels, if the Prop.2 hypothesis
      holds globally (x−x ∩ y−y = {0} for all distinct nonzero x,y), then
      |x+y| must equal d.
    """
    d = h.d
    nonzero = list(range(h.r))
    # Precompute x − x for each nonzero label
    diff = {x: _difference_set(h, x) for x in nonzero}
    # Hypothesis: pairwise intersections of (x−x)\\{0} empty? Prop says
    # (xQ−xQ) ∩ (yQ−yQ) = {0} for the relevant pairs — i.e. only zero in common.
    hyp_ok = True
    for i, x in enumerate(nonzero):
        for y in nonzero[i + 1 :]:
            inter = diff[x] & diff[y]
            if inter - {ZERO}:
                hyp_ok = False
                break
        if not hyp_ok:
            break

    if not hyp_ok:
        return CheckResult(
            name="Massouros Prop.2 (conditional)",
            passed=True,
            detail="hypothesis (x−x)∩(y−y)={0} fails globally; Prop.2 not applicable",
        )

    # −1 label for opposite test
    if h.q % 2 == 0:
        opp = {x: x for x in nonzero}  # char 2
    else:
        neg_one = ((h.q - 1) // 2) % h.r
        opp = {x: (x + neg_one) % h.r for x in nonzero}

    failures: list[str] = []
    for i, x in enumerate(nonzero):
        for y in nonzero:
            if x == y:
                continue
            if opp[x] == y:
                continue  # opposite
            s = h.hyperadd(x, y)
            if len(s) != d:
                failures.append(f"|[{x}]+[{y}]|={len(s)} ≠ d={d}")
    return CheckResult(
        name="Massouros Prop.2 (when hyp. holds)",
        passed=len(failures) == 0,
        detail="hypothesis holds; checked non-opposite unequal pairs have |sum|=d",
        failures=failures[:20],
    )


def verify_massouros_criteria(h: QuotientHyperfield) -> list[CheckResult]:
    """Run Massouros Propositions 1–2 checks."""
    return [massouros_prop1_sum_bound(h), massouros_prop2_when_applicable(h)]


# ---------------------------------------------------------------------------
# Linzi et al.: characteristic bounds from |G|
# ---------------------------------------------------------------------------


def linzi_prop7_char_bound(h: QuotientHyperfield) -> CheckResult:
    """
    Proposition 7 (char paper).

    If n > 1 divides |G| = d, then char(K_G) ≤ n.
    Equivalently: char ≤ every divisor n>1 of d, so char ≤ smallest prime factor
    of d when d > 1. If d = 1, G = {1} and the quotient is the field; char = p.
    """
    d = h.d
    ch = h.char
    if math.isinf(ch):
        return CheckResult(
            name="Linzi Prop.7 (char ≤ n for n|d, n>1)",
            passed=False,
            detail="characteristic is infinite (unexpected for finite-field quotients)",
        )
    if d == 1:
        # field case: char should equal p
        ok = ch == float(h.p)
        return CheckResult(
            name="Linzi Prop.7 (d=1 field case)",
            passed=ok,
            detail=f"d=1 ⇒ field F_{h.q}; char={ch}, p={h.p}",
        )
    # Every n|d with n>1 gives char ≤ n. Tightest: min such n = least prime factor.
    factors = prime_factors(d)
    if not factors:
        return CheckResult(
            name="Linzi Prop.7",
            passed=False,
            detail=f"no prime factors for d={d}",
        )
    bound = min(factors)
    # Also char ≤ d (since d>1 divides d)
    bound = min(bound, d)
    # Verify all n>1 dividing d
    failures = []
    # check against all divisors > 1 up to sqrt for thoroughness
    n = 2
    while n * n <= d:
        if d % n == 0:
            if ch > n:
                failures.append(f"char={ch} > n={n} | d={d}")
            m = d // n
            if m > 1 and ch > m:
                failures.append(f"char={ch} > n={m} | d={d}")
        n += 1
    if ch > d:
        failures.append(f"char={ch} > d={d}")
    return CheckResult(
        name="Linzi Prop.7 (char ≤ n for n|d, n>1)",
        passed=len(failures) == 0,
        detail=f"char={ch}, d={d}, least prime factor bound={bound}",
        failures=failures[:20],
    )


def linzi_corollary1_char2_iff_even_d(h: QuotientHyperfield) -> CheckResult:
    """
    Corollary 1 (char paper).

    If char(K) ≠ 2 (i.e. underlying field prime p odd, so 1 ≠ −1), then
    char(K_G) = 2  ⇔  |G| is even.
    """
    if h.p == 2:
        return CheckResult(
            name="Linzi Cor.1 (char=2 ⇔ |G| even)",
            passed=True,
            detail="p=2: corollary assumes 1≠−1 in K; skipped",
        )
    ch = h.char
    d_even = h.d % 2 == 0
    if math.isinf(ch):
        return CheckResult(
            name="Linzi Cor.1 (char=2 ⇔ |G| even)",
            passed=False,
            detail="infinite char",
        )
    char_is_2 = ch == 2.0
    ok = char_is_2 == d_even
    return CheckResult(
        name="Linzi Cor.1 (char=2 ⇔ |G| even)",
        passed=ok,
        detail=f"p={h.p} odd, char={ch}, d={h.d} even={d_even}",
        failures=[]
        if ok
        else [f"char_is_2={char_is_2} but d_even={d_even}"],
    )


def linzi_prop9_cchar_bound(h: QuotientHyperfield) -> CheckResult:
    """
    Proposition 9 (char paper).

    If n > 1 divides |G|, then C-char(K_G) ≤ n (same bound as Prop.7 for C-char).
    """
    d = h.d
    cc = h.c_char
    if math.isinf(cc):
        return CheckResult(
            name="Linzi Prop.9 (C-char ≤ n for n|d)",
            passed=False,
            detail="C-char infinite",
        )
    if d == 1:
        ok = cc == float(h.p)
        return CheckResult(
            name="Linzi Prop.9 (d=1)",
            passed=ok,
            detail=f"C-char={cc}, p={h.p}",
        )
    failures = []
    n = 2
    while n * n <= d:
        if d % n == 0:
            if cc > n:
                failures.append(f"C-char={cc} > n={n}")
            m = d // n
            if m > 1 and cc > m:
                failures.append(f"C-char={cc} > n={m}")
        n += 1
    if cc > d:
        failures.append(f"C-char={cc} > d={d}")
    return CheckResult(
        name="Linzi Prop.9 (C-char ≤ n for n|d)",
        passed=len(failures) == 0,
        detail=f"C-char={cc}, d={d}",
        failures=failures[:20],
    )


def verify_linzi_criteria(h: QuotientHyperfield) -> list[CheckResult]:
    """Run Linzi et al. Prop.7, Cor.1, Prop.9 checks."""
    return [
        linzi_prop7_char_bound(h),
        linzi_corollary1_char2_iff_even_d(h),
        linzi_prop9_cchar_bound(h),
    ]


def verify_all_paper_criteria(h: QuotientHyperfield) -> list[CheckResult]:
    """Massouros + Linzi criteria."""
    return verify_massouros_criteria(h) + verify_linzi_criteria(h)
