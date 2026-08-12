"""
Core object: Krasner quotient hyperfield F_q / G_d.

Nonzero elements are labeled by discrete logs in Z/rZ, where r = (q-1)/d is the
index of G_d in F_q^x. Zero is the sentinel label -1.

Backends (v0.2+)
----------------
- **prime** (``k = 1``): pure modular integer arithmetic; no ``galois`` field
  object. This is the fast path for the bulk of paper scans.
- **extension** (``k > 1``): ``galois`` finite field, with cosets precomputed
  so hyperaddition is only |G|^2 field additions per pair.

Mathematical definitions (hyperaddition, char / C-char) are unchanged from
v0.1.0 / arXiv:2608.03625. Large-``q`` char shortcuts apply Baker–Jin + the
stable-characteristic theorem when ``use_stable=True`` (default).
"""

from __future__ import annotations

from functools import cached_property
from math import inf
from typing import Iterable, Literal

import pandas as pd

from quotient_hyperfields.primes import (
    factor_prime_power,
    is_prime,
    primitive_root_mod_p,
)


# Sentinel for the zero element of the hyperfield
ZERO = -1

Backend = Literal["prime", "extension"]


class QuotientHyperfield:
    """
    Krasner quotient hyperfield K = F_q / G, with |K| = r + 1.

    Parameters may be given as (p, k, d) or as (q, d) / (q, r).
    """

    def __init__(self, p: int, k: int, d: int):
        if not is_prime(p):
            raise ValueError(f"p must be prime, got {p}")
        if k < 1:
            raise ValueError(f"k must be a positive integer, got {k}")
        if d < 1:
            raise ValueError(f"d must be a positive integer, got {d}")

        q = p**k
        if (q - 1) % d != 0:
            raise ValueError(f"d={d} must divide q-1={q - 1}")

        self.p = p
        self.k = k
        self.q = q
        self.d = d
        self.r = (q - 1) // d  # index of G; |K^x| = r; |K| = r + 1

        # Cache of hyperaddition results: (i, j) -> frozenset of labels
        self._hyperadd_cache: dict[tuple[int, int], frozenset[int]] = {}
        # Cache of S_n = n x [1]
        self._n_ones_cache: dict[int, frozenset[int]] = {}
        self._fingerprint_cache: tuple | None = None

        if k == 1:
            self._init_prime_backend()
        else:
            self._init_extension_backend()

    # ------------------------------------------------------------------
    # Backends
    # ------------------------------------------------------------------

    def _init_prime_backend(self) -> None:
        """F_p with modular integer arithmetic (no galois)."""
        p = self.p
        r = self.r
        d = self.d
        alpha = primitive_root_mod_p(p)
        self.backend: Backend = "prime"
        self._GF = None
        self._primitive = alpha

        # G = <alpha^r> has order d
        h = pow(alpha, r, p)
        subgroup = [pow(h, i, p) for i in range(d)]
        self._subgroup_list: list[int] = subgroup
        self._subgroup_keys: frozenset[int] = frozenset(subgroup)

        # Coset of alpha^i is labeled i (mod r); representative alpha^i
        reps = [pow(alpha, i, p) for i in range(r)]
        self._rep_keys: list[int] = reps

        # label[field_int] for field_int in 1..p-1
        label = [0] * p
        cosets: list[list[int]] = []
        for i, rep in enumerate(reps):
            cells = [(rep * g) % p for g in subgroup]
            cosets.append(cells)
            for x in cells:
                label[x] = i
        self._element_to_label_list: list[int] = label
        self._element_to_label: dict[int, int] = {
            x: label[x] for x in range(1, p)
        }
        self._cosets_int: list[list[int]] = cosets

    def _init_extension_backend(self) -> None:
        """F_{p^k} via galois, with precomputed cosets for faster hyperadd."""
        import galois

        q = self.q
        r = self.r
        d = self.d
        self.backend = "extension"
        GF = galois.GF(q)
        self._GF = GF
        alpha = GF.primitive_element
        self._primitive = int(alpha)

        h = alpha**r
        subgroup = [int(h**i) for i in range(d)]
        self._subgroup_list = subgroup
        self._subgroup_keys = frozenset(subgroup)

        reps = [int(alpha**i) for i in range(r)]
        self._rep_keys = reps

        element_to_label: dict[int, int] = {}
        cosets_int: list[list[int]] = []
        cosets_gf: list[list] = []
        for i, rep in enumerate(reps):
            cells_int: list[int] = []
            cells_gf = []
            rep_el = GF(rep)
            for gk in subgroup:
                val = rep_el * GF(gk)
                key = int(val)
                cells_int.append(key)
                cells_gf.append(val)
                element_to_label[key] = i
            cosets_int.append(cells_int)
            cosets_gf.append(cells_gf)

        self._element_to_label = element_to_label
        self._element_to_label_list = None
        self._cosets_int = cosets_int
        self._cosets_gf = cosets_gf

    # ------------------------------------------------------------------
    # Constructors
    # ------------------------------------------------------------------

    @classmethod
    def from_params(cls, p: int, k: int, d: int) -> QuotientHyperfield:
        """Build from prime p, extension degree k, and subgroup order d."""
        return cls(p, k, d)

    @classmethod
    def from_q_d(cls, q: int, d: int) -> QuotientHyperfield:
        """Build from field order q and subgroup order d."""
        p, k = factor_prime_power(q)
        return cls(p, k, d)

    @classmethod
    def from_q_r(cls, q: int, r: int) -> QuotientHyperfield:
        """Build from field order q and subgroup index r (so |K| = r + 1)."""
        if r < 1:
            raise ValueError(f"r must be a positive integer, got {r}")
        if (q - 1) % r != 0:
            raise ValueError(f"r={r} must divide q-1={q - 1}")
        d = (q - 1) // r
        return cls.from_q_d(q, d)

    # ------------------------------------------------------------------
    # Basic properties
    # ------------------------------------------------------------------

    @property
    def order(self) -> int:
        """Number of elements of the hyperfield (|K| = r + 1)."""
        return self.r + 1

    @property
    def index(self) -> int:
        """Index r = [F_q^x : G] = |K^x|."""
        return self.r

    def __repr__(self) -> str:
        return (
            f"QuotientHyperfield(q={self.q}, d={self.d}, r={self.r}, "
            f"backend={self.backend!r})"
        )

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, QuotientHyperfield):
            return NotImplemented
        return self.q == other.q and self.d == other.d

    def __hash__(self) -> int:
        return hash((self.q, self.d))

    # ------------------------------------------------------------------
    # Hyperaddition
    # ------------------------------------------------------------------

    def hyperadd(self, i: int, j: int) -> frozenset[int]:
        """
        Hyperaddition of two elements labeled in {-1} ∪ {0,...,r-1}.

        -1 denotes zero. For nonzero labels, addition is the Krasner coset sum
        (same definition as v0.1.0 / arXiv:2608.03625).
        """
        if i == ZERO and j == ZERO:
            return frozenset({ZERO})
        if i == ZERO:
            return frozenset({j})
        if j == ZERO:
            return frozenset({i})

        if not (0 <= i < self.r and 0 <= j < self.r):
            raise ValueError(f"Labels must be -1 or in 0..{self.r - 1}, got {i}, {j}")

        key = (i, j) if i <= j else (j, i)
        cached = self._hyperadd_cache.get(key)
        if cached is not None:
            return cached

        if self.backend == "prime":
            out = self._hyperadd_prime(i, j)
        else:
            out = self._hyperadd_extension(i, j)

        self._hyperadd_cache[key] = out
        return out

    def _hyperadd_prime(self, i: int, j: int) -> frozenset[int]:
        p = self.p
        label = self._element_to_label_list
        assert label is not None
        result: set[int] = set()
        for a in self._cosets_int[i]:
            for b in self._cosets_int[j]:
                s = a + b
                if s >= p:
                    s -= p
                if s == 0:
                    result.add(ZERO)
                else:
                    result.add(label[s])
        return frozenset(result)

    def _hyperadd_extension(self, i: int, j: int) -> frozenset[int]:
        result: set[int] = set()
        etol = self._element_to_label
        for a in self._cosets_gf[i]:
            for b in self._cosets_gf[j]:
                s = a + b
                if s == 0:
                    result.add(ZERO)
                else:
                    result.add(etol[int(s)])
        return frozenset(result)

    def one_plus(self, x: int) -> frozenset[int]:
        """Return [1] ⊞ x. Label of [1] is 0."""
        return self.hyperadd(0, x)

    def n_times_one(self, n: int) -> frozenset[int]:
        """
        Hyper-sum of n copies of [1], written n ×_K [1].

        Convention: 1 × [1] = {[1]} = {0} as a label set.
        """
        if n < 1:
            raise ValueError("n must be >= 1")
        if n in self._n_ones_cache:
            return self._n_ones_cache[n]

        if n == 1:
            out = frozenset({0})
            self._n_ones_cache[1] = out
            return out

        prev = self.n_times_one(n - 1)
        result: set[int] = set()
        for a in prev:
            result |= self.hyperadd(a, 0)  # a ⊞ [1]
        out = frozenset(result)
        self._n_ones_cache[n] = out
        return out

    # ------------------------------------------------------------------
    # Invariants: characteristic and C-characteristic
    # (Kędzierski–Linzi–Stojałowska definitions)
    # ------------------------------------------------------------------

    def _stable_char_pair(self) -> tuple[float, float] | None:
        """
        (char, C-char) from the stable-class theorem when q is large.

        Returns None if the Baker–Jin large regime does not apply (or r < 2),
        so the iterative definition must be used.
        """
        # Local import avoids circular dependency at module load.
        from quotient_hyperfields.isomorphism import is_baker_jin_large

        r = self.r
        if r < 2:
            return None
        if not is_baker_jin_large(self.q, r):
            return None
        # Stable-characteristic theorem (arXiv:2608.03625):
        # odd r → (2, 1); even r, q ≡ 1 (mod 2r) → (2, 1);
        # even r, q ≡ r+1 (mod 2r) → (3, 1).
        if r % 2 == 1:
            return 2.0, 1.0
        mod = self.q % (2 * r)
        if mod == 1:
            return 2.0, 1.0
        if mod == r + 1:
            return 3.0, 1.0
        return None

    def characteristic(
        self,
        max_terms: int | None = None,
        *,
        use_stable: bool = True,
    ) -> float:
        """
        Minimal n >= 1 such that 0 ∈ n × [1], or +∞ if none exists up to bound.

        For finite quotients the characteristic is always finite and at most q
        (in practice much smaller). Default search bound is max(2000, 2*q).

        Parameters
        ----------
        use_stable :
            If True (default) and q is in the Baker–Jin large regime, return
            the value from the stable-characteristic theorem without iterating.
            Pass False to force the definitional loop (e.g. for verification).
        """
        if use_stable and max_terms is None:
            pair = self._stable_char_pair()
            if pair is not None:
                return pair[0]

        bound = max_terms if max_terms is not None else max(2000, 2 * self.q)
        for n in range(1, bound + 1):
            if ZERO in self.n_times_one(n):
                return float(n)
        return inf

    def c_characteristic(
        self,
        max_terms: int | None = None,
        *,
        use_stable: bool = True,
    ) -> float:
        """
        Minimal n >= 1 such that 1 ∈ (n+1) × [1], or +∞ if none exists up to bound.

        Always C-char <= char when both are finite.

        Parameters
        ----------
        use_stable :
            See :meth:`characteristic`.
        """
        if use_stable and max_terms is None:
            pair = self._stable_char_pair()
            if pair is not None:
                return pair[1]

        bound = max_terms if max_terms is not None else max(2000, 2 * self.q)
        for n in range(1, bound + 1):
            if 0 in self.n_times_one(n + 1):  # label 0 is [1]
                return float(n)
        return inf

    @cached_property
    def char(self) -> float:
        return self.characteristic()

    @cached_property
    def c_char(self) -> float:
        return self.c_characteristic()

    # ------------------------------------------------------------------
    # Addition table
    # ------------------------------------------------------------------

    def addition_table(self, as_dataframe: bool = True):
        """
        Full hyperaddition table including zero.

        Labels: '0' for zero, and integer strings for discrete-log labels 0..r-1
        (where '0' as a nonzero label is written '1' meaning [1] = g^0).
        """
        labels = [ZERO] + list(range(self.r))
        table = [[self.hyperadd(i, j) for j in labels] for i in labels]

        def fmt_label(x: int) -> str:
            if x == ZERO:
                return "0"
            if x == 0:
                return "1"
            return f"g^{x}"

        def fmt_cell(s: Iterable[int]) -> str:
            parts = sorted(s, key=lambda t: (-1 if t == ZERO else t))
            return "{" + ",".join(f"[{fmt_label(t)}]" for t in parts) + "}"

        headers = [fmt_label(x) for x in labels]
        formatted = [[fmt_cell(cell) for cell in row] for row in table]

        if not as_dataframe:
            return formatted

        return pd.DataFrame(
            formatted,
            index=[f"[{h}]" for h in headers],
            columns=[f"[{h}]" for h in headers],
        )

    def one_plus_table_masks(self) -> dict[int, int]:
        """
        ``1 ⊞ x`` for all labels ``x``, encoded as bitmasks (Talotti-style).

        See :mod:`quotient_hyperfields.bitsets`.
        """
        from quotient_hyperfields.bitsets import labels_to_mask

        domain = [ZERO] + list(range(self.r))
        return {x: labels_to_mask(self.one_plus(x)) for x in domain}

    # ------------------------------------------------------------------
    # Baker–Jin residue class (for even r)
    # ------------------------------------------------------------------

    @property
    def baker_jin_class(self) -> str | None:
        """
        Iso class label under Baker–Jin for large q.

        - odd r: 'odd' (unique class when q >= N_r)
        - even r: 'even_0' if q ≡ 1 (mod 2r), 'even_1' if q ≡ r+1 (mod 2r)
        - None if parameters are invalid for the theorem statement (r < 2)
        """
        if self.r < 2:
            return "r1"
        if self.r % 2 == 1:
            return "odd"
        mod = self.q % (2 * self.r)
        if mod == 1:
            return "even_0"
        if mod == self.r + 1:
            return "even_1"
        # Should not happen when r | (q-1)
        return f"even_mod_{mod}"
