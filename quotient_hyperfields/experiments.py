"""
Computational probes for Baker–Jin open questions (Proc. AMS 149 (2021), §5).

Open questions (paraphrased)
----------------------------
(1) Classify all hyperfields of order 5 (quotients of finite / infinite fields).
(2) Algorithm: is a given finite hyperfield a quotient of some *infinite* field?
(3) True growth of N_r (minimal large-q threshold in Theorem 1.1).
(4) Growth of Q_r / H_r (quotient classes vs all hyperfield classes).

This module focuses on what finite-field quotients can answer directly:
  - (3) empirical N_r by fingerprinting F_q / G_r^q
  - (4) Q_r^{fin}: number of iso classes among *finite-field* quotients of order r+1
    (a lower bound on the paper's Q_r, which also allows infinite-field quotients)

Questions (1)–(2) need abstract hyperfield enumeration and are out of scope of the
finite-field quotient core; see module docstrings for status.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from functools import lru_cache

from quotient_hyperfields.hyperfield import ZERO, QuotientHyperfield
from quotient_hyperfields.isomorphism import (
    baker_jin_bound_safe,
    baker_jin_bound_sharp,
    is_baker_jin_large,
    weil_excess,
)
from quotient_hyperfields.primes import factor_prime_power, gcd, is_prime


# ---------------------------------------------------------------------------
# Prime powers
# ---------------------------------------------------------------------------


def is_prime_power(n: int) -> bool:
    try:
        factor_prime_power(n)
        return True
    except ValueError:
        return False


def prime_powers_up_to(limit: int) -> list[int]:
    """All prime powers q with 2 <= q <= limit."""
    if limit < 2:
        return []
    out: list[int] = []
    # primes
    for p in range(2, limit + 1):
        if not is_prime(p):
            continue
        pk = p
        while pk <= limit:
            out.append(pk)
            if pk > limit // p:
                break
            pk *= p
    return sorted(set(out))


def prime_powers_with_index(r: int, limit: int) -> list[int]:
    """Prime powers q <= limit with r | (q - 1)."""
    if r < 1:
        raise ValueError("r must be positive")
    return [q for q in prime_powers_up_to(limit) if (q - 1) % r == 0]


# ---------------------------------------------------------------------------
# Fingerprints (Aut-canonical 1⊞x table)
# ---------------------------------------------------------------------------


def structure_fingerprint(h: QuotientHyperfield) -> tuple:
    """
    Isomorphism invariant for finite-field quotients of fixed r.

    Canonical form of the multiset of 1⊞x sets under Aut(C_r): choose the
    lexicographically minimal encoding among multipliers k coprime to r.
    """
    r = h.r
    domain = [ZERO] + list(range(r))
    raw = {x: h.one_plus(x) for x in domain}

    def encode(k: int) -> tuple:
        def phi(x: int) -> int:
            if x == ZERO:
                return ZERO
            return (k * x) % r

        # Map each domain element to sorted image of 1⊞x under phi, keyed by phi(x)
        rows = []
        for x in domain:
            key = phi(x)
            vals = tuple(sorted(phi(y) for y in raw[x]))
            rows.append((key, vals))
        rows.sort()
        return tuple(rows)

    if r == 1:
        return encode(0)

    candidates = [encode(k) for k in range(r) if gcd(k, r) == 1]
    return min(candidates)


def baker_jin_residue(q: int, r: int) -> str:
    """Predicted stable class label for large q (Theorem 1.1)."""
    if r < 2:
        return "r1"
    if r % 2 == 1:
        return "odd"
    mod = q % (2 * r)
    if mod == 1:
        return "even_0"
    if mod == r + 1:
        return "even_1"
    return f"even_other_{mod}"


# ---------------------------------------------------------------------------
# Open question (3): empirical N_r
# ---------------------------------------------------------------------------


@dataclass
class EmpiricalNrResult:
    """Probe of the true Baker–Jin threshold for a fixed index r."""

    r: int
    max_q: int
    n_r_remark_12: int
    n_r_safe: int
    n_r_lower_paper: int | None
    n_r_empirical: int | None
    prime_powers: list[int]
    fingerprints: dict[int, tuple] = field(repr=False)
    stable_by_residue: dict[str, tuple] = field(repr=False)
    exceptions: list[tuple[int, str]]  # (q, residue) where fp ≠ stable
    residue_classes: dict[str, list[int]]
    verified_up_to: int
    notes: list[str] = field(default_factory=list)

    def summary(self) -> str:
        lines = [
            f"r = {self.r}",
            f"  Remark 1.2 N_r     = {self.n_r_remark_12}",
            f"  Safe r^4 N_r       = {self.n_r_safe}",
            f"  Paper lower bound  = {self.n_r_lower_paper}",
            f"  Empirical N_r      = {self.n_r_empirical}",
            f"  # prime powers ≤ {self.max_q} with r | (q-1): {len(self.prime_powers)}",
            f"  # exceptions (non-stable fps): {len(self.exceptions)}",
            f"  verified up to     = {self.verified_up_to}",
        ]
        if self.exceptions:
            ex = ", ".join(f"q={q}({res})" for q, res in self.exceptions[:12])
            if len(self.exceptions) > 12:
                ex += ", ..."
            lines.append(f"  exception qs       = {ex}")
        for note in self.notes:
            lines.append(f"  note: {note}")
        return "\n".join(lines)


def paper_lower_bound_nr(r: int) -> int | None:
    """
    Baker–Jin §5: N_r > (r-1)^2 whenever r-1 is prime
    (via H = F_{p^2} / F_p^x with p = r-1).

    Returns the integer lower bound L such that N_r >= L, i.e. (r-1)^2 + 1
    when r-1 is prime, else None if that construction does not apply.
    """
    p = r - 1
    if p >= 2 and is_prime(p):
        return p * p + 1  # N_r > p^2
    return None


def required_stable_residues(r: int) -> list[str]:
    """Baker–Jin stable class labels that must be witnessed for large q."""
    if r % 2 == 1:
        return ["odd"]
    return ["even_0", "even_1"]


def has_large_witnesses(r: int, qs: list[int]) -> bool:
    """True if every stable residue has at least one large prime power in qs."""
    n12 = baker_jin_bound_sharp(r)
    large = [q for q in qs if q >= n12]
    by_res: dict[str, list[int]] = {}
    for q in large:
        by_res.setdefault(baker_jin_residue(q, r), []).append(q)
    return all(res in by_res for res in required_stable_residues(r))


def resolve_scan_limit(
    r: int,
    max_q: int | None = None,
    *,
    max_q_cap: int = 8000,
) -> int:
    """
    Choose a scan ceiling that includes Remark 1.2 and large stable witnesses.

    Starts at max(max_q or N_r, N_r) and grows (up to max_q_cap) until each
    Baker–Jin stable residue has a prime power q >= N_r with r | (q-1).
    """
    n12 = baker_jin_bound_sharp(r)
    limit = n12 if max_q is None else max(max_q, n12)
    limit = min(limit, max_q_cap)

    while limit <= max_q_cap:
        qs = prime_powers_with_index(r, limit)
        if has_large_witnesses(r, qs):
            return limit
        # grow: need room for the next prime powers past N_r
        nxt = min(max_q_cap, max(limit + max(50, 10 * r), int(limit * 1.5) + 1))
        if nxt == limit:
            break
        limit = nxt
    return limit


def empirical_nr(
    r: int,
    max_q: int | None = None,
    *,
    max_q_cap: int = 8000,
    progress: bool = False,
) -> EmpiricalNrResult:
    """
    Estimate the true N_r for Theorem 1.1 by scanning prime powers.

    Method
    ------
    1. Resolve a scan limit that includes large stable witnesses (q >= N_r^{1.2}).
    2. Fingerprint each F_q / G_r (Aut-canonical 1⊞x table).
    3. Stable type per residue = fingerprint of the largest large q in that class
       (justified by Theorem 1.1 once a large witness exists).
    4. Exceptions: prime powers whose fingerprint differs from their class stable type.
    5. Empirical N_r = 1 + max(exception q). If there are no exceptions,
       Remark 1.2 is not sharp on this range and we report an upper bound equal
       to the smallest large prime power (true N_r is at most that).

    Notes
    -----
    This is a computational lower/upper sandwich relative to the scan, not a
    proof of the exact minimal N_r for all larger q (that still rests on Thm 1.1
    for q past Remark 1.2).
    """
    if r < 2:
        raise ValueError("empirical N_r is for r >= 2")

    n12 = baker_jin_bound_sharp(r)
    n_safe = baker_jin_bound_safe(r)
    max_q = resolve_scan_limit(r, max_q, max_q_cap=max_q_cap)

    qs = prime_powers_with_index(r, max_q)
    notes: list[str] = []
    if not has_large_witnesses(r, qs):
        notes.append(
            f"WARNING: incomplete large witnesses up to max_q={max_q} "
            f"(cap={max_q_cap}); stable types may be provisional"
        )

    fps: dict[int, tuple] = {}
    for i, q in enumerate(qs):
        if progress:
            print(f"  [r={r}] fingerprint {i + 1}/{len(qs)} q={q}", flush=True)
        h = QuotientHyperfield.from_q_r(q, r)
        fps[q] = structure_fingerprint(h)

    by_res: dict[str, list[int]] = {}
    for q in qs:
        by_res.setdefault(baker_jin_residue(q, r), []).append(q)

    stable: dict[str, tuple] = {}
    for res, group in by_res.items():
        group_sorted = sorted(group)
        large = [q for q in group_sorted if q >= n12]
        if large:
            witness = large[-1]
            stable[res] = fps[witness]
            # Sanity: all large witnesses in the class share the fingerprint
            for q in large:
                if fps[q] != stable[res]:
                    notes.append(
                        f"WARNING: large q={q} in residue {res} disagrees with "
                        f"stable witness {witness} (contradicts Thm 1.1 / bug)"
                    )
        else:
            witness = group_sorted[-1]
            stable[res] = fps[witness]
            notes.append(
                f"residue {res}: no large q (>= {n12}); "
                f"provisional stable from q={witness}"
            )

    if r % 2 == 0 and "even_0" in stable and "even_1" in stable:
        if stable["even_0"] == stable["even_1"]:
            notes.append(
                "WARNING: even_0 and even_1 stable fingerprints coincide "
                "(unexpected for large q)"
            )

    exceptions: list[tuple[int, str]] = []
    for q in qs:
        res = baker_jin_residue(q, r)
        if res not in stable:
            continue
        # Only compare when we have a trustworthy stable type for this residue
        large_in_res = [x for x in by_res[res] if x >= n12]
        if not large_in_res and res in required_stable_residues(r):
            continue  # skip exception logic without a real stable type
        if fps[q] != stable[res]:
            exceptions.append((q, res))

    lower = paper_lower_bound_nr(r)
    large_all = sorted(q for q in qs if q >= n12)

    if exceptions:
        n_emp = max(q for q, _ in exceptions) + 1
        notes.append(
            "empirical N_r = 1 + max exception q "
            "(minimal threshold consistent with scanned counterexamples)"
        )
    elif large_all:
        # No exceptions: every small q already matches its large stable type.
        # True N_r is at most the first large prime power (upper bound).
        n_emp = large_all[0]
        notes.append(
            "no exceptions below the large regime; empirical value is an "
            f"upper bound (first large prime power q={n_emp})"
        )
    else:
        n_emp = None
        notes.append("no large prime powers in range")

    if lower is not None and n_emp is not None and n_emp < lower:
        notes.append(
            f"raising empirical N_r from {n_emp} to paper lower bound {lower}"
        )
        n_emp = lower

    # How sharp is Remark 1.2 on this range?
    if n_emp is not None:
        if n_emp >= n12:
            notes.append(
                f"Remark 1.2 bound {n12} is sharp or nearly sharp on this scan "
                f"(empirical {n_emp})"
            )
        else:
            notes.append(
                f"Remark 1.2 bound {n12} is not sharp on this scan "
                f"(empirical {n_emp} < {n12})"
            )

    return EmpiricalNrResult(
        r=r,
        max_q=max_q,
        n_r_remark_12=n12,
        n_r_safe=n_safe,
        n_r_lower_paper=lower,
        n_r_empirical=n_emp,
        prime_powers=qs,
        fingerprints=fps,
        stable_by_residue=stable,
        exceptions=exceptions,
        residue_classes={k: sorted(v) for k, v in by_res.items()},
        verified_up_to=max_q,
        notes=notes,
    )


def scan_empirical_nr(
    r_values: list[int] | range,
    max_q: int | None = None,
    *,
    max_q_cap: int = 8000,
    progress: bool = False,
) -> list[EmpiricalNrResult]:
    """Run empirical_nr for several r with resolved scan limits."""
    results = []
    for r in r_values:
        if progress:
            print(f"=== empirical N_r for r={r} ===", flush=True)
        results.append(
            empirical_nr(r, max_q, max_q_cap=max_q_cap, progress=progress)
        )
    return results


# ---------------------------------------------------------------------------
# Open question (4): Q_r among finite-field quotients
# ---------------------------------------------------------------------------


@dataclass
class QuotientClassCount:
    """Number of iso classes of finite-field quotients of order r+1."""

    r: int
    order: int  # r + 1
    max_q: int
    q_r_fin: int
    class_sizes: dict[str, int]  # fingerprint-id -> count of prime powers
    n_r_remark_12: int
    theoretically_complete: bool
    notes: list[str] = field(default_factory=list)

    def summary(self) -> str:
        lines = [
            f"r = {self.r} (hyperfield order {self.order})",
            f"  Q_r^{{fin}} (classes among F_q/G) = {self.q_r_fin}",
            f"  max_q scanned = {self.max_q}",
            f"  Remark 1.2 N_r = {self.n_r_remark_12}",
            f"  theoretically complete = {self.theoretically_complete}",
        ]
        for note in self.notes:
            lines.append(f"  note: {note}")
        return "\n".join(lines)


def count_finite_quotient_classes(
    r: int,
    max_q: int | None = None,
) -> QuotientClassCount:
    """
    Count distinct iso classes of F_q / G_r among prime powers q <= max_q.

    If max_q >= Remark 1.2 N_r, Theorem 1.1 implies every larger q falls into a
    class already represented by some q' with N_r <= q' <= max_q (for the
    stable residues). Sporadic classes only appear for q < N_r. Hence the count
    is the exact Q_r^{fin} whenever max_q >= N_r and every residue that occurs
    for large q is witnessed in range.
    """
    if r < 1:
        raise ValueError("r must be positive")
    n12 = baker_jin_bound_sharp(r) if r >= 2 else 0
    if max_q is None:
        max_q = max(n12, 32)

    qs = prime_powers_with_index(r, max_q)
    fps = []
    for q in qs:
        h = QuotientHyperfield.from_q_r(q, r)
        fps.append(structure_fingerprint(h))

    # Count unique fingerprints
    unique = set(fps)
    sizes: dict[str, int] = {}
    for i, fp in enumerate(unique):
        sizes[f"class_{i}"] = sum(1 for f in fps if f == fp)

    notes = []
    complete = r >= 2 and max_q >= n12
    if not complete:
        notes.append(
            "scan below Remark 1.2 bound; Q_r^{fin} may be incomplete"
        )
    else:
        notes.append(
            "max_q >= Remark 1.2 N_r: by Thm 1.1 this is the full Q_r^{fin} "
            "(finite-field quotients only; infinite-field quotients may add more)"
        )

    # Stable class count predicted by Baker–Jin
    if r >= 2 and r % 2 == 1:
        notes.append("Baker–Jin predicts 1 stable class (odd r) + sporadics")
    elif r >= 2:
        notes.append("Baker–Jin predicts 2 stable classes (even r) + sporadics")

    return QuotientClassCount(
        r=r,
        order=r + 1,
        max_q=max_q,
        q_r_fin=len(unique),
        class_sizes=sizes,
        n_r_remark_12=n12,
        theoretically_complete=complete,
        notes=notes,
    )


# ---------------------------------------------------------------------------
# CLI-style report
# ---------------------------------------------------------------------------


def report_open_question_3(
    r_max: int = 8,
    max_q_cap: int = 8000,
    *,
    progress: bool = True,
) -> str:
    """Human-readable report for open question (3)."""
    lines = [
        "Baker–Jin open question (3): true growth of N_r",
        "=" * 60,
        "Empirical N_r = 1 + max{q : fingerprint differs from large stable type},",
        "with scan limit extended past Remark 1.2 until large witnesses exist",
        f"(cap={max_q_cap}).",
        "",
    ]
    rows = scan_empirical_nr(
        range(2, r_max + 1),
        max_q_cap=max_q_cap,
        progress=progress,
    )
    for res in rows:
        lines.append(res.summary())
        lines.append("")

    lines.append("Comparison table")
    lines.append("-" * 60)
    lines.append(
        f"{'r':>3} {'lower':>8} {'empirical':>10} {'Remark1.2':>10} "
        f"{'r^4':>10} {'#exc':>5} {'max_q':>7}"
    )
    for res in rows:
        lines.append(
            f"{res.r:>3} {str(res.n_r_lower_paper):>8} "
            f"{str(res.n_r_empirical):>10} {res.n_r_remark_12:>10} "
            f"{res.n_r_safe:>10} {len(res.exceptions):>5} {res.max_q:>7}"
        )
    return "\n".join(lines)


def report_open_question_4(
    r_max: int = 8,
    max_q_cap: int = 8000,
    *,
    progress: bool = True,
) -> str:
    """Human-readable report for open question (4), finite-field part."""
    lines = [
        "Baker–Jin open question (4): Q_r among finite-field quotients",
        "=" * 60,
        "Q_r^{fin} = # iso classes of F_q/G of order r+1.",
        "Paper Q_r also counts infinite-field quotients; Q_r^{fin} <= Q_r.",
        "H_r (all hyperfields) requires full enumeration — not computed here.",
        "",
    ]
    lines.append(
        f"{'r':>3} {'order':>6} {'Q_r^fin':>8} {'N_r(1.2)':>10} "
        f"{'max_q':>7} {'complete':>9}"
    )
    for r in range(1, r_max + 1):
        if r >= 2:
            limit = resolve_scan_limit(r, max_q_cap=max_q_cap)
        else:
            limit = min(64, max_q_cap)
        if progress:
            print(f"=== Q_r^fin for r={r} (max_q={limit}) ===", flush=True)
        c = count_finite_quotient_classes(r, limit)
        lines.append(
            f"{c.r:>3} {c.order:>6} {c.q_r_fin:>8} {c.n_r_remark_12:>10} "
            f"{c.max_q:>7} {str(c.theoretically_complete):>9}"
        )
    return "\n".join(lines)


def main() -> None:
    import os
    from datetime import datetime, timezone

    r_max = int(os.environ.get("QH_R_MAX", "8"))
    max_q_cap = int(os.environ.get("QH_MAX_Q_CAP", "8000"))
    out_dir = os.environ.get("QH_OUT_DIR", "experiments_output")
    os.makedirs(out_dir, exist_ok=True)

    print(
        f"Running open-question probes: r_max={r_max}, max_q_cap={max_q_cap}",
        flush=True,
    )
    q3 = report_open_question_3(r_max=r_max, max_q_cap=max_q_cap, progress=True)
    print(q3, flush=True)
    q4 = report_open_question_4(r_max=r_max, max_q_cap=max_q_cap, progress=True)
    print(q4, flush=True)

    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    path = os.path.join(out_dir, f"open_questions_{stamp}.txt")
    latest = os.path.join(out_dir, "open_questions_latest.txt")
    text = (
        f"Generated (UTC): {stamp}\n"
        f"r_max={r_max} max_q_cap={max_q_cap}\n\n"
        f"{q3}\n\n{q4}\n"
    )
    with open(path, "w", encoding="utf-8") as f:
        f.write(text)
    with open(latest, "w", encoding="utf-8") as f:
        f.write(text)
    print(f"\nWrote {path}", flush=True)
    print(f"Wrote {latest}", flush=True)


if __name__ == "__main__":
    main()
