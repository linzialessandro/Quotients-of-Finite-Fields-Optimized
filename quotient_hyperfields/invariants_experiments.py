"""
Experiments on characteristic and C-characteristic of finite-field quotients.

A. Stable (char, C-char) by Baker–Jin class for each index r
B. Gap census: char - C-char over scanned prime powers
D. Parity check: even hyperfield order ⇒ char = 2
   (Kędzierski–Linzi–Stojałowska, finite hyperfields of even cardinality)

Uses the same scan limits as Baker–Jin open-question probes where helpful.
"""

from __future__ import annotations

import math
import os
from collections import Counter
from dataclasses import dataclass, field
from datetime import datetime, timezone

from quotient_hyperfields.experiments import (
    baker_jin_residue,
    prime_powers_with_index,
    required_stable_residues,
    resolve_scan_limit,
)
from quotient_hyperfields.hyperfield import QuotientHyperfield
from quotient_hyperfields.isomorphism import baker_jin_bound_sharp
from quotient_hyperfields.primes import factor_prime_power


def _fmt_inv(x: float) -> str:
    if isinstance(x, float) and math.isinf(x):
        return "∞"
    return str(int(x)) if float(x).is_integer() else str(x)


# ---------------------------------------------------------------------------
# A. Stable invariants by Baker–Jin class
# ---------------------------------------------------------------------------


@dataclass
class StableInvariantRow:
    r: int
    residue: str
    witness_q: int
    char: float
    c_char: float
    field_prime: int
    consistent: bool  # all large witnesses agree
    sample_qs: list[int] = field(default_factory=list)
    notes: list[str] = field(default_factory=list)


def stable_invariants_for_r(
    r: int,
    max_q: int | None = None,
    *,
    max_q_cap: int = 4000,
    progress: bool = False,
) -> list[StableInvariantRow]:
    """
    For each Baker–Jin stable residue of index r, compute (char, C-char)
    on large witnesses (q >= Remark 1.2 N_r) and check constancy.
    """
    if r < 2:
        raise ValueError("r >= 2 required")
    n12 = baker_jin_bound_sharp(r)
    limit = resolve_scan_limit(r, max_q, max_q_cap=max_q_cap)
    qs = prime_powers_with_index(r, limit)
    large_by_res: dict[str, list[int]] = {}
    for q in qs:
        if q < n12:
            continue
        res = baker_jin_residue(q, r)
        large_by_res.setdefault(res, []).append(q)

    rows: list[StableInvariantRow] = []
    for res in required_stable_residues(r):
        group = sorted(large_by_res.get(res, []))
        if not group:
            rows.append(
                StableInvariantRow(
                    r=r,
                    residue=res,
                    witness_q=-1,
                    char=math.inf,
                    c_char=math.inf,
                    field_prime=-1,
                    consistent=False,
                    notes=[f"no large q >= {n12} for residue {res} up to {limit}"],
                )
            )
            continue

        values: list[tuple[float, float]] = []
        for q in group:
            if progress:
                print(f"  [A r={r} {res}] q={q}", flush=True)
            h = QuotientHyperfield.from_q_r(q, r)
            values.append((h.char, h.c_char))

        witness = group[-1]
        h_w = QuotientHyperfield.from_q_r(witness, r)
        p, _ = factor_prime_power(witness)
        consistent = all(v == values[-1] for v in values)
        notes = []
        if not consistent:
            notes.append(f"INCONSISTENT large invariants: {set(values)}")
        if h_w.c_char > h_w.char:
            notes.append("WARNING: C-char > char (should be impossible)")

        rows.append(
            StableInvariantRow(
                r=r,
                residue=res,
                witness_q=witness,
                char=h_w.char,
                c_char=h_w.c_char,
                field_prime=p,
                consistent=consistent,
                sample_qs=group,
                notes=notes,
            )
        )
    return rows


def report_stable_invariants(
    r_max: int = 10,
    max_q_cap: int = 4000,
    *,
    progress: bool = True,
) -> str:
    lines = [
        "Experiment A: stable (char, C-char) by Baker–Jin class",
        "=" * 60,
        "Large q means q >= Remark 1.2 N_r. Invariants should be constant",
        "on each stable class (iso invariants).",
        "",
        f"{'r':>3} {'residue':>10} {'char':>6} {'C-char':>6} "
        f"{'p(wit)':>7} {'witness_q':>10} {'ok':>4} {'#large':>6}",
    ]
    for r in range(2, r_max + 1):
        if progress:
            print(f"=== A stable invariants r={r} ===", flush=True)
        for row in stable_invariants_for_r(r, max_q_cap=max_q_cap, progress=progress):
            lines.append(
                f"{row.r:>3} {row.residue:>10} {_fmt_inv(row.char):>6} "
                f"{_fmt_inv(row.c_char):>6} {row.field_prime:>7} "
                f"{row.witness_q:>10} {str(row.consistent):>4} "
                f"{len(row.sample_qs):>6}"
            )
            for n in row.notes:
                lines.append(f"       note: {n}")
    return "\n".join(lines)


# ---------------------------------------------------------------------------
# B. Gap census
# ---------------------------------------------------------------------------


@dataclass
class GapCensus:
    r: int
    max_q: int
    n_samples: int
    gap_counts: Counter  # gap -> count
    cchar_counts: Counter
    char_counts: Counter
    equal_count: int  # char == C-char
    cchar_one_count: int
    pairs: Counter  # (char, c_char) -> count

    def summary_lines(self) -> list[str]:
        lines = [
            f"r = {self.r}  (n={self.n_samples} prime powers, max_q={self.max_q})",
            f"  equal char=C-char: {self.equal_count}/{self.n_samples}",
            f"  C-char = 1:        {self.cchar_one_count}/{self.n_samples}",
            f"  gap histogram (char - C-char): "
            + ", ".join(
                f"{g}:{c}" for g, c in sorted(self.gap_counts.items())
            ),
            f"  top (char, C-char) pairs: "
            + ", ".join(
                f"({_fmt_inv(a)},{_fmt_inv(b)}):{n}"
                for (a, b), n in self.pairs.most_common(8)
            ),
        ]
        return lines


def gap_census_for_r(
    r: int,
    max_q: int | None = None,
    *,
    max_q_cap: int = 4000,
    progress: bool = False,
) -> GapCensus:
    """Scan all prime powers with r | (q-1) up to resolved limit; record gaps."""
    if r < 2:
        raise ValueError("r >= 2 required")
    limit = resolve_scan_limit(r, max_q, max_q_cap=max_q_cap)
    qs = prime_powers_with_index(r, limit)

    gap_counts: Counter = Counter()
    cchar_counts: Counter = Counter()
    char_counts: Counter = Counter()
    pairs: Counter = Counter()
    equal = 0
    cchar_one = 0

    for i, q in enumerate(qs):
        if progress:
            print(f"  [B r={r}] {i + 1}/{len(qs)} q={q}", flush=True)
        h = QuotientHyperfield.from_q_r(q, r)
        ch, cc = h.char, h.c_char
        if math.isinf(ch) or math.isinf(cc):
            gap = None
        else:
            gap = int(ch - cc)
            gap_counts[gap] += 1
        char_counts[ch] += 1
        cchar_counts[cc] += 1
        pairs[(ch, cc)] += 1
        if ch == cc:
            equal += 1
        if cc == 1.0:
            cchar_one += 1

    return GapCensus(
        r=r,
        max_q=limit,
        n_samples=len(qs),
        gap_counts=gap_counts,
        cchar_counts=cchar_counts,
        char_counts=char_counts,
        equal_count=equal,
        cchar_one_count=cchar_one,
        pairs=pairs,
    )


def report_gap_census(
    r_max: int = 8,
    max_q_cap: int = 4000,
    *,
    progress: bool = True,
) -> str:
    lines = [
        "Experiment B: gap census (char − C-char)",
        "=" * 60,
        "Over all prime powers q with r | (q-1) up to the usual scan limit.",
        "",
    ]
    for r in range(2, r_max + 1):
        if progress:
            print(f"=== B gap census r={r} ===", flush=True)
        census = gap_census_for_r(r, max_q_cap=max_q_cap, progress=progress)
        lines.extend(census.summary_lines())
        lines.append("")
    return "\n".join(lines)


# ---------------------------------------------------------------------------
# D. Parity theorem: even order ⇒ char = 2
# ---------------------------------------------------------------------------


@dataclass
class ParityCheckResult:
    r: int  # odd r ⇒ order r+1 even
    order: int
    max_q: int
    n_checked: int
    n_ok: int
    failures: list[tuple[int, float]]  # (q, char)

    @property
    def passed(self) -> bool:
        return self.n_checked > 0 and len(self.failures) == 0


def parity_char_check_for_r(
    r: int,
    max_q: int | None = None,
    *,
    max_q_cap: int = 4000,
    progress: bool = False,
) -> ParityCheckResult:
    """
    For odd r (hyperfield order even), verify char == 2 on all scanned quotients.
    """
    if r % 2 == 0:
        raise ValueError("parity check applies to odd r (even hyperfield order)")
    limit = resolve_scan_limit(r, max_q, max_q_cap=max_q_cap)
    qs = prime_powers_with_index(r, limit)
    failures: list[tuple[int, float]] = []
    for i, q in enumerate(qs):
        if progress:
            print(f"  [D r={r}] {i + 1}/{len(qs)} q={q}", flush=True)
        h = QuotientHyperfield.from_q_r(q, r)
        if h.char != 2.0:
            failures.append((q, h.char))
    return ParityCheckResult(
        r=r,
        order=r + 1,
        max_q=limit,
        n_checked=len(qs),
        n_ok=len(qs) - len(failures),
        failures=failures,
    )


def report_parity_check(
    r_max: int = 11,
    max_q_cap: int = 4000,
    *,
    progress: bool = True,
) -> str:
    lines = [
        "Experiment D: even order ⇒ char = 2",
        "=" * 60,
        "For odd index r, |K|=r+1 is even. Theorem (char paper): finite",
        "hyperfields of even cardinality have characteristic 2.",
        "",
        f"{'r':>3} {'order':>6} {'checked':>8} {'ok':>6} {'fail':>6} {'status':>8}",
    ]
    for r in range(3, r_max + 1, 2):  # odd r only
        if progress:
            print(f"=== D parity r={r} ===", flush=True)
        res = parity_char_check_for_r(r, max_q_cap=max_q_cap, progress=progress)
        status = "PASS" if res.passed else "FAIL"
        lines.append(
            f"{res.r:>3} {res.order:>6} {res.n_checked:>8} {res.n_ok:>6} "
            f"{len(res.failures):>6} {status:>8}"
        )
        for q, ch in res.failures[:10]:
            lines.append(f"       failure: q={q} char={_fmt_inv(ch)}")
    return "\n".join(lines)


# ---------------------------------------------------------------------------
# Combined report + CLI
# ---------------------------------------------------------------------------


def full_invariants_report(
    r_max_stable: int = 10,
    r_max_gap: int = 8,
    r_max_parity: int = 11,
    max_q_cap: int = 4000,
    *,
    progress: bool = True,
) -> str:
    parts = [
        report_stable_invariants(
            r_max=r_max_stable, max_q_cap=max_q_cap, progress=progress
        ),
        "",
        report_gap_census(r_max=r_max_gap, max_q_cap=max_q_cap, progress=progress),
        report_parity_check(
            r_max=r_max_parity, max_q_cap=max_q_cap, progress=progress
        ),
    ]
    return "\n".join(parts)


def main() -> None:
    r_stable = int(os.environ.get("QH_R_MAX_STABLE", "10"))
    r_gap = int(os.environ.get("QH_R_MAX_GAP", "8"))
    r_parity = int(os.environ.get("QH_R_MAX_PARITY", "11"))
    max_q_cap = int(os.environ.get("QH_MAX_Q_CAP", "4000"))
    out_dir = os.environ.get("QH_OUT_DIR", "experiments_output")
    os.makedirs(out_dir, exist_ok=True)

    print(
        f"Running invariant experiments: "
        f"stable r<= {r_stable}, gap r<= {r_gap}, parity r<= {r_parity}, "
        f"max_q_cap={max_q_cap}",
        flush=True,
    )
    text = full_invariants_report(
        r_max_stable=r_stable,
        r_max_gap=r_gap,
        r_max_parity=r_parity,
        max_q_cap=max_q_cap,
        progress=True,
    )
    print(text, flush=True)

    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    path = os.path.join(out_dir, f"invariants_{stamp}.txt")
    latest = os.path.join(out_dir, "invariants_latest.txt")
    header = (
        f"Generated (UTC): {stamp}\n"
        f"r_stable={r_stable} r_gap={r_gap} r_parity={r_parity} "
        f"max_q_cap={max_q_cap}\n\n"
    )
    body = header + text + "\n"
    with open(path, "w", encoding="utf-8") as f:
        f.write(body)
    with open(latest, "w", encoding="utf-8") as f:
        f.write(body)
    print(f"\nWrote {path}", flush=True)
    print(f"Wrote {latest}", flush=True)


if __name__ == "__main__":
    main()
