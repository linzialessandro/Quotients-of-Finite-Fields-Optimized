"""
Paper-driven experiments (tracks 1–5).

===========================================================================
Documentation map
===========================================================================

1. Stable (char, C-char) formula check  — Baker–Jin + Linzi
   ---------------------------------------------------------
   Conjecture (supported by prior invariant runs for r≤8):
     For large F_q / G_r (q ≥ N_r from Remark 1.2):
       - odd r:  (char, C-char) = (2, 1)
       - even r, q ≡ 1 (mod 2r):     (2, 1)
       - even r, q ≡ r+1 (mod 2r):   (3, 1)
   This experiment verifies the pattern on every large prime power in range.

2. Order-7 finite-field quotient atlas  — Massouros + Baker–Jin
   -----------------------------------------------------------
   Massouros: H_7 = 277 hyperfields of order 7.
   We compute Q_6^{fin} = # iso classes of F_q/G_6 and list each class
   (char, C-char, sample q). Ratio Q_6^{fin}/277 is a data point for
   Baker–Jin open question (4).

3. Ameri H_n vs Q^{fin} for orders 2–6  — Ameri Table 1
   ----------------------------------------------------
   Compare literature counts of *all* hyperfields to our finite-field
   quotient class counts (lower bound on field-quotients).

4. Massouros sum-cardinality criteria  — Massouros Props 1–2
   ---------------------------------------------------------
   On a sample of quotients, verify |x+y| ≤ |G| always, and Prop.2 when
   the difference hypothesis holds.

5. Linzi Prop.7 / Cor.1 / Prop.9        — char paper
   -------------------------------------------------
   char ≤ n for every n>1 dividing |G|; C-char same; for p odd,
   char=2 ⇔ |G| even.

Real quotients / positive cones are intentionally out of scope.

CLI
---
  qh-papers
  python -c "from quotient_hyperfields.papers_experiments import main; main()"

Environment
-----------
  QH_MAX_Q_CAP   scan cap (default 3500)
  QH_R_MAX_FORMULA  max r for experiment 1 (default 10)
  QH_OUT_DIR     output directory (default experiments_output)
"""

from __future__ import annotations

import math
import os
from dataclasses import dataclass, field
from datetime import datetime, timezone

from quotient_hyperfields.atlas import OrderAtlas, build_order_atlas
from quotient_hyperfields.criteria import (
    CheckResult,
    verify_all_paper_criteria,
    verify_linzi_criteria,
    verify_massouros_criteria,
)
from quotient_hyperfields.experiments import (
    baker_jin_residue,
    prime_powers_with_index,
    resolve_scan_limit,
)
from quotient_hyperfields.hyperfield import QuotientHyperfield
from quotient_hyperfields.isomorphism import baker_jin_bound_sharp
from quotient_hyperfields.literature import (
    AMERI_H_BY_ORDER,
    BAKER_JIN_ORDER_LE4_NOTES,
    MASSOUROS_H7,
    MASSOUROS_ORDER7_QUOTIENT_EXAMPLES,
)


# ---------------------------------------------------------------------------
# 1. Stable char / C-char formula
# ---------------------------------------------------------------------------


def expected_stable_invariants(r: int, q: int) -> tuple[float, float]:
    """
    Predicted (char, C-char) for large F_q/G_r under the experimental formula.

    See module docstring, experiment 1.
    """
    if r < 2:
        raise ValueError("r >= 2")
    if r % 2 == 1:
        return (2.0, 1.0)
    # even r: class by q mod 2r
    mod = q % (2 * r)
    if mod == 1:
        return (2.0, 1.0)
    if mod == r + 1:
        return (3.0, 1.0)
    # should not occur when r | (q-1)
    return (float("nan"), float("nan"))


@dataclass
class FormulaCheckResult:
    r: int
    max_q: int
    n_large: int
    n_ok: int
    failures: list[str] = field(default_factory=list)

    @property
    def passed(self) -> bool:
        return self.n_large > 0 and len(self.failures) == 0


def experiment1_stable_formula(
    r_max: int = 10,
    max_q_cap: int = 3500,
    *,
    progress: bool = True,
) -> str:
    """
    Experiment 1 — verify stable (char, C-char) formula for large quotients.
    """
    lines = [
        "EXPERIMENT 1 — Stable (char, C-char) formula",
        "=" * 60,
        "Conjecture: for q ≥ N_r (Remark 1.2),",
        "  odd r              → (char, C-char) = (2, 1)",
        "  even r, q≡1(2r)    → (2, 1)",
        "  even r, q≡r+1(2r)  → (3, 1)",
        "References: Baker–Jin Thm 1.1 (class structure); Linzi et al. (invariants).",
        "",
        f"{'r':>3} {'#large':>7} {'ok':>5} {'fail':>5} {'status':>8}",
    ]
    all_pass = True
    for r in range(2, r_max + 1):
        if progress:
            print(f"=== Exp1 formula r={r} ===", flush=True)
        n12 = baker_jin_bound_sharp(r)
        limit = resolve_scan_limit(r, max_q_cap=max_q_cap)
        qs = [q for q in prime_powers_with_index(r, limit) if q >= n12]
        failures: list[str] = []
        for q in qs:
            if progress:
                print(f"  [1] r={r} q={q}", flush=True)
            h = QuotientHyperfield.from_q_r(q, r)
            exp_ch, exp_cc = expected_stable_invariants(r, q)
            if (h.char, h.c_char) != (exp_ch, exp_cc):
                failures.append(
                    f"q={q} res={baker_jin_residue(q, r)}: "
                    f"got ({h.char},{h.c_char}) expected ({exp_ch},{exp_cc})"
                )
        status = "PASS" if qs and not failures else ("FAIL" if failures else "EMPTY")
        if status != "PASS":
            all_pass = False
        lines.append(
            f"{r:>3} {len(qs):>7} {len(qs) - len(failures):>5} "
            f"{len(failures):>5} {status:>8}"
        )
        for f in failures[:5]:
            lines.append(f"       {f}")
    lines.append("")
    lines.append(
        "Overall: "
        + ("formula holds on all tested large quotients." if all_pass else "FAILURES present.")
    )
    return "\n".join(lines)


# ---------------------------------------------------------------------------
# 2. Order-7 atlas vs Massouros 277
# ---------------------------------------------------------------------------


def experiment2_order7_atlas(
    max_q_cap: int = 3500,
    *,
    progress: bool = True,
) -> str:
    """
    Experiment 2 — finite-field quotients of order 7 vs Massouros H_7=277.
    """
    lines = [
        "EXPERIMENT 2 — Order-7 finite-field quotient atlas (Massouros)",
        "=" * 60,
        f"Massouros–Massouros: H_7 = {MASSOUROS_H7} hyperfields of order 7.",
        "We list iso classes of F_q / G_6 (index r=6 ⇒ order 7).",
        "References: arXiv:2412.11331; Baker–Jin open question (4).",
        "",
    ]
    if progress:
        print("=== Exp2 order-7 atlas ===", flush=True)
    atlas = build_order_atlas(7, max_q_cap=max_q_cap, progress=progress)
    lines.append(atlas.summary())
    lines.append("")
    ratio = atlas.q_fin / MASSOUROS_H7
    lines.append(
        f"Q_6^{{fin}} / H_7 = {atlas.q_fin}/{MASSOUROS_H7} ≈ {ratio:.4f} "
        f"({100 * ratio:.2f}%)"
    )
    lines.append(
        "Interpretation: only finite-field quotients; infinite-field quotients "
        "could add more to the paper's Q_6. Still Q_6^{fin}/H_7 is a concrete "
        "lower bound on 'quotient share' at order 7."
    )
    lines.append("")
    lines.append("Massouros narrative examples (q with index-6 subgroup):")
    # Check which example q fall into our atlas
    fp_by_q = {}
    for c in atlas.classes:
        for q in c.sample_qs:
            fp_by_q[q] = c.class_id
    for q, note in MASSOUROS_ORDER7_QUOTIENT_EXAMPLES:
        if q in fp_by_q:
            lines.append(f"  q={q}: in our class {fp_by_q[q]}  ({note})")
        else:
            # may be outside scan limit
            lines.append(f"  q={q}: not in scan (max_q={atlas.max_q})  ({note})")
    n_stable = sum(1 for c in atlas.classes if c.is_stable_large)
    n_sp = atlas.q_fin - n_stable
    lines.append("")
    lines.append(
        f"Stable large classes: {n_stable}  |  sporadic classes: {n_sp}  "
        f"(Baker–Jin predicts 2 stable for even r=6)"
    )
    return "\n".join(lines)


# ---------------------------------------------------------------------------
# 3. Ameri H_n vs Q^{fin}
# ---------------------------------------------------------------------------


def experiment3_ameri_comparison(
    max_q_cap: int = 3500,
    *,
    progress: bool = True,
) -> str:
    """
    Experiment 3 — Ameri enumeration vs finite-field quotient counts.
    """
    lines = [
        "EXPERIMENT 3 — Ameri H_n vs Q^{fin} (finite-field quotients)",
        "=" * 60,
        "Ameri et al. Table 1: all hyperfields of order n ≤ 6.",
        "Q^{fin}: iso classes of F_q/G with |K|=n (this library).",
        "Baker–Jin Thm 1.3 notes for order ≤ 4 (any field, not only finite).",
        "",
        f"{'n':>3} {'r':>3} {'H_n(Ameri)':>12} {'Q^fin':>7} {'Q/H':>8} "
        f"{'complete':>9}",
    ]
    for n in sorted(AMERI_H_BY_ORDER):
        if progress:
            print(f"=== Exp3 Ameri n={n} ===", flush=True)
        h_n = AMERI_H_BY_ORDER[n]
        atlas = build_order_atlas(n, max_q_cap=max_q_cap, progress=progress)
        ratio = atlas.q_fin / h_n if h_n else float("nan")
        lines.append(
            f"{n:>3} {atlas.r:>3} {h_n:>12} {atlas.q_fin:>7} {ratio:>8.3f} "
            f"{str(atlas.theoretically_complete):>9}"
        )
    lines.append("")
    lines.append("Baker–Jin Thm 1.3 (order ≤ 4, narrative):")
    for n, note in sorted(BAKER_JIN_ORDER_LE4_NOTES.items()):
        lines.append(f"  n={n}: {note}")
    lines.append("")
    lines.append(
        "Note: Q^fin counts only finite-field quotients. Baker–Jin/Ameri "
        "also allow infinite-field quotients (e.g. signs at order 3), so "
        "Q_full ≥ Q^fin and Q^fin/H_n underestimates the full quotient share."
    )
    return "\n".join(lines)


# ---------------------------------------------------------------------------
# 4. Massouros criteria battery
# ---------------------------------------------------------------------------


def experiment4_massouros_criteria(
    sample: list[tuple[int, int]] | None = None,
    *,
    progress: bool = True,
) -> str:
    """
    Experiment 4 — Massouros Props 1–2 on sample quotients.

    sample: list of (q, r) pairs. Default: a mix of small and large.
    """
    if sample is None:
        sample = [
            (17, 2),
            (19, 2),
            (13, 3),
            (31, 3),
            (53, 4),
            (61, 4),
            (11, 5),
            (181, 5),
            (97, 6),
            (157, 6),
            (7, 1),
            (5, 1),  # field as hyperfield, r=4 actually for F_5 d=1 → r=4
        ]
        # fix: F_5 with d=1 means r=4
        sample = [
            (17, 2),
            (19, 2),
            (13, 3),
            (31, 3),
            (53, 4),
            (101, 4),
            (181, 5),
            (97, 6),
            (157, 6),
            (29, 1),  # Krasner-type order 2
            (5, 4),   # F_5 as field (d=1, r=4)
            (7, 6),   # F_7 order 7
        ]

    lines = [
        "EXPERIMENT 4 — Massouros Propositions 1–2 (sum cardinalities)",
        "=" * 60,
        "Prop.1: |x+y| ≤ |G| for all x,y in a quotient F/G.",
        "Prop.2: if difference hypothesis holds, |x+y|=|G| for non-opposite unequal pairs.",
        "Reference: Massouros–Massouros arXiv:2412.11331.",
        "",
        f"{'q':>6} {'r':>4} {'d':>6} {'Prop1':>7} {'Prop2':>7}",
    ]
    n_ok = 0
    for q, r in sample:
        if progress:
            print(f"=== Exp4 Massouros q={q} r={r} ===", flush=True)
        try:
            h = QuotientHyperfield.from_q_r(q, r)
        except Exception as e:
            lines.append(f"{q:>6} {r:>4} {'—':>6} ERROR:{e}")
            continue
        results = verify_massouros_criteria(h)
        p1, p2 = results[0], results[1]
        if p1.passed and p2.passed:
            n_ok += 1
        lines.append(
            f"{q:>6} {r:>4} {h.d:>6} "
            f"{'PASS' if p1.passed else 'FAIL':>7} "
            f"{'PASS' if p2.passed else 'FAIL':>7}"
        )
        for res in results:
            if not res.passed:
                for f in res.failures[:3]:
                    lines.append(f"         {res.name}: {f}")
    lines.append("")
    lines.append(f"Samples fully passing both checks: {n_ok}/{len(sample)}")
    return "\n".join(lines)


# ---------------------------------------------------------------------------
# 5. Linzi Prop 7 / Cor 1 / Prop 9
# ---------------------------------------------------------------------------


def experiment5_linzi_bounds(
    sample: list[tuple[int, int]] | None = None,
    *,
    progress: bool = True,
) -> str:
    """
    Experiment 5 — Linzi Prop.7, Cor.1, Prop.9 on sample quotients.
    """
    if sample is None:
        sample = [
            (17, 2),
            (19, 2),
            (13, 3),
            (31, 3),
            (5, 4),  # field F_5
            (25, 4),
            (11, 5),
            (181, 5),
            (97, 6),
            (13, 6),
            (29, 1),
            (8, 7),  # 2^3
            (9, 2),
            (49, 6),
        ]

    lines = [
        "EXPERIMENT 5 — Linzi et al. Prop.7 / Cor.1 / Prop.9",
        "=" * 60,
        "Prop.7: n|d, n>1 ⇒ char ≤ n.",
        "Cor.1:  p odd ⇒ (char=2 ⇔ d even).",
        "Prop.9: n|d, n>1 ⇒ C-char ≤ n.",
        "Reference: Kędzierski–Linzi–Stojałowska, Mathematics 2023.",
        "",
        f"{'q':>6} {'r':>4} {'d':>6} {'char':>5} {'Cc':>5} "
        f"{'P7':>5} {'C1':>5} {'P9':>5}",
    ]
    n_all = 0
    for q, r in sample:
        if progress:
            print(f"=== Exp5 Linzi q={q} r={r} ===", flush=True)
        try:
            h = QuotientHyperfield.from_q_r(q, r)
        except Exception as e:
            lines.append(f"{q:>6} {r:>4} ERROR {e}")
            continue
        results = verify_linzi_criteria(h)
        p7, c1, p9 = results
        if p7.passed and c1.passed and p9.passed:
            n_all += 1
        ch = "∞" if math.isinf(h.char) else str(int(h.char))
        cc = "∞" if math.isinf(h.c_char) else str(int(h.c_char))
        lines.append(
            f"{q:>6} {r:>4} {h.d:>6} {ch:>5} {cc:>5} "
            f"{'Y' if p7.passed else 'N':>5} "
            f"{'Y' if c1.passed else 'N':>5} "
            f"{'Y' if p9.passed else 'N':>5}"
        )
        for res in results:
            if not res.passed:
                lines.append(f"         FAIL {res.name}: {res.detail}")
                for f in res.failures[:3]:
                    lines.append(f"           {f}")
    lines.append("")
    lines.append(f"Samples passing all three checks: {n_all}/{len(sample)}")
    return "\n".join(lines)


# ---------------------------------------------------------------------------
# Combined runner
# ---------------------------------------------------------------------------


def full_papers_report(
    max_q_cap: int = 3500,
    r_max_formula: int = 10,
    *,
    progress: bool = True,
) -> str:
    parts = [
        experiment1_stable_formula(
            r_max=r_max_formula, max_q_cap=max_q_cap, progress=progress
        ),
        "",
        experiment2_order7_atlas(max_q_cap=max_q_cap, progress=progress),
        "",
        experiment3_ameri_comparison(max_q_cap=max_q_cap, progress=progress),
        "",
        experiment4_massouros_criteria(progress=progress),
        "",
        experiment5_linzi_bounds(progress=progress),
    ]
    return "\n".join(parts)


def main() -> None:
    max_q_cap = int(os.environ.get("QH_MAX_Q_CAP", "3500"))
    r_formula = int(os.environ.get("QH_R_MAX_FORMULA", "8"))
    out_dir = os.environ.get("QH_OUT_DIR", "experiments_output")
    os.makedirs(out_dir, exist_ok=True)

    print(
        f"Paper experiments 1–5: r_formula≤{r_formula}, max_q_cap={max_q_cap}",
        flush=True,
    )
    text = full_papers_report(
        max_q_cap=max_q_cap,
        r_max_formula=r_formula,
        progress=True,
    )
    print(text, flush=True)

    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    path = os.path.join(out_dir, f"papers_experiments_{stamp}.txt")
    latest = os.path.join(out_dir, "papers_experiments_latest.txt")
    header = (
        f"Generated (UTC): {stamp}\n"
        f"r_max_formula={r_formula} max_q_cap={max_q_cap}\n"
        f"Tracks: 1 stable formula | 2 order-7 atlas | 3 Ameri vs Qfin | "
        f"4 Massouros criteria | 5 Linzi bounds\n\n"
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
