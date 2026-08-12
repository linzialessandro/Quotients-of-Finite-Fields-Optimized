"""
Extended computational probes beyond the arXiv:2608.03625 table defaults.

Purpose
-------
v0.2 made the paper algorithms fast enough to push **past** the published
range (r ≤ 8, n ≤ 7) and produce **new** finite-field data:

  1. Empirical N_r and sharpness ratios for higher index r
  2. Q_r^{fin} (iso classes of F_q/G_r) for higher r / orders
  3. Order atlases n ≥ 8 (no full H_n in the literature — pure Q^{fin})
  4. Stable (char, C-char) formula checks on the extended range

Paper tables (r ≤ 8, n ≤ 7) are regression targets; this module is for
exploration that the speedup was meant to enable.

Environment
-----------
  QH_R_MAX         max index r (default 12)
  QH_MAX_Q_CAP     scan cap (default 15000)
  QH_N_ATLAS_MAX   max hyperfield order for atlases (default r_max+1)
  QH_DO_ATLAS      0/1 (default 1)
  QH_OUT_DIR       output directory (default experiments_output)

CLI
---
  qh-extended
  python -m quotient_hyperfields.extended_experiments

Checkpoints
-----------
After each phase, writes under QH_OUT_DIR:

  extended_partial_<stamp>.txt   — cumulative text so far
  extended_nr_<stamp>.csv        — N_r table (phase 1)
  extended_qfin_<stamp>.csv      — Q_r^{fin} table (phase 2)
  extended_atlas_<stamp>.csv     — order atlases (phase 3)
  extended_<stamp>.txt           — final full report
  extended_latest.txt            — copy of final report
"""

from __future__ import annotations

import csv
import os
import time
from datetime import datetime, timezone
from pathlib import Path

from quotient_hyperfields.atlas import build_order_atlas
from quotient_hyperfields.experiments import (
    count_finite_quotient_classes,
    resolve_scan_limit,
    scan_empirical_nr,
)
from quotient_hyperfields.isomorphism import baker_jin_bound_sharp
from quotient_hyperfields.literature import AMERI_H_BY_ORDER, MASSOUROS_H7
from quotient_hyperfields.papers_experiments import experiment1_stable_formula


def _write_text(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8")


def _write_csv(path: Path, rows: list[dict]) -> None:
    if not rows:
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)


def format_nr_sharpness(rows) -> str:
    """Compact N_r table with sharpness ratios N_emp / N_{1.2}."""
    lines = [
        "N_r sharpness summary (extended range)",
        "=" * 60,
        "NEW beyond arXiv table (r≤8): any r≥9 rows below.",
        f"{'r':>3} {'N_emp':>8} {'N_1.2':>8} {'emp/1.2':>8} "
        f"{'r^4':>8} {'#exc':>5} {'max_q':>7}",
    ]
    for res in rows:
        emp = res.n_r_empirical
        ratio = (
            f"{emp / res.n_r_remark_12:.3f}"
            if emp is not None and res.n_r_remark_12
            else "—"
        )
        tag = " *" if res.r > 8 else ""
        lines.append(
            f"{res.r:>3} {str(emp):>8} {res.n_r_remark_12:>8} {ratio:>8} "
            f"{res.n_r_safe:>8} {len(res.exceptions):>5} {res.max_q:>7}{tag}"
        )
    lines.append("")
    lines.append("* = beyond published paper range (r≤8).")
    lines.append(
        "emp/1.2 < 1 ⇒ Remark 1.2 is not sharp on this scan "
        "(empirical threshold below the Weil bound)."
    )
    return "\n".join(lines)


def format_nr_detail(rows) -> str:
    lines = [
        "Baker–Jin open question (3): true growth of N_r (detail)",
        "=" * 60,
    ]
    for res in rows:
        lines.append(res.summary())
        lines.append("")
    return "\n".join(lines)


def nr_rows_to_csv(rows) -> list[dict]:
    out = []
    for res in rows:
        emp = res.n_r_empirical
        ratio = (
            round(emp / res.n_r_remark_12, 6)
            if emp is not None and res.n_r_remark_12
            else ""
        )
        out.append(
            {
                "r": res.r,
                "N_emp": emp if emp is not None else "",
                "N_remark12": res.n_r_remark_12,
                "N_safe_r4": res.n_r_safe,
                "N_lower": res.n_r_lower_paper if res.n_r_lower_paper else "",
                "emp_over_remark12": ratio,
                "n_exceptions": len(res.exceptions),
                "max_q": res.max_q,
                "beyond_paper": res.r > 8,
            }
        )
    return out


def report_qfin_table(
    r_max: int,
    max_q_cap: int,
    *,
    progress: bool = True,
) -> tuple[str, list[dict]]:
    """Q_r^{fin} table + CSV rows (new content for r>8 / order>7)."""
    lines = [
        "Baker–Jin open question (4): Q_r^{fin} (finite-field quotients)",
        "=" * 60,
        "Q_r^{fin} = # iso classes of F_q/G of order r+1.",
        "NEW: rows with r>8 (order>9) are beyond the arXiv census.",
        "",
        f"{'r':>3} {'order':>6} {'Q_r^fin':>8} {'N_r(1.2)':>10} "
        f"{'max_q':>7} {'complete':>9}",
    ]
    csv_rows: list[dict] = []
    for r in range(1, r_max + 1):
        if r >= 2:
            limit = resolve_scan_limit(r, max_q_cap=max_q_cap)
        else:
            limit = min(64, max_q_cap)
        if progress:
            print(f"=== Q_r^fin for r={r} (max_q={limit}) ===", flush=True)
        c = count_finite_quotient_classes(r, limit)
        tag = " *" if r > 8 else ""
        lines.append(
            f"{c.r:>3} {c.order:>6} {c.q_r_fin:>8} {c.n_r_remark_12:>10} "
            f"{c.max_q:>7} {str(c.theoretically_complete):>9}{tag}"
        )
        csv_rows.append(
            {
                "r": c.r,
                "order": c.order,
                "Q_fin": c.q_r_fin,
                "N_remark12": c.n_r_remark_12,
                "max_q": c.max_q,
                "complete": c.theoretically_complete,
                "beyond_paper": r > 8,
            }
        )
    lines.append("")
    lines.append("* = beyond published paper range (r≤8).")
    return "\n".join(lines), csv_rows


def report_extended_atlases(
    n_max: int = 13,
    max_q_cap: int = 15000,
    *,
    progress: bool = True,
) -> tuple[str, list[dict]]:
    """Finite-field Q^{fin} atlases for orders 2..n_max."""
    lines = [
        "Extended finite-field quotient atlases",
        "=" * 60,
        f"Orders n=2..{n_max}, max_q_cap={max_q_cap}.",
        "For n≥8 there is no published H_n in our literature table —",
        "these Q_fin values are new computational census data.",
        "",
        f"{'n':>3} {'r':>3} {'Q_fin':>6} {'H_lit':>7} {'ratio':>8} "
        f"{'complete':>9} {'N_r(1.2)':>10} {'max_q':>7}",
    ]
    csv_rows: list[dict] = []
    for n in range(2, n_max + 1):
        r = n - 1
        if progress:
            print(f"=== atlas n={n} (r={r}) ===", flush=True)
        t0 = time.perf_counter()
        atlas = build_order_atlas(n, max_q_cap=max_q_cap, progress=False)
        dt = time.perf_counter() - t0
        if n in AMERI_H_BY_ORDER:
            H: int | str = AMERI_H_BY_ORDER[n]
        elif n == 7:
            H = MASSOUROS_H7
        else:
            H = "—"
        if isinstance(H, int):
            ratio = f"{atlas.q_fin / H:.4f}"
            H_s = str(H)
            ratio_num: float | str = round(atlas.q_fin / H, 6)
        else:
            ratio = "—"
            H_s = "—"
            ratio_num = ""
        n12 = baker_jin_bound_sharp(r) if r >= 2 else 0
        tag = " *" if n > 7 else ""
        lines.append(
            f"{n:>3} {r:>3} {atlas.q_fin:>6} {H_s:>7} {ratio:>8} "
            f"{str(atlas.theoretically_complete):>9} {n12:>10} "
            f"{atlas.max_q:>7}  ({dt:.1f}s){tag}"
        )
        if progress:
            print(
                f"  Q_fin={atlas.q_fin} complete={atlas.theoretically_complete} "
                f"in {dt:.1f}s",
                flush=True,
            )
        csv_rows.append(
            {
                "order_n": n,
                "r": r,
                "Q_fin": atlas.q_fin,
                "H_literature": H if isinstance(H, int) else "",
                "ratio": ratio_num,
                "complete": atlas.theoretically_complete,
                "N_remark12": n12,
                "max_q": atlas.max_q,
                "beyond_paper": n > 7,
            }
        )
    lines.append("")
    lines.append("* = order beyond Massouros/Ameri comparison in the paper (n≤7).")
    return "\n".join(lines), csv_rows


def main() -> None:
    r_max = int(os.environ.get("QH_R_MAX", "12"))
    max_q_cap = int(os.environ.get("QH_MAX_Q_CAP", "15000"))
    n_atlas = int(os.environ.get("QH_N_ATLAS_MAX", str(r_max + 1)))
    do_atlas = os.environ.get("QH_DO_ATLAS", "1") not in ("0", "false", "False")
    out_dir = Path(os.environ.get("QH_OUT_DIR", "experiments_output"))
    out_dir.mkdir(parents=True, exist_ok=True)

    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    partial_path = out_dir / f"extended_partial_{stamp}.txt"
    print(
        f"Extended experiments: r_max={r_max}, max_q_cap={max_q_cap}, "
        f"n_atlas≤{n_atlas}, atlas={do_atlas}",
        flush=True,
    )
    print(f"Checkpoints → {out_dir} (stamp {stamp})", flush=True)
    t_all = time.perf_counter()

    parts: list[str] = []
    header = (
        f"Generated (UTC): {stamp}\n"
        f"software: quotient-hyperfields extended suite (v0.2+)\n"
        f"r_max={r_max} max_q_cap={max_q_cap} n_atlas_max={n_atlas} "
        f"atlas={do_atlas}\n"
        "Goal: new data beyond arXiv:2608.03625 tables (r≤8, n≤7).\n"
    )

    def checkpoint() -> None:
        body = header + f"wall_so_far={time.perf_counter() - t_all:.1f}s\n\n"
        body += "\n".join(parts) + "\n"
        _write_text(partial_path, body)
        print(f"  [checkpoint] {partial_path}", flush=True)

    # ---- Phase 1: N_r (single pass) ----
    print("\n[1/4] Empirical N_r (single pass) …", flush=True)
    t0 = time.perf_counter()
    nr_rows = scan_empirical_nr(
        range(2, r_max + 1), max_q_cap=max_q_cap, progress=True
    )
    parts.append(format_nr_sharpness(nr_rows))
    parts.append("")
    parts.append(format_nr_detail(nr_rows))
    _write_csv(out_dir / f"extended_nr_{stamp}.csv", nr_rows_to_csv(nr_rows))
    print(f"  done in {time.perf_counter() - t0:.1f}s", flush=True)
    checkpoint()

    # ---- Phase 2: Q_r^{fin} ----
    print("\n[2/4] Q_r^{fin} …", flush=True)
    t0 = time.perf_counter()
    qfin_text, qfin_csv = report_qfin_table(
        r_max=r_max, max_q_cap=max_q_cap, progress=True
    )
    parts.append("")
    parts.append(qfin_text)
    _write_csv(out_dir / f"extended_qfin_{stamp}.csv", qfin_csv)
    print(f"  done in {time.perf_counter() - t0:.1f}s", flush=True)
    checkpoint()

    # ---- Phase 3: atlases ----
    if do_atlas:
        print("\n[3/4] Extended atlases …", flush=True)
        t0 = time.perf_counter()
        atlas_text, atlas_csv = report_extended_atlases(
            n_max=n_atlas, max_q_cap=max_q_cap, progress=True
        )
        parts.append("")
        parts.append(atlas_text)
        _write_csv(out_dir / f"extended_atlas_{stamp}.csv", atlas_csv)
        print(f"  done in {time.perf_counter() - t0:.1f}s", flush=True)
        checkpoint()
    else:
        parts.append("")
        parts.append("(atlases skipped: QH_DO_ATLAS=0)")
        checkpoint()

    # ---- Phase 4: stable formula ----
    print("\n[4/4] Stable formula check …", flush=True)
    t0 = time.perf_counter()
    parts.append("")
    parts.append(
        experiment1_stable_formula(
            r_max=r_max, max_q_cap=max_q_cap, progress=True
        )
    )
    print(f"  done in {time.perf_counter() - t0:.1f}s", flush=True)

    wall = time.perf_counter() - t_all
    body = (
        header
        + f"wall_seconds={wall:.1f}\n"
        + "status=complete\n\n"
        + "\n".join(parts)
        + "\n"
    )
    final_path = out_dir / f"extended_{stamp}.txt"
    latest = out_dir / "extended_latest.txt"
    _write_text(final_path, body)
    _write_text(latest, body)
    # Keep partial as completed snapshot too
    _write_text(partial_path, body)

    print(body, flush=True)
    print(f"\nWrote {final_path}", flush=True)
    print(f"Wrote {latest}", flush=True)
    print(f"CSVs: extended_nr/qfin/atlas_{stamp}.csv", flush=True)
    print(f"Total wall time: {wall:.1f}s", flush=True)


if __name__ == "__main__":
    main()
