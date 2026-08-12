# Extended computational results (v0.2+)

Data produced with the accelerated library (prime-field backend, large-q shortcuts).
**Beyond arXiv:2608.03625 tables** (paper: empirical N_r for r≤8, Q_fin for n≤7).

## Files

| File | Content |
|------|---------|
| `combined_nr_r2_to_12.csv` | Empirical N_r vs Remark 1.2 for r=2..12 |
| `combined_qfin_r1_to_12.csv` | Q_r^{fin} for r=1..12 (order n=r+1) |
| `extended_nr_r10.csv` / `extended_r10_report.txt` | Full suite r≤10 |
| `extended_nr_r11_12.csv` / `extended_r11_12_report.txt` | Partial suite r=11..12 |
| `extended_atlas_*.csv` | Order atlases |

## Headline findings

- Remark 1.2 Weil bound is sharp for r=2,3 and increasingly loose for larger r
  (emp/N_1.2 ≈ 0.18 at r=9, ≈ 0.14 at r=12).
- Finite-field quotient class counts stay modest: Q_fin ∈ {9..36} for orders 8..13.
- Stable (char, C-char) formula holds on all large witnesses checked through r=12.

Protocol: `qh-extended` / `quotient_hyperfields.extended_experiments`.
