# Quotient Hyperfields of Finite Fields

Research-grade Python library for **Krasner quotient hyperfields** of the form

\[
K = \mathbb{F}_q / G_d,
\]

where \(G_d \le \mathbb{F}_q^\times\) is the unique subgroup of order \(d\), the index is \(r = (q-1)/d\), and \(|K| = r+1\).

## Features

- Construct \(K\) from \((p,k,d)\), \((q,d)\), or \((q,r)\)
- Hyperaddition and full additive tables
- **Layered isomorphism**
  - **Baker–Jin** (arXiv:1912.11496, Thm 1.1) for large \(q\), with **Remark 1.2** bound \(N_r\) (sharp Weil threshold; not the coarser \(r^4\))
  - **General** Aut\((C_r)\) check of \(1 \boxplus x\) tables (gold standard / sporadics)
  - **Auto** policy combining both
- **Characteristic** and **C-characteristic** (Kędzierski–Linzi–Stojałowska)
- Modest census: list triples and classify iso classes for fixed order \(n\)
- Thin interactive CLI
- Computational probes for Baker–Jin open questions (§5)

## Install

```bash
python3 -m venv .venv
source .venv/bin/activate
pip install -e ".[dev]"
```

Requires Python 3.10+.

## Quick start

```python
from quotient_hyperfields import QuotientHyperfield, are_isomorphic

# F_17 with index r=2 (order-3 hyperfield)
h1 = QuotientHyperfield.from_q_r(17, 2)
h2 = QuotientHyperfield.from_q_r(29, 2)

print(h1, h1.char, h1.c_char)
print(are_isomorphic(h1, h2))  # Baker–Jin: same even class

print(h1.addition_table())
```

## CLI

```bash
quotient-hyperfields
# or
python -m quotient_hyperfields.cli
# or
python main.py
```

### Baker–Jin open-question probes

Heavy experiment helpers live in `quotient_hyperfields.experiments` and are **not**
re-exported from the package root (avoids double-import warnings).

```bash
qh-open-questions
# or
python -c "from quotient_hyperfields.experiments import main; main()"
```

Reports are written under `experiments_output/` (gitignored; regenerate anytime).

## Mathematical references (in `papers/`)

| Paper | Role |
|-------|------|
| Baker–Jin, arXiv:1912.11496 | Iso criterion for large finite-field quotients (Remark 1.2 \(N_r\); safe \(r^4\) available) |
| Kędzierski–Linzi–Stojałowska | Char / C-char definitions and motivation |
| Ameri–Eyvazi–Hoskova | Enumeration order ≤ 6; iso via \(1+x\) |
| Massouros, arXiv:2412.11331 | Full hyperfields of order 7; quotient examples |

## Project layout

```
quotient_hyperfields/   # installable package
  hyperfield.py         # QuotientHyperfield
  isomorphism.py        # Baker–Jin / general / auto
  census.py             # triples + classification
  experiments.py        # open-question probes (Q3 / Q_r^fin)
  cli.py                # thin interactive UI
  visualization.py      # census plots
tests/                  # paper-backed tests
papers/                 # reference PDFs
main.py                 # CLI entry (compat)
```

## Tests

```bash
pytest -q
```

## Scope (v1)

**In scope:** finite-field quotient hyperfields only; laptop-scale experiments; certified iso; char/C-char; empirical \(N_r\) / \(Q_r^{\mathrm{fin}}\) probes.

**Out of scope:** full enumeration of all finite hyperfields; GUIs; non-Python cores.
