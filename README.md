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
- Char / C-char experiments (stable classes, gaps, parity theorem)
- Paper-driven experiments (stable formula, order-7 atlas, Ameri comparison,
  Massouros/Linzi criteria)

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

### Char / C-char invariant probes

```bash
qh-invariants
# or
python -c "from quotient_hyperfields.invariants_experiments import main; main()"
```

Runs experiments **A** (stable char/C-char by Baker–Jin class), **B** (gap census),
and **D** (even order ⇒ char 2). Same `experiments_output/` directory.

### Paper-driven experiments (tracks 1–5)

```bash
qh-papers
# or
python -c "from quotient_hyperfields.papers_experiments import main; main()"
```

| Track | Content | Primary paper |
|------:|---------|---------------|
| 1 | Stable (char, C-char) formula for large \(\mathbb{F}_q/G_r\) | Baker–Jin + Linzi |
| 2 | Order-7 finite-field quotient atlas vs \(H_7=277\) | Massouros |
| 3 | Ameri \(H_n\) vs \(Q^{\mathrm{fin}}\) for \(n\le 6\) | Ameri |
| 4 | Massouros Props 1–2 sum-cardinality checks | Massouros |
| 5 | Linzi Prop.7 / Cor.1 / Prop.9 char bounds | Linzi et al. |

Module docs: `quotient_hyperfields/papers_experiments.py` (header map).

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
  invariants_experiments.py  # char / C-char experiments A/B/D
  papers_experiments.py # paper tracks 1–5 (formula, atlas, criteria)
  criteria.py           # Massouros / Linzi checks
  atlas.py              # finite-field quotient class atlas
  literature.py         # Ameri / Massouros reference counts
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
