# Quotient Hyperfields of Finite Fields

[![arXiv](https://img.shields.io/badge/arXiv-2608.03625-b31b1b.svg)](https://arxiv.org/abs/2608.03625)
[![DOI](https://img.shields.io/badge/DOI-10.48550/arXiv.2608.03625-blue.svg)](https://doi.org/10.48550/arXiv.2608.03625)
[![Version](https://img.shields.io/badge/version-0.2.0-green.svg)](https://github.com/linzialessandro/Quotients-of-Finite-Fields-Optimized/releases/tag/v0.2.0)

Research-grade Python library for **Krasner quotient hyperfields** of the form

$$
K = \mathbb{F}_q / G_d,
$$

where $G_d \le \mathbb{F}_q^\times$ is the unique subgroup of order $d$, the index is $r = (q-1)/d$, and $|K| = r+1$.

## Paper

Companion software for:

> Alessandro Linzi, *Finite-field Krasner quotients: isomorphism thresholds, characteristics, and censuses*, arXiv:[2608.03625](https://arxiv.org/abs/2608.03625) [math.RA], 2026.

- Abstract: https://arxiv.org/abs/2608.03625  
- PDF: https://arxiv.org/pdf/2608.03625  
- arXiv DOI: https://doi.org/10.48550/arXiv.2608.03625  
- Ancillary CSV tables on arXiv (same id) under **Ancillary files**  
- This repository is the living code; **v0.1.0** is the release used for the paper tables
- **v0.2.0** is a performance release: same algorithms and intended table values, faster prime-field and large-`q` scans (see [CHANGELOG](CHANGELOG.md))

To cite the paper and/or this software, see [`CITATION.cff`](CITATION.cff) or:

```bibtex
@misc{linzi2026quotients,
  title         = {Finite-field {K}rasner quotients: isomorphism thresholds,
                   characteristics, and censuses},
  author        = {Linzi, Alessandro},
  year          = {2026},
  eprint        = {2608.03625},
  archivePrefix = {arXiv},
  primaryClass  = {math.RA},
  doi           = {10.48550/arXiv.2608.03625},
  url           = {https://arxiv.org/abs/2608.03625}
}
```

## Features

- Construct $K$ from $(p,k,d)$, $(q,d)$, or $(q,r)$
- Hyperaddition and full additive tables
- **Fast backends (v0.2)**
  - **Prime fields** $\mathbb{F}_p$: modular integer arithmetic (no per-object `galois` field)
  - **Extension fields** $\mathbb{F}_{p^k}$: `galois` with precomputed cosets
  - Large-$q$ **stable char / C-char** shortcuts (paper theorem; optional iterative check)
  - Scan-level **fingerprint cache** for Baker–Jin residues (atlases / $N_r^{\mathrm{emp}}$)
  - **Bitmask** encoding of table cells (Talotti-style; see below)
- **Layered isomorphism**
  - **Baker–Jin** (Proc. AMS 2021, Thm 1.1) for large $q$, with **Remark 1.2** bound $N_r$ (sharp Weil threshold; not the coarser $r^4$)
  - **General** $\mathrm{Aut}(C_r)$ check of $1 \boxplus x$ tables (gold standard / sporadics)
  - **Auto** policy combining both
- **Characteristic** and **C-characteristic** (Kędzierski–Linzi–Stojałowska)
- Modest census: list triples and classify iso classes for fixed order $n$ (fingerprint buckets by default)
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
| 1 | Stable (char, C-char) formula for large $\mathbb{F}_q/G_r$ | Baker–Jin + Linzi |
| 2 | Order-7 finite-field quotient atlas vs $H_7=277$ | Massouros |
| 3 | Ameri $H_n$ vs $Q^{\mathrm{fin}}$ for $n\le 6$ | Ameri |
| 4 | Massouros Props 1–2 sum-cardinality checks | Massouros |
| 5 | Linzi Prop.7 / Cor.1 / Prop.9 char bounds | Linzi et al. |

Module docs: `quotient_hyperfields/papers_experiments.py` (header map).

### Extended experiments (beyond arXiv tables)

```bash
qh-extended
# or
QH_R_MAX=12 QH_MAX_Q_CAP=15000 qh-extended
```

Pushes empirical \(N_r\), \(Q_r^{\mathrm{fin}}\), and order atlases **past** the paper range (\(r\le 8\), \(n\le 7\)). Writes checkpoints and CSVs under `experiments_output/` (see `experiments_output/README.md`). This is the scientific payoff of the v0.2 speedups.

## Mathematical references (in `papers/`)

| Paper | Role |
|-------|------|
| Baker–Jin, Proc. AMS 149 (2021) | Iso criterion for large finite-field quotients (Remark 1.2 $N_r$; safe $r^4$) |
| Kędzierski–Linzi–Stojałowska, Mathematics 11 (2023) | Char / C-char definitions and motivation |
| Ameri–Eyvazi–Hošková-Mayerová, AIMS Math. 5 (2020) | Enumeration order ≤ 6 |
| Massouros–Massouros, AIMS Math. 10 (2025) | Full hyperfields of order 7 (277 classes) |
| Talotti, *Hyperstructures in computer science* (2025) | Bitstring / integer encoding of finite magma table cells (implementation idea for `bitsets`; not a source of hyperfield theorems) |

Local PDF: [`papers/hyperstructure_comp_science.pdf`](papers/hyperstructure_comp_science.pdf). Related code by Talotti: [github.com/enh11/hyperstructures](https://github.com/enh11/hyperstructures).

## Performance (v0.2)

The companion paper freezes the **mathematical** procedures (coset hyperaddition, Baker–Jin, $\mathrm{Aut}(C_r)$ fingerprints, char via $n\times[1]$). Version **0.2.0** only optimises **how** those procedures run:

| Lever | Effect |
|-------|--------|
| Prime-field backend | Much faster construct + hyperadd for $\mathbb{F}_p$ (bulk of scans) |
| Precomputed cosets | Fewer multiplications per $\boxplus$ on extension fields |
| Stable char shortcut | $O(1)$ char/C-char for large $q$ (theorem of the paper) |
| Residue fingerprint cache | Avoids rebuilding $1\boxplus x$ for every large $q$ in a scan |
| Fingerprint census | $O(m)$ class partition instead of pairwise iso |

Bitmask table cells follow the representation style of Talotti (subsets as integers); the algebra remains finite-field Krasner quotients only.

Re-run experiments after upgrading:

```bash
qh-open-questions
qh-invariants
qh-papers
```

## Project layout

```
quotient_hyperfields/   # installable package
  hyperfield.py         # QuotientHyperfield (prime / extension backends)
  bitsets.py            # Talotti-style bitmask encodings of table cells
  isomorphism.py        # Baker–Jin / general / auto
  census.py             # triples + classification
  experiments.py        # open-question probes (Q3 / Q_r^fin)
  invariants_experiments.py  # char / C-char experiments A/B/D
  papers_experiments.py # paper tracks 1–5 (formula, atlas, criteria)
  extended_experiments.py # beyond-paper scans (N_r, Q_fin, atlases)
  criteria.py           # Massouros / Linzi checks
  atlas.py              # finite-field quotient class atlas
  literature.py         # Ameri / Massouros reference counts
  primes.py             # primality, prime powers, primitive roots
  cli.py                # thin interactive UI
  visualization.py      # census plots
tests/                  # paper-backed tests
papers/                 # reference PDFs
CHANGELOG.md            # version history
main.py                 # CLI entry (compat)
```

## Tests

```bash
pytest -q
```

## Scope

**In scope:** finite-field quotient hyperfields only; laptop-scale experiments; certified iso; char/C-char; empirical $N_r$ / $Q_r^{\mathrm{fin}}$ probes.

**Out of scope:** full enumeration of all finite hyperfields; general magma/hypergroup CAS; GUIs; non-Python cores.

Paper tables: use **v0.1.0** for byte-for-byte historical reproduction; **v0.2.0+** for faster re-runs with the same algorithms.

## AI use disclosure

Development of this repository was assisted by **Grok Build CLI** with the model **Grok 4.5** (xAI), primarily for software engineering: implementation, tests, packaging, documentation, and repository maintenance.

- Generative AI tools are **not** authors of the companion paper or of this software.
- Mathematical claims, proofs, experimental design, numerical tables, and released code were reviewed and verified by the human author (**Alessandro Linzi**), who remains solely responsible for correctness and for any errors.
- The peer-facing disclosure also appears in the **Acknowledgements** of [arXiv:2608.03625](https://arxiv.org/abs/2608.03625).

This wording follows common academic practice (name the tool and model; state the role; affirm human responsibility; do not list AI as an author).
