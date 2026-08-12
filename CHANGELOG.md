# Changelog

## 0.2.0 — 2026-08-12

Performance release for the computational laboratory of arXiv:2608.03625.
**Mathematical definitions and paper algorithms are unchanged**; results of the
v0.1.0 / arXiv table suite are intended to be identical under the same scan
limits. Frozen release **v0.1.0** remains the artifact cited in the preprint.

### Faster core

- **Prime-field backend** (`k = 1`): construct and hyperadd with modular
  integers and a primitive root — no `galois` field object per prime `q`.
- **Extension backend** (`k > 1`): still uses `galois`, but cosets are
  precomputed so each hyperaddition is only `|G|^2` field additions.
- **Stable char / C-char shortcut**: when `q` is in the Baker–Jin large regime,
  `char` / `c_char` use the stable-characteristic theorem (default
  `use_stable=True`); pass `use_stable=False` for the definitional loop.
- **Scan fingerprint cache**: atlases, empirical `N_r`, and `Q_r^{fin}` reuse
  one Aut-canonical `1 ⊞ x` fingerprint per Baker–Jin residue for large `q`.
- **Census classification**: default path buckets by structure fingerprint
  (linear in the number of triples); pairwise iso remains available.
- **Bitmask cells** (`bitsets` module, Talotti-style): optional integer encoding
  of subset-valued table entries; used in the general isomorphism check.

### Documentation

- README: v0.2 performance notes and Talotti citation.
- Module docs describe backends and shortcuts.
- `CITATION.cff` version bump.

### Tests

- `tests/test_performance_v2.py`: backend selection, stable vs iterative char,
  fingerprint cache agreement, census fingerprint vs pairwise, regression on
  Massouros pairs.

## 0.1.0 — 2026-08-04

Initial public release accompanying arXiv:2608.03625 (tables and ancillary
software snapshot).
