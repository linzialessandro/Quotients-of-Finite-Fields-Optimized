"""
Reference data from the project papers (counts and stated facts).

Sources
-------
- Ameri–Eyvazi–Hoskova-Mayerova, AIMS Math. 2020, Table 1
  (enumeration of all finite hyperfields of order ≤ 6).
- Baker–Jin, Proc. Amer. Math. Soc. 149 (2021), Theorem 1.3
  (which hyperfields of order ≤ 4 are quotients of fields).
- Massouros–Massouros, AIMS Mathematics 10 (2025), 21287–21421
  (DOI 10.3934/math.2025951; 277 hyperfields of order 7).

These constants are used by :mod:`papers_experiments` for comparison tables.
They are *not* recomputed here; they are literature inputs.
"""

from __future__ import annotations

# Ameri et al., Table 1: number of isomorphism classes of hyperfields of order n
# (all Krasner hyperfields, not only quotients).
AMERI_H_BY_ORDER: dict[int, int] = {
    2: 2,
    3: 5,
    4: 7,
    5: 27,  # 11 + 16 by mult. group type in their table
    6: 16,
}

# Baker–Jin Theorem 1.3 (order ≤ 4): among all hyperfields,
# how many are quotients of *some* field (finite or infinite), when stated.
# Values below are for orientation; finite-field-only counts come from our atlas.
BAKER_JIN_ORDER_LE4_NOTES: dict[int, str] = {
    2: "2 classes total; both are quotients of finite fields (K and F_2).",
    3: "5 classes total; 4 are finite-field quotients; signs ≅ R/R>0 not finite-field.",
    4: "7 classes total; 4 finite-field quotients; 1 infinite-field only; 2 non-quotient.",
}

# Massouros order 7: total hyperfields up to isomorphism
MASSOUROS_H7: int = 277

# Massouros: index-6 quotients of prime fields they explicitly identify
# (non-exhaustive list of *examples* from the paper narrative).
# Format: (q, note)
MASSOUROS_ORDER7_QUOTIENT_EXAMPLES: list[tuple[int, str]] = [
    (19, "Z_19/G ≅ HF_7^9 in their labelling discussion"),
    (31, "Z_31/G"),
    (43, "Z_43/G"),
    (13, "Z_13/G"),
    (37, "Z_37/G"),
    (61, "Z_61/G"),
    (73, "Z_73/G"),
    (97, "Z_97/G ≅ Z_157/G"),
    (157, "Z_157/G ≅ Z_97/G"),
    (109, "Z_109/G"),
    (67, "Z_67/G ≅ Z_79/G ≅ Z_139/G"),
    (79, "Z_79/G"),
    (139, "Z_139/G"),
    (103, "Z_103/G and related"),
    (181, "Z_181/G"),
]

# Hyperfield order n ↔ subgroup index r = n - 1
def order_to_index(n: int) -> int:
    return n - 1


def index_to_order(r: int) -> int:
    return r + 1
