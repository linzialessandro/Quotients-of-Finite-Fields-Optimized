"""
Isomorphism checking functions for hyperfields.

This module contains functions to check if two quotient hyperfields are isomorphic.
It includes both a general method and an optimized criterion.
"""

import numpy as np

from hyperfield_modules.core import (
    build_field_and_subgroup,
    gf_coset_representatives,
    build_element_to_coset_index,
    gf_hyperaddition,
    gf_index_of_one,
)

def are_isomorphic_general(p1: int, k1: int, d1: int, p2: int, k2: int, d2: int) -> bool:
    """
    Check if GF(p1^k1)/G_d1 is isomorphic to GF(p2^k2)/G_d2 using the general method
    of comparing the multisets of sums involving the element [1].

    Args:
        p1, p2: Prime numbers
        k1, k2: Positive integers (field extension degrees)
        d1, d2: Positive integers dividing (p1^k1 - 1) and (p2^k2 - 1) respectively

    Returns:
        True if the hyperfields are isomorphic, False otherwise.
    """
    try:
        GF1, subgroup_keys1 = build_field_and_subgroup(p1, k1, d1)
        reps1 = gf_coset_representatives(GF1, subgroup_keys1)
        m1 = len(reps1)

        GF2, subgroup_keys2 = build_field_and_subgroup(p2, k2, d2)
        reps2 = gf_coset_representatives(GF2, subgroup_keys2)
        m2 = len(reps2)
    except Exception:
        return False

    if m1 != m2:
        return False

    element_to_idx1 = build_element_to_coset_index(GF1, subgroup_keys1, reps1)
    element_to_idx2 = build_element_to_coset_index(GF2, subgroup_keys2, reps2)

    table1_indices = [[set() for _ in range(m1 + 1)] for _ in range(m1 + 1)]
    for i in range(m1 + 1):
        table1_indices[0][i] = {i - 1} if i > 0 else {-1}
        table1_indices[i][0] = {i - 1} if i > 0 else {-1}
    for i in range(m1):
        for j in range(m1):
            table1_indices[i + 1][j + 1] = gf_hyperaddition(GF1, subgroup_keys1, reps1, element_to_idx1, i, j)

    table2_indices = [[set() for _ in range(m2 + 1)] for _ in range(m2 + 1)]
    for i in range(m2 + 1):
        table2_indices[0][i] = {i - 1} if i > 0 else {-1}
        table2_indices[i][0] = {i - 1} if i > 0 else {-1}
    for i in range(m2):
        for j in range(m2):
            table2_indices[i + 1][j + 1] = gf_hyperaddition(GF2, subgroup_keys2, reps2, element_to_idx2, i, j)

    try:
        idx_1_tab1 = gf_index_of_one(GF1, subgroup_keys1, reps1) + 1
        idx_1_tab2 = gf_index_of_one(GF2, subgroup_keys2, reps2) + 1
    except ValueError:
        return False

    inv_1_tab_idx1 = -1
    for j in range(m1 + 1):
        if -1 in table1_indices[idx_1_tab1][j]:
            inv_1_tab_idx1 = j
            break

    inv_1_tab_idx2 = -1
    for j in range(m2 + 1):
        if -1 in table2_indices[idx_1_tab2][j]:
            inv_1_tab_idx2 = j
            break

    if inv_1_tab_idx1 == -1 or inv_1_tab_idx2 == -1:
        return False

    if inv_1_tab_idx1 != inv_1_tab_idx2:
        return False

    multiset1 = []
    for i in range(m1 + 1):
        multiset1.append(frozenset(table1_indices[idx_1_tab1][i]))

    multiset2 = []
    for i in range(m2 + 1):
        multiset2.append(frozenset(table2_indices[idx_1_tab2][i]))

    multiset1.sort()
    multiset2.sort()

    return multiset1 == multiset2

def are_isomorphic_optimized(q1: int, d1: int, q2: int, d2: int) -> bool:
    """
    Check if GF(q1)/G_d1 is isomorphic to GF(q2)/G_d2 using an optimized criterion.
    Falls back to the general method for sporadic cases.

    Args:
        q1, q2: Orders of the finite fields (prime powers).
        d1, d2: Orders of the multiplicative subgroups.

    Returns:
        True if the hyperfields are isomorphic, False otherwise.
    """
    if (q1 - 1) % d1 != 0 or (q2 - 1) % d2 != 0:
        return False

    r1 = (q1 - 1) // d1
    r2 = (q2 - 1) // d2

    if r1 != r2:
        return False

    r = r1

    def get_p_k(q):
        if q < 2: return None, None
        p = q
        k = 1
        for i in range(2, int(q**0.5) + 1):
            if q % i == 0:
                p = i
                break
        if p == q:
            return q, 1
        k = round(np.log(q) / np.log(p))
        if p**k != q:
            return None, None
        return p, k

    p1, k1 = get_p_k(q1)
    p2, k2 = get_p_k(q2)

    if p1 is None or p2 is None:
        return False

    if q1 < r**4 or q2 < r**4:
        print(f"Sporadic case detected (q1={q1}, q2={q2}, r={r}). Falling back to general isomorphism check.")
        return are_isomorphic_general(p1, k1, d1, p2, k2, d2)

    if r % 2 != 0:
        return True
    else:
        return (q1 % (2 * r)) == (q2 % (2 * r))
