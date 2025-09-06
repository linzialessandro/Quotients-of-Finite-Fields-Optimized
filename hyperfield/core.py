"""
Core hyperfield functions.

This module contains the core logic for hyperfield generation, including building
the finite field and subgroup, computing coset representatives, and performing
hyperoperations.
"""

import galois
import pandas as pd

from utils import is_prime

def validate_inputs(p, k, d):
    """Validate the inputs for the hyperfield generation."""
    if not is_prime(p):
        raise ValueError("p must be prime")
    if k <= 0:
        raise ValueError("k must be a positive integer")
    if d <= 0:
        raise ValueError("d must be a positive integer")

def build_field_and_subgroup(p, k, d):
    """Build the finite field GF(p^k) and the multiplicative subgroup of order d."""
    if galois is None:
        raise ImportError("The 'galois' library is required. Install with 'pip install galois'.")
    validate_inputs(p, k, d)
    q = p ** k
    if (q - 1) % d != 0:
        raise ValueError(f"d must divide p^k - 1 = {q - 1}")
    GF = galois.GF(q)
    alpha = GF.primitive_element
    step = (q - 1) // d
    h = alpha ** step
    subgroup_keys = set(int(h ** i) for i in range(d))
    return GF, subgroup_keys

def gf_coset_representatives(GF, subgroup_keys):
    """Determine the coset representatives for GF(p^k)^*/G_d."""
    q = GF.order
    seen_keys = set()
    reps_keys = []
    order = list(range(1, q))
    order.sort(key=lambda t: 0 if t == 1 else 1)
    for key in order:
        if key in seen_keys:
            continue
        coset_keys = set(int(GF(key) * GF(gk)) for gk in subgroup_keys)
        if not any(ck in seen_keys for ck in coset_keys):
            reps_keys.append(key)
            seen_keys |= coset_keys
    return [GF(k) for k in reps_keys]

def build_element_to_coset_index(GF, subgroup_keys, reps):
    """Build a mapping from field element keys to coset indices."""
    mapping = {}
    for idx, r in enumerate(reps):
        r_key = int(r)
        for gk in subgroup_keys:
            elem_key = int(GF(r_key) * GF(gk))
            mapping[elem_key] = idx
    return mapping

def gf_hyperaddition(GF, subgroup_keys, reps, element_to_idx, i, j):
    """Perform hyperaddition of two cosets."""
    zero_key = 0
    r_i_key = int(reps[i])
    r_j_key = int(reps[j])
    A_keys = set(int(GF(r_i_key) * GF(gk)) for gk in subgroup_keys)
    B_keys = set(int(GF(r_j_key) * GF(gk)) for gk in subgroup_keys)
    result_indices = set()
    for a_key in A_keys:
        for b_key in B_keys:
            s_key = int(GF(a_key) + GF(b_key))
            if s_key == zero_key:
                result_indices.add(-1)
            else:
                result_indices.add(element_to_idx[s_key])
    return result_indices

def gf_index_of_one(GF, subgroup_keys, reps):
    """Get the coset index of the element [1]."""
    element_to_idx = build_element_to_coset_index(GF, subgroup_keys, reps)
    one_key = int(GF(1))
    if one_key not in element_to_idx:
        raise ValueError("Element [1] (coset containing 1) not found in representatives.")
    return element_to_idx[one_key]

def gf_hyper_sum_of_ones(GF, subgroup_keys, reps, n_terms):
    """Compute the hyper-sum of n_terms ones."""
    element_to_idx = build_element_to_coset_index(GF, subgroup_keys, reps)
    idx_1 = gf_index_of_one(GF, subgroup_keys, reps)
    S = {idx_1}
    for _ in range(n_terms - 1):
        T = set()
        for a in S:
            if a == -1:
                T.add(idx_1)  # 0 + [1] = [1]
            else:
                T |= gf_hyperaddition(GF, subgroup_keys, reps, element_to_idx, a, idx_1)
        S = T
    return S

def gf_characteristic(GF, subgroup_keys, reps, max_terms=2000):
    """Calculate the characteristic of the hyperfield."""
    for n in range(1, max_terms + 1):
        if -1 in gf_hyper_sum_of_ones(GF, subgroup_keys, reps, n):
            return n
    return 0

def gf_c_characteristic(GF, subgroup_keys, reps, max_terms=2000):
    """Calculate the c-characteristic of the hyperfield."""
    idx_1 = gf_index_of_one(GF, subgroup_keys, reps)
    for n in range(1, max_terms + 1):
        if idx_1 in gf_hyper_sum_of_ones(GF, subgroup_keys, reps, n + 1):
            return n
    return 0

def build_addition_table_gf(GF, subgroup_keys, reps):
    """Build the addition table for the hyperfield."""
    m = len(reps)
    element_to_idx = build_element_to_coset_index(GF, subgroup_keys, reps)
    table = [[set() for _ in range(m + 1)] for _ in range(m + 1)]
    rep_keys = [int(r) for r in reps]

    def label_of_key(k):
        if k == 0:
            return "0"
        if k == 1:
            return "1"
        return str(GF(k))

    headers_labels = ["0"] + [label_of_key(k) for k in rep_keys]
    for i in range(m + 1):
        for j in range(m + 1):
            if i == 0 and j == 0:
                table[i][j] = {0}
            elif i == 0:
                table[i][j] = {rep_keys[j - 1]}
            elif j == 0:
                table[i][j] = {rep_keys[i - 1]}
            else:
                result_indices = gf_hyperaddition(GF, subgroup_keys, reps, element_to_idx, i - 1, j - 1)
                result_keys = set(rep_keys[idx] if idx >= 0 else 0 for idx in result_indices)
                table[i][j] = result_keys

    df_table = pd.DataFrame(table, index=[f"[{h}]" for h in headers_labels], columns=[f"[{h}]" for h in headers_labels])

    def format_cell(cell_set):
        return "{" + ",".join([f"[{label_of_key(k)}]" for k in sorted(cell_set)]) + "}"

    df_table = df_table.map(format_cell)
    return df_table
