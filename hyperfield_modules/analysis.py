"""
Analysis pipeline for hyperfields.

This module contains functions for generating triples (p, k, d), classifying
hyperfields into isomorphism classes, and analyzing their characteristics.
"""

from utils import is_prime
from hyperfield_modules.isomorphism import are_isomorphic_optimized
from hyperfield_modules.core import (
    build_field_and_subgroup,
    gf_coset_representatives,
    gf_characteristic,
    gf_c_characteristic,
)

def generate_triples(n, max_p, max_k):
    """
    Generates valid triples (p, k, d) such that n = (p^k - 1) / d + 1,
    p is prime and <= max_p, and 1 <= k <= max_k.
    """
    valid_triples = []
    if n <= 1:
        print("Warning: n must be greater than 1 for a valid quotient hyperfield.")
        return valid_triples

    for p in range(2, max_p + 1):
        if is_prime(p):
            for k in range(1, max_k + 1):
                q = p ** k
                if (n - 1) > 0 and (q - 1) % (n - 1) == 0:
                    d = (q - 1) // (n - 1)
                    valid_triples.append((p, k, d))
    return valid_triples

def classify_hyperfields(n, max_p, max_k):
    """
    Generates and classifies hyperfields based on (p, k, d) triples
    into isomorphism classes using the optimized check.
    """
    triples = generate_triples(n, max_p, max_k)
    isomorphism_classes = []
    assigned_triples = set()

    q_d_triples = [(p**k, d) for p, k, d in triples]
    original_triples_map = {(p**k, d): (p, k, d) for p, k, d in triples}

    for q1, d1 in q_d_triples:
        original_triple1 = original_triples_map[(q1, d1)]
        if original_triple1 in assigned_triples:
            continue

        current_class = [original_triple1]
        assigned_triples.add(original_triple1)

        for q2, d2 in q_d_triples:
            original_triple2 = original_triples_map[(q2, d2)]
            if original_triple2 not in assigned_triples:
                if are_isomorphic_optimized(q1, d1, q2, d2):
                    current_class.append(original_triple2)
                    assigned_triples.add(original_triple2)

        isomorphism_classes.append(current_class)

    return isomorphism_classes

def analyze_characteristics(isomorphism_classes):
    """
    Calculates the characteristic and c-characteristic for each isomorphism class representative.
    """
    characteristics = []
    c_characteristics = []

    for iso_class in isomorphism_classes:
        p, k, d = iso_class[0]
        try:
            GF, subgroup_keys = build_field_and_subgroup(p, k, d)
            reps = gf_coset_representatives(GF, subgroup_keys)
            char = gf_characteristic(GF, subgroup_keys, reps)
            cchar = gf_c_characteristic(GF, subgroup_keys, reps)
            characteristics.append(char)
            c_characteristics.append(cchar)
        except Exception as e:
            print(f"Could not analyze characteristics for triple ({p}, {k}, {d}): {e}")
            characteristics.append(None)
            c_characteristics.append(None)

    return characteristics, c_characteristics
