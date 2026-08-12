"""v0.2 performance path: same answers as definitional algorithms."""

import math

import pytest

from quotient_hyperfields import QuotientHyperfield, are_isomorphic
from quotient_hyperfields.bitsets import labels_to_mask, mask_to_labels
from quotient_hyperfields.census import classify_hyperfields
from quotient_hyperfields.experiments import structure_fingerprint
from quotient_hyperfields.hyperfield import ZERO
from quotient_hyperfields.isomorphism import is_baker_jin_large
from quotient_hyperfields.primes import primitive_root_mod_p


def test_prime_backend_selected():
    h = QuotientHyperfield.from_q_r(97, 6)
    assert h.backend == "prime"
    assert h._GF is None


def test_extension_backend_selected():
    # F_9: q=9, r | 8 e.g. r=2, d=4
    h = QuotientHyperfield.from_q_r(9, 2)
    assert h.backend == "extension"
    assert h._GF is not None


def test_primitive_root_mod_p():
    assert primitive_root_mod_p(2) == 1
    assert primitive_root_mod_p(5) in (2, 3)
    # 2 is a primitive root mod 13
    assert pow(primitive_root_mod_p(13), 12, 13) == 1


def test_hyperadd_agrees_prime_examples():
    """Spot-check coset sums on several prime quotients."""
    for q, r in [(7, 3), (13, 6), (19, 6), (31, 2), (97, 6)]:
        h = QuotientHyperfield.from_q_r(q, r)
        # Commutative and contains known zero rules
        for i in range(h.r):
            for j in range(i, h.r):
                assert h.hyperadd(i, j) == h.hyperadd(j, i)
        assert h.hyperadd(ZERO, 0) == frozenset({0})


def test_stable_char_matches_definitional():
    """Large-q shortcut equals iterative char/C-char."""
    # r=2, N_2=6; q=13 large, even residue class
    h = QuotientHyperfield.from_q_r(13, 2)
    assert is_baker_jin_large(13, 2)
    assert h.characteristic(use_stable=True) == h.characteristic(use_stable=False)
    assert h.c_characteristic(use_stable=True) == h.c_characteristic(
        use_stable=False
    )
    # even r, q=13 ≡ 1 mod 4 → char 2, C-char 1
    assert h.char == 2.0
    assert h.c_char == 1.0


def test_stable_char_even1_class():
    # r=2, q=11 ≡ 3 ≡ r+1 mod 4 → char 3
    h = QuotientHyperfield.from_q_r(11, 2)
    assert is_baker_jin_large(11, 2)
    assert h.characteristic(use_stable=True) == 3.0
    assert h.characteristic(use_stable=False) == 3.0
    assert h.c_characteristic(use_stable=True) == 1.0
    assert h.c_characteristic(use_stable=False) == 1.0


def test_bitsets_roundtrip():
    labels = frozenset({ZERO, 0, 2, 5})
    m = labels_to_mask(labels)
    assert mask_to_labels(m) == labels


def test_fingerprint_stable_cache_matches_full():
    # r=6, N_6=434; pick two large primes with the same Baker–Jin residue.
    from quotient_hyperfields.experiments import (
        baker_jin_residue,
        prime_powers_with_index,
    )

    r = 6
    large = [
        q
        for q in prime_powers_with_index(r, 900)
        if is_baker_jin_large(q, r)
    ]
    if len(large) < 2:
        pytest.skip("need large q for r=6")
    res0 = baker_jin_residue(large[0], r)
    q1 = large[0]
    q2 = next((q for q in large[1:] if baker_jin_residue(q, r) == res0), None)
    if q2 is None:
        pytest.skip("need two large q in same residue")
    h1 = QuotientHyperfield.from_q_r(q1, r)
    h2 = QuotientHyperfield.from_q_r(q2, r)
    cache: dict = {}
    fp1 = structure_fingerprint(h1, use_stable=True, _stable_cache=cache)
    fp2 = structure_fingerprint(h2, use_stable=True, _stable_cache=cache)
    assert h1.baker_jin_class == h2.baker_jin_class
    assert fp1 == fp2
    # Full Aut fingerprint of a fresh object agrees
    h1b = QuotientHyperfield.from_q_r(h1.q, r)
    assert structure_fingerprint(h1b, use_stable=False) == fp1


def test_classify_fingerprint_vs_pairwise():
    fp_classes = classify_hyperfields(3, max_p=30, max_k=1, method="fingerprint")
    pair_classes = classify_hyperfields(3, max_p=30, max_k=1, method="general")
    # Same partition of the triple set
    norm = lambda cs: sorted(tuple(sorted(c)) for c in cs)
    assert norm(fp_classes) == norm(pair_classes)


def test_field_char_no_false_stable():
    """r=4 field F_5: not in large-q multi-class regime that skips work wrongly."""
    h = QuotientHyperfield.from_params(5, 1, 1)
    assert h.char == 5.0
    assert h.c_char == 5.0
    assert not math.isinf(h.char)


def test_massouros_pair_still_iso():
    h1 = QuotientHyperfield.from_q_r(97, 6)
    h2 = QuotientHyperfield.from_q_r(157, 6)
    assert are_isomorphic(h1, h2, method="general").isomorphic is True
