"""Paper-backed isomorphism tests (Baker–Jin, Massouros pairs)."""

import pytest

from quotient_hyperfields import (
    QuotientHyperfield,
    are_isomorphic,
    are_isomorphic_baker_jin,
    are_isomorphic_general,
)
from quotient_hyperfields.isomorphism import baker_jin_bound


def test_different_r_not_iso():
    h1 = QuotientHyperfield.from_q_r(13, 3)
    h2 = QuotientHyperfield.from_q_r(13, 6)
    r = are_isomorphic(h1, h2)
    assert r.isomorphic is False


def test_baker_jin_odd_r_large():
    # r=3 odd, N_3 = 81 with r^4; use q >= 81
    # F_81 = 3^4, r=3 | 80? 80/3 no. Need q ≡ 1 (mod 3), q >= 81.
    # F_103: 103-1=102, 102/3=34, d=34. 103 > 81.
    # F_109: 108/3=36.
    h1 = QuotientHyperfield.from_q_r(103, 3)
    h2 = QuotientHyperfield.from_q_r(109, 3)
    bj = are_isomorphic_baker_jin(h1, h2)
    assert bj.sporadic is False
    assert bj.isomorphic is True
    gen = are_isomorphic_general(h1, h2)
    assert gen.isomorphic is True
    auto = are_isomorphic(h1, h2, method="auto")
    assert auto.isomorphic is True


def test_baker_jin_even_r_same_class():
    # r=2, N_2 = 16. q ≡ 1 mod 4 vs q ≡ 3 mod 4 are the two classes.
    # F_17: 17 ≡ 1 mod 4; F_29: 29 ≡ 1 mod 4. Both large, same class.
    h1 = QuotientHyperfield.from_q_r(17, 2)
    h2 = QuotientHyperfield.from_q_r(29, 2)
    bj = are_isomorphic_baker_jin(h1, h2)
    assert bj.isomorphic is True
    assert are_isomorphic_general(h1, h2).isomorphic is True


def test_baker_jin_even_r_different_class():
    # F_17 ≡ 1 mod 4, F_19 ≡ 3 mod 4, r=2, both > 16
    h1 = QuotientHyperfield.from_q_r(17, 2)
    h2 = QuotientHyperfield.from_q_r(19, 2)
    bj = are_isomorphic_baker_jin(h1, h2)
    assert bj.isomorphic is False
    assert are_isomorphic_general(h1, h2).isomorphic is False


def test_sporadic_uses_general():
    # r=3, N_r=81; F_7 and F_13 are sporadic (q < 81)
    h1 = QuotientHyperfield.from_q_r(7, 3)
    h2 = QuotientHyperfield.from_q_r(13, 3)
    bj = are_isomorphic_baker_jin(h1, h2)
    assert bj.sporadic is True
    assert bj.isomorphic is None
    auto = are_isomorphic(h1, h2, method="auto")
    assert auto.sporadic is True
    assert auto.isomorphic is not None  # general decided


def test_self_isomorphic():
    h = QuotientHyperfield.from_params(11, 1, 5)
    assert are_isomorphic(h, h).isomorphic is True


def test_massouros_order7_z97_z157():
    """
    Massouros (arXiv:2412.11331): Z_97/G ≅ Z_157/G for index-6 subgroups.
    r=6, N_6 = 6^4 = 1296, so both q are sporadic under r^4 and use general.
    """
    h1 = QuotientHyperfield.from_q_r(97, 6)
    h2 = QuotientHyperfield.from_q_r(157, 6)
    assert h1.order == 7 and h2.order == 7
    r = are_isomorphic(h1, h2, method="general")
    assert r.isomorphic is True


def test_massouros_order7_distinct_classes():
    """Z_13/G and Z_19/G (index 6) are different classes in Massouros tables."""
    h1 = QuotientHyperfield.from_q_r(13, 6)
    h2 = QuotientHyperfield.from_q_r(19, 6)
    r = are_isomorphic(h1, h2, method="general")
    assert r.isomorphic is False


def test_bound_r4():
    assert baker_jin_bound(2) == 16
    assert baker_jin_bound(3) == 81
    assert baker_jin_bound(6) == 1296


def test_auto_matches_general_on_medium():
    pairs = [
        ((11, 2), (19, 2)),
        ((13, 3), (19, 3)),
        ((17, 4), (41, 4)),
    ]
    for (q1, r), (q2, r2) in pairs:
        assert r == r2
        h1 = QuotientHyperfield.from_q_r(q1, r)
        h2 = QuotientHyperfield.from_q_r(q2, r)
        auto = are_isomorphic(h1, h2, method="auto")
        gen = are_isomorphic_general(h1, h2)
        assert auto.isomorphic == gen.isomorphic
