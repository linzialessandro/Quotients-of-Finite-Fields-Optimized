"""Massouros and Linzi paper criteria on quotient hyperfields."""

from quotient_hyperfields.criteria import (
    linzi_corollary1_char2_iff_even_d,
    linzi_prop7_char_bound,
    massouros_prop1_sum_bound,
    verify_all_paper_criteria,
)
from quotient_hyperfields.hyperfield import QuotientHyperfield
from quotient_hyperfields.papers_experiments import expected_stable_invariants


def test_massouros_prop1_on_samples():
    for q, r in [(17, 2), (19, 2), (31, 3), (97, 6), (5, 4)]:
        h = QuotientHyperfield.from_q_r(q, r)
        assert massouros_prop1_sum_bound(h).passed


def test_linzi_prop7_field():
    h = QuotientHyperfield.from_params(7, 1, 1)  # F_7
    assert linzi_prop7_char_bound(h).passed
    assert h.char == 7


def test_linzi_cor1_odd_p():
    # F_17, r=2: d=8 even ⇒ char should be 2
    h = QuotientHyperfield.from_q_r(17, 2)
    assert h.p == 17
    assert h.d % 2 == 0
    assert h.char == 2.0
    assert linzi_corollary1_char2_iff_even_d(h).passed

    # F_19, r=2: d=9 odd ⇒ char should not be 2 (we expect 3)
    h2 = QuotientHyperfield.from_q_r(19, 2)
    assert h2.d % 2 == 1
    assert h2.char != 2.0
    assert linzi_corollary1_char2_iff_even_d(h2).passed


def test_verify_all_on_large_odd_r():
    h = QuotientHyperfield.from_q_r(31, 3)
    results = verify_all_paper_criteria(h)
    assert all(r.passed for r in results)


def test_expected_stable_formula():
    assert expected_stable_invariants(3, 31) == (2.0, 1.0)
    assert expected_stable_invariants(2, 17) == (2.0, 1.0)  # 17 % 4 == 1
    assert expected_stable_invariants(2, 19) == (3.0, 1.0)  # 19 % 4 == 3
