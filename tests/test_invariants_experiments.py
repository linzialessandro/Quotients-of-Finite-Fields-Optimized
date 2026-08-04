"""Tests for char / C-char experiments A, B, D."""

import math

from quotient_hyperfields.hyperfield import QuotientHyperfield
from quotient_hyperfields.invariants_experiments import (
    gap_census_for_r,
    parity_char_check_for_r,
    stable_invariants_for_r,
)


def test_stable_invariants_r2_separates_classes():
    """r=2 large classes: char 2 vs 3, both C-char 1 (observed pattern)."""
    rows = stable_invariants_for_r(2, max_q_cap=200)
    by_res = {row.residue: row for row in rows}
    assert "even_0" in by_res and "even_1" in by_res
    assert by_res["even_0"].consistent and by_res["even_1"].consistent
    # Different chars separate the two Baker–Jin classes
    assert by_res["even_0"].char != by_res["even_1"].char
    assert by_res["even_0"].c_char == 1.0
    assert by_res["even_1"].c_char == 1.0


def test_stable_invariants_r3_consistent():
    rows = stable_invariants_for_r(3, max_q_cap=200)
    assert len(rows) == 1
    assert rows[0].residue == "odd"
    assert rows[0].consistent
    assert rows[0].c_char <= rows[0].char


def test_gap_census_r2():
    c = gap_census_for_r(2, max_q_cap=200)
    assert c.n_samples >= 5
    assert c.equal_count + sum(
        n for g, n in c.gap_counts.items() if g and g > 0
    ) <= c.n_samples + c.equal_count
    # Some order-3 quotients have C-char 1
    assert c.cchar_one_count >= 1


def test_parity_odd_r_char_2():
    res = parity_char_check_for_r(3, max_q_cap=200)
    assert res.order == 4  # even
    assert res.passed
    assert res.n_checked >= 1


def test_field_as_hyperfield_equal_invariants():
    h = QuotientHyperfield.from_params(7, 1, 1)
    assert h.char == h.c_char == 7.0
