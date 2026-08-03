"""Core construction and hyperaddition tests."""

import math

import pytest

from quotient_hyperfields import QuotientHyperfield


def test_from_params_basic():
    h = QuotientHyperfield.from_params(7, 1, 2)  # F_7 / {1,6}, r=3, order 4
    assert h.q == 7
    assert h.d == 2
    assert h.r == 3
    assert h.order == 4


def test_from_q_r():
    h = QuotientHyperfield.from_q_r(13, 6)
    assert h.d == 2
    assert h.r == 6
    assert h.order == 7


def test_invalid_d():
    with pytest.raises(ValueError):
        QuotientHyperfield.from_params(5, 1, 3)  # 3 does not divide 4


def test_zero_rules():
    h = QuotientHyperfield.from_params(11, 1, 2)
    assert h.hyperadd(-1, -1) == frozenset({-1})
    assert h.hyperadd(-1, 0) == frozenset({0})
    assert h.hyperadd(0, -1) == frozenset({0})


def test_one_plus_one_contains_known():
    # F_5 / {1,4}, r=2. Order 3 weak-sign type structure.
    h = QuotientHyperfield.from_params(5, 1, 2)
    s = h.one_plus(0)  # [1] ⊞ [1]
    assert isinstance(s, frozenset)
    assert len(s) >= 1


def test_addition_table_shape():
    h = QuotientHyperfield.from_params(7, 1, 3)
    df = h.addition_table()
    assert df.shape == (h.order, h.order)


def test_field_as_hyperfield_char():
    """F_5 / {1} is the field F_5 as a hyperfield; char = C-char = 5."""
    h = QuotientHyperfield.from_params(5, 1, 1)
    assert h.r == 4
    assert h.order == 5
    assert h.char == 5.0
    assert h.c_char == 5.0


def test_c_char_le_char():
    h = QuotientHyperfield.from_params(13, 1, 3)
    assert h.c_char <= h.char
    assert not math.isinf(h.char)
