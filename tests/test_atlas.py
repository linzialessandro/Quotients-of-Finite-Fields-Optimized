"""Finite-field quotient atlas tests."""

from quotient_hyperfields.atlas import build_order_atlas
from quotient_hyperfields.literature import AMERI_H_BY_ORDER, MASSOUROS_H7


def test_atlas_order3():
    # r=2, order 3: Baker–Jin says 2 stable + 2 sporadics (q=3,5) = 4 finite-field
    atlas = build_order_atlas(3, max_q_cap=200)
    assert atlas.theoretically_complete
    assert atlas.q_fin == 4
    assert atlas.q_fin <= AMERI_H_BY_ORDER[3]


def test_atlas_order4():
    atlas = build_order_atlas(4, max_q_cap=200)
    assert atlas.q_fin >= 1
    assert atlas.q_fin <= AMERI_H_BY_ORDER[4]


def test_massouros_ratio_bound():
    # Q_fin for order 7 should be far below 277
    atlas = build_order_atlas(7, max_q_cap=800)
    assert atlas.q_fin < MASSOUROS_H7
    assert atlas.q_fin >= 2  # at least two stable classes for r=6 even
