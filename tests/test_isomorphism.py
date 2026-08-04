"""Paper-backed isomorphism tests (Baker–Jin, Massouros pairs)."""

from quotient_hyperfields import (
    QuotientHyperfield,
    are_isomorphic,
    are_isomorphic_baker_jin,
    are_isomorphic_general,
)
from quotient_hyperfields.isomorphism import (
    baker_jin_bound,
    baker_jin_bound_safe,
    baker_jin_bound_sharp,
    is_baker_jin_large,
    weil_excess,
)


def test_different_r_not_iso():
    h1 = QuotientHyperfield.from_q_r(13, 3)
    h2 = QuotientHyperfield.from_q_r(13, 6)
    r = are_isomorphic(h1, h2)
    assert r.isomorphic is False


def test_remark_12_paper_values():
    """Baker–Jin Remark 1.2: N_2 = 6, N_3 = 17."""
    assert baker_jin_bound_sharp(2) == 6
    assert baker_jin_bound_sharp(3) == 17
    assert baker_jin_bound(2) == 6
    assert baker_jin_bound(3) == 17

    assert weil_excess(5, 2) <= 0
    assert weil_excess(6, 2) > 0
    assert weil_excess(16, 3) <= 0
    assert weil_excess(17, 3) > 0


def test_remark_12_threshold_monotonic():
    """Once q is large for index r, every larger order stays large."""
    for r in range(2, 12):
        n_r = baker_jin_bound_sharp(r)
        assert weil_excess(n_r, r) > 0
        if n_r > 1:
            assert weil_excess(n_r - 1, r) <= 0
        assert is_baker_jin_large(n_r, r)
        assert not is_baker_jin_large(n_r - 1, r) if n_r > 1 else True
        # safe bound dominates
        assert n_r <= baker_jin_bound_safe(r)
        assert baker_jin_bound_safe(r) == r**4


def test_baker_jin_odd_r_large():
    # r=3 odd, N_3 = 17 (Remark 1.2). F_19 and F_31 both >= 17.
    h1 = QuotientHyperfield.from_q_r(19, 3)
    h2 = QuotientHyperfield.from_q_r(31, 3)
    assert is_baker_jin_large(19, 3) and is_baker_jin_large(31, 3)
    bj = are_isomorphic_baker_jin(h1, h2)
    assert bj.sporadic is False
    assert bj.isomorphic is True
    gen = are_isomorphic_general(h1, h2)
    assert gen.isomorphic is True
    auto = are_isomorphic(h1, h2, method="auto")
    assert auto.isomorphic is True


def test_baker_jin_even_r_same_class():
    # r=2, N_2 = 6. F_13 ≡ 1 mod 4, F_17 ≡ 1 mod 4.
    h1 = QuotientHyperfield.from_q_r(13, 2)
    h2 = QuotientHyperfield.from_q_r(17, 2)
    bj = are_isomorphic_baker_jin(h1, h2)
    assert bj.sporadic is False
    assert bj.isomorphic is True
    assert are_isomorphic_general(h1, h2).isomorphic is True


def test_baker_jin_even_r_different_class():
    # F_13 ≡ 1 mod 4, F_11 ≡ 3 mod 4, r=2, both >= 6
    h1 = QuotientHyperfield.from_q_r(13, 2)
    h2 = QuotientHyperfield.from_q_r(11, 2)
    bj = are_isomorphic_baker_jin(h1, h2)
    assert bj.sporadic is False
    assert bj.isomorphic is False
    assert are_isomorphic_general(h1, h2).isomorphic is False


def test_sporadic_uses_general():
    # r=3, N_r=17; F_7 and F_13 are sporadic
    h1 = QuotientHyperfield.from_q_r(7, 3)
    h2 = QuotientHyperfield.from_q_r(13, 3)
    assert not is_baker_jin_large(7, 3)
    assert not is_baker_jin_large(13, 3)
    bj = are_isomorphic_baker_jin(h1, h2)
    assert bj.sporadic is True
    assert bj.isomorphic is None
    auto = are_isomorphic(h1, h2, method="auto")
    assert auto.sporadic is True
    assert auto.isomorphic is not None  # general decided


def test_safe_bound_override():
    """Forcing r^4 marks some Remark-1.2-large fields as sporadic."""
    # r=3: 19 >= 17 (sharp) but 19 < 81 (safe)
    h1 = QuotientHyperfield.from_q_r(19, 3)
    h2 = QuotientHyperfield.from_q_r(31, 3)
    sharp = are_isomorphic_baker_jin(h1, h2)
    assert sharp.sporadic is False
    safe = are_isomorphic_baker_jin(h1, h2, n_r=baker_jin_bound_safe(3))
    assert safe.sporadic is True
    assert safe.isomorphic is None


def test_self_isomorphic():
    h = QuotientHyperfield.from_params(11, 1, 5)
    assert are_isomorphic(h, h).isomorphic is True


def test_massouros_order7_z97_z157():
    """
    Massouros–Massouros (AIMS Math. 2025): Z_97/G ≅ Z_157/G for index-6 subgroups.
    r=6, N_6 (Remark 1.2) is still > 157, so general path applies.
    """
    n6 = baker_jin_bound_sharp(6)
    assert n6 > 157
    h1 = QuotientHyperfield.from_q_r(97, 6)
    h2 = QuotientHyperfield.from_q_r(157, 6)
    assert h1.order == 7 and h2.order == 7
    bj = are_isomorphic_baker_jin(h1, h2)
    assert bj.sporadic is True
    r = are_isomorphic(h1, h2, method="general")
    assert r.isomorphic is True
    assert are_isomorphic(h1, h2, method="auto").isomorphic is True


def test_massouros_order7_distinct_classes():
    """Z_13/G and Z_19/G (index 6) are different classes in Massouros tables."""
    h1 = QuotientHyperfield.from_q_r(13, 6)
    h2 = QuotientHyperfield.from_q_r(19, 6)
    r = are_isomorphic(h1, h2, method="general")
    assert r.isomorphic is False


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
