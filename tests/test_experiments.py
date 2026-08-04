"""Tests for Baker–Jin open-question experiments."""

from quotient_hyperfields.experiments import (
    baker_jin_residue,
    count_finite_quotient_classes,
    empirical_nr,
    has_large_witnesses,
    paper_lower_bound_nr,
    prime_powers_with_index,
    resolve_scan_limit,
    structure_fingerprint,
)
from quotient_hyperfields.hyperfield import QuotientHyperfield
from quotient_hyperfields.isomorphism import baker_jin_bound_sharp


def test_prime_powers_with_index_r2():
    qs = prime_powers_with_index(2, 30)
    assert 3 in qs and 9 in qs and 25 in qs
    assert all((q - 1) % 2 == 0 for q in qs)


def test_fingerprint_iso_invariant():
    h1 = QuotientHyperfield.from_q_r(13, 2)
    h2 = QuotientHyperfield.from_q_r(17, 2)
    assert structure_fingerprint(h1) == structure_fingerprint(h2)


def test_fingerprint_distinguishes_even_classes():
    h0 = QuotientHyperfield.from_q_r(13, 2)  # 13 ≡ 1 mod 4
    h1 = QuotientHyperfield.from_q_r(11, 2)  # 11 ≡ 3 mod 4
    assert structure_fingerprint(h0) != structure_fingerprint(h1)


def test_paper_lower_bound():
    assert paper_lower_bound_nr(3) == 2**2 + 1
    assert paper_lower_bound_nr(4) == 3**2 + 1
    assert paper_lower_bound_nr(6) == 5**2 + 1
    assert paper_lower_bound_nr(5) is None


def test_resolve_scan_limit_has_witnesses():
    for r in (2, 3, 4):
        limit = resolve_scan_limit(r, max_q_cap=500)
        qs = prime_powers_with_index(r, limit)
        assert has_large_witnesses(r, qs)
        assert limit >= baker_jin_bound_sharp(r)


def test_empirical_nr_r2():
    """Paper: Remark 1.2 bound is sharp for r=2 (N_2=6); sporadics q=3,5."""
    res = empirical_nr(2, max_q_cap=200)
    assert res.n_r_remark_12 == 6
    assert res.n_r_empirical is not None
    # Exceptions at small q should force empirical near 6
    assert res.n_r_empirical >= 6
    assert res.n_r_empirical <= 16  # well below safe bound
    assert has_large_witnesses(2, res.prime_powers)


def test_empirical_nr_r3():
    """Paper: Remark 1.2 bound sharp for r=3 (N_3=17)."""
    res = empirical_nr(3, max_q_cap=200)
    assert res.n_r_remark_12 == 17
    assert res.n_r_empirical is not None
    assert res.n_r_empirical <= 17
    # Should see some exceptions below the large regime
    assert any(q < 17 for q, _ in res.exceptions)


def test_baker_jin_residue():
    assert baker_jin_residue(17, 2) == "even_0"
    assert baker_jin_residue(19, 2) == "even_1"
    assert baker_jin_residue(19, 3) == "odd"


def test_count_qr_fin_r2():
    limit = resolve_scan_limit(2, max_q_cap=200)
    c = count_finite_quotient_classes(2, limit)
    assert c.theoretically_complete
    assert c.q_r_fin >= 2
