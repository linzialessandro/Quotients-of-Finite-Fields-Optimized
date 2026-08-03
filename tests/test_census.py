"""Census helpers."""

from quotient_hyperfields.census import classify_hyperfields, generate_triples


def test_generate_triples_order3():
    # order 3 => r=2
    triples = generate_triples(3, max_p=20, max_k=2)
    assert len(triples) >= 1
    for p, k, d in triples:
        q = p**k
        assert (q - 1) // d + 1 == 3


def test_classify_small():
    classes = classify_hyperfields(3, max_p=20, max_k=1)
    assert len(classes) >= 1
    # All triples appear exactly once
    triples = generate_triples(3, 20, 1)
    flat = [t for c in classes for t in c]
    assert sorted(flat) == sorted(triples)
