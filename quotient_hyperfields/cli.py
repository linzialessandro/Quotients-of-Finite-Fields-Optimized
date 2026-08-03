"""Thin interactive CLI over the research library."""

from __future__ import annotations

from quotient_hyperfields.census import (
    analyze_characteristics,
    classify_hyperfields,
    generate_triples,
)
from quotient_hyperfields.hyperfield import QuotientHyperfield
from quotient_hyperfields.isomorphism import are_isomorphic
from quotient_hyperfields.primes import is_prime
from quotient_hyperfields.visualization import visualize_characteristics


def _display(obj) -> None:
    try:
        from IPython.display import display as ipy_display

        ipy_display(obj)
    except Exception:
        try:
            print(obj.to_string())
        except Exception:
            print(obj)


def generate_and_display_table() -> None:
    print("\nGenerate Additive Table for a Hyperfield")
    print("=" * 40)
    try:
        p = int(input("Enter a prime number p: "))
        k = int(input("Enter an integer k >= 1: "))
        q = p**k
        print(f"p^k - 1 = {q - 1}")
        d = int(input(f"Enter divisor d of {q - 1}: "))
    except (ValueError, KeyboardInterrupt):
        print("\nCancelled.")
        return

    if not is_prime(p) or k < 1 or d < 1 or (q - 1) % d != 0:
        print("Invalid parameters.")
        return

    h = QuotientHyperfield.from_params(p, k, d)
    print(f"{h}: order {h.order}, char={h.char}, C-char={h.c_char}")
    print("\nAdditive hyperoperation table (+):\n")
    _display(h.addition_table())


def run_isomorphism_check() -> None:
    print("Isomorphism Check for Quotient Hyperfields F_q / G_d")
    print("=" * 55)
    try:
        print("\nFirst hyperfield:")
        p1 = int(input("  prime p1: "))
        k1 = int(input("  degree k1: "))
        d1 = int(input(f"  d1 | {p1**k1 - 1}: "))
        print("\nSecond hyperfield:")
        p2 = int(input("  prime p2: "))
        k2 = int(input("  degree k2: "))
        d2 = int(input(f"  d2 | {p2**k2 - 1}: "))
    except (ValueError, KeyboardInterrupt):
        print("\nCancelled.")
        return

    h1 = QuotientHyperfield.from_params(p1, k1, d1)
    h2 = QuotientHyperfield.from_params(p2, k2, d2)
    result = are_isomorphic(h1, h2, method="auto")
    print(f"\n{h1}")
    print(f"{h2}")
    print(f"\nResult: {result}")
    if result.isomorphic:
        print("The hyperfields are isomorphic.")
    else:
        print("The hyperfields are not isomorphic.")


def run_analysis() -> None:
    print("Hyperfield Analysis and Visualization")
    print("=" * 40)
    try:
        n = int(input("Hyperfield order n: "))
        max_p = int(input("max_p: "))
        max_k = int(input("max_k: "))
    except (ValueError, KeyboardInterrupt):
        print("\nCancelled.")
        return

    triples = generate_triples(n, max_p, max_k)
    print(f"Found {len(triples)} triples.")
    if not triples:
        return

    classes = classify_hyperfields(n, max_p, max_k)
    print(f"Found {len(classes)} isomorphism classes.")
    for i, c in enumerate(classes):
        print(f"  class {i}: {c}")

    chars, cchars = analyze_characteristics(classes)
    print("Characteristics:", chars)
    print("C-characteristics:", cchars)
    visualize_characteristics(chars, cchars)


def main() -> None:
    while True:
        print("\nSelect an option:")
        print("1. Generate additive table for a specific hyperfield")
        print("2. Check isomorphism between two hyperfields")
        print("3. Run hyperfield analysis and visualization")
        print("4. Exit")
        choice = input("Enter your choice (1-4): ").strip()
        if choice == "1":
            generate_and_display_table()
        elif choice == "2":
            run_isomorphism_check()
        elif choice == "3":
            run_analysis()
        elif choice == "4":
            print("Exiting...")
            break
        else:
            print("Invalid choice.")


if __name__ == "__main__":
    main()
