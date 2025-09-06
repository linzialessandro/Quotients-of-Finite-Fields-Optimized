"""
Main entry point for the hyperfield analysis tool.

This module provides a command-line interface for generating and analyzing
quotient hyperfields.
"""

from hyperfield.analysis import (
    generate_triples,
    classify_hyperfields,
    analyze_characteristics,
)
from hyperfield.isomorphism import are_isomorphic_optimized
from visualization import visualize_characteristics
from utils import is_prime

def run_analysis():
    """
    Orchestrates the process of generating triples, classifying hyperfields,
    analyzing characteristics, and visualizing the results based on user input.
    """
    print("Hyperfield Analysis and Visualization (Optimized Isomorphism Check)")
    print("=" * 55)

    while True:
        try:
            n = int(input("Enter the desired number of elements in the hyperfield (n): "))
            if n <= 1:
                print("n must be greater than 1.")
                continue
            max_p = int(input("Enter the maximum value for the prime p (max_p): "))
            if max_p < 2:
                print("max_p must be at least 2.")
                continue
            max_k = int(input("Enter the maximum value for the extension degree k (max_k): "))
            if max_k < 1:
                print("max_k must be at least 1.")
                continue
            break
        except ValueError:
            print("Invalid input. Please enter integer values.")

    print(f"\nGenerating triples for n={n}, max_p={max_p}, max_k={max_k}...")
    triples = generate_triples(n, max_p, max_k)
    print(f"Found {len(triples)} valid triples.")

    if not triples:
        print("No valid triples found for the given inputs. Cannot proceed with classification or analysis.")
        return

    print("\nClassifying hyperfields into isomorphism classes (using optimized check)...")
    isomorphism_classes = classify_hyperfields(n, max_p, max_k)
    print(f"Found {len(isomorphism_classes)} distinct isomorphism classes.")

    if not isomorphism_classes:
        print("No isomorphism classes found. Cannot proceed with analysis or visualization.")
        return

    print("\nAnalyzing characteristics for each isomorphism class representative...")
    characteristics, c_characteristics = analyze_characteristics(isomorphism_classes)
    print("Characteristic values:", characteristics)
    print("C-characteristic values:", c_characteristics)

    print("\nVisualizing characteristic distributions...")
    visualize_characteristics(characteristics, c_characteristics)

def main_isomorphism_check():
    """Interactive main function for testing the generalized isomorphism check."""
    print("Generalized Isomorphism Check for Quotient Hyperfields GF(p^k)/G_d")
    print("=" * 70)

    while True:
        try:
            print("\nFirst hyperfield GF(p1^k1)/G_d1:")
            p1 = int(input("Enter prime p1: "))
            k1 = int(input("Enter extension degree k1: "))
            q1 = p1 ** k1
            print(f"p1^k1 - 1 = {q1 - 1}")
            d1 = int(input(f"Enter divisor d1 of {q1 - 1}: "))

            print("\nSecond hyperfield GF(p2^k2)/G_d2:")
            p2 = int(input("Enter prime p2: "))
            k2 = int(input("Enter extension degree k2: "))
            q2 = p2 ** k2
            print(f"p2^k2 - 1 = {q2 - 1}")
            d2 = int(input(f"Enter divisor d2 of {q2 - 1}: "))

            if not is_prime(p1) or not is_prime(p2):
                print("Error: Both p1 and p2 must be prime numbers.")
                continue
            if k1 <= 0 or k2 <= 0:
                print("Error: Both k1 and k2 must be positive integers.")
                continue
            if d1 <= 0 or (q1 - 1) % d1 != 0:
                print(f"Error: d1 must divide {q1 - 1}.")
                continue
            if d2 <= 0 or (q2 - 1) % d2 != 0:
                print(f"Error: d2 must divide {q2 - 1}.")
                continue

            break
        except ValueError:
            print("Invalid input, please enter integer values.")
        except KeyboardInterrupt:
            print("\nExiting...")
            return

    print(f"\nChecking isomorphism between GF({p1}^{k1})/G_{d1} and GF({p2}^{k2})/G_{d2}...")
    result = are_isomorphic_optimized(p1**k1, d1, p2**k2, d2)
    print(f"\nResult: {'The generated hyperfields are isomorphic.' if result else 'The generated hyperfields are not isomorphic.'}")

if __name__ == "__main__":
    while True:
        print("\nSelect an option:")
        print("1. Run hyperfield analysis and visualization")
        print("2. Check isomorphism between two hyperfields")
        print("3. Exit")

        choice = input("Enter your choice (1, 2, or 3): ")

        if choice == '1':
            run_analysis()
        elif choice == '2':
            main_isomorphism_check()
        elif choice == '3':
            print("Exiting...")
            break
        else:
            print("Invalid choice. Please try again.")