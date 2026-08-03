"""Visualization helpers for census experiments."""

from __future__ import annotations

import math

import matplotlib.pyplot as plt
import numpy as np


def visualize_characteristics(characteristics, c_characteristics) -> None:
    """Histograms of char / C-char (finite values only)."""

    def finite_ints(values):
        out = []
        for c in values:
            if c is None:
                continue
            if isinstance(c, float) and math.isinf(c):
                continue
            out.append(int(c))
        return out

    valid_char = finite_ints(characteristics)
    valid_cchar = finite_ints(c_characteristics)

    if not valid_char and not valid_cchar:
        print("No valid characteristic data to visualize.")
        return

    plt.figure(figsize=(12, 5))

    plt.subplot(1, 2, 1)
    if valid_char:
        bins = np.arange(min(valid_char) - 0.5, max(valid_char) + 1.5, 1)
        plt.hist(valid_char, bins=bins, edgecolor="black")
        plt.xticks(np.arange(min(valid_char), max(valid_char) + 1, 1))
    plt.title("Distribution of Characteristics")
    plt.xlabel("Characteristic")
    plt.ylabel("Frequency")

    plt.subplot(1, 2, 2)
    if valid_cchar:
        bins = np.arange(min(valid_cchar) - 0.5, max(valid_cchar) + 1.5, 1)
        plt.hist(valid_cchar, bins=bins, edgecolor="black")
        plt.xticks(np.arange(min(valid_cchar), max(valid_cchar) + 1, 1))
    plt.title("Distribution of C-characteristics")
    plt.xlabel("C-characteristic")
    plt.ylabel("Frequency")

    plt.tight_layout()
    plt.show()
