"""
Visualization functions for hyperfield analysis.

This module contains functions for visualizing the distribution of characteristics
and c-characteristics of hyperfields.
"""

import matplotlib.pyplot as plt
import numpy as np

def visualize_characteristics(characteristics, c_characteristics):
    """
    Visualizes the distribution of characteristics and c-characteristics using histograms.
    Adjusts x-axis ticks to be integers for clarity and y-axis ticks to be integers for frequency.
    Filters out None values before plotting.
    """
    valid_characteristics = [c for c in characteristics if c is not None]
    valid_c_characteristics = [cc for cc in c_characteristics if cc is not None]

    if not valid_characteristics and not valid_c_characteristics:
        print("No valid characteristic data to visualize.")
        return

    plt.figure(figsize=(12, 5))

    # Histogram for Characteristics
    plt.subplot(1, 2, 1)
    if valid_characteristics:
        bins_char = np.arange(min(valid_characteristics)-0.5, max(valid_characteristics)+1.5, 1)
        plt.hist(valid_characteristics, bins=bins_char, edgecolor='black')
        plt.xticks(np.arange(min(valid_characteristics), max(valid_characteristics)+1, 1))
        plt.yticks(np.arange(0, max(plt.yticks()[0]) + 1, 1))
    plt.title("Distribution of Characteristics")
    plt.xlabel("Characteristic")
    plt.ylabel("Frequency")

    # Histogram for C-characteristics
    plt.subplot(1, 2, 2)
    if valid_c_characteristics:
        bins_cchar = np.arange(min(valid_c_characteristics)-0.5, max(valid_c_characteristics)+1.5, 1)
        plt.hist(valid_c_characteristics, bins=bins_cchar, edgecolor='black')
        plt.xticks(np.arange(min(valid_c_characteristics), max(valid_c_characteristics)+1, 1))
        plt.yticks(np.arange(0, max(plt.yticks()[0]) + 1, 1))
    plt.title("Distribution of C-characteristics")
    plt.xlabel("C-characteristic")
    plt.ylabel("Frequency")

    plt.tight_layout()
    plt.show()
