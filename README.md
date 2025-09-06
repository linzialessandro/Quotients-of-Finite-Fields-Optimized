# Hyperfield Additive Table Generator and Isomorphism Analyzer

Welcome to this repository! This project provides tools to explore quotient hyperfields of the form GF(p^k) / G_d, where G_d is a multiplicative subgroup of GF(p^k) of order `d`.

## Features

*   **Generate Additive Tables:** You can input a prime `p`, an integer `k`, and a divisor `d` of `p^k - 1` to see the additive hyperoperation table for the resulting hyperfield.
*   **Check for Isomorphism:** You can input parameters for two hyperfields and the tool can tell you if they are isomorphic.
*   **Analyze and Visualize:** You can specify parameters to generate a range of potential hyperfields, classify them into isomorphism classes, and visualize the distribution of their characteristics and c-characteristics.

## Requirements

To run this project, you'll need the following Python libraries:

*   Python 3.6+
*   pandas
*   galois
*   matplotlib
*   numpy

You can install the dependencies using pip:

```bash
pip install -r requirements.txt
```

## How it Works

Here's a brief overview of the process implemented:

1.  The `galois` library is used to construct finite fields GF(p^k).
2.  A primitive element is identified to construct the multiplicative subgroup of order `d`.
3.  The coset representatives of GF(p^k)^*/G_d are determined.
4.  The additive hyperoperation is computed based on Krasner's quotient construction for pairs of cosets (including the zero element).

For isomorphism checking, the characteristics and c-characteristics of the hyperfields are compared.

## Project Structure

The project is organized into several modules for better readability and maintainability:

*   `main.py`: The main entry point of the application, providing a command-line interface.
*   `utils.py`: Contains general utility functions like primality testing and prime factorization.
*   `visualization.py`: Handles the visualization of hyperfield characteristics.
*   `hyperfield_modules/`:
    *   `core.py`: Implements the core logic for hyperfield construction and operations.
    *   `isomorphism.py`: Contains functions for checking isomorphism between hyperfields.
    *   `analysis.py`: Provides functions for generating and classifying hyperfields for analysis.

## Usage

To run the application, execute `python main.py`:

```bash
python main.py
```

You will be presented with a menu to choose from:

1.  **Generate additive table for a specific hyperfield:** This option prompts for `p`, `k`, and `d` and displays the additive hyperoperation table for the resulting hyperfield.
2.  **Check isomorphism between two hyperfields:** This option lets you input parameters for two hyperfields and determine if they are isomorphic.
3.  **Run hyperfield analysis and visualization:** This option allows you to generate, classify, analyze, and visualize hyperfields based on user-defined parameters (n, max_p, max_k).
4.  **Exit:** Exits the application.
