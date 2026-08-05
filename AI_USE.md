# AI use disclosure

This document records how generative AI was used in this project, for transparency toward readers, referees, and future journal policies.

## Tools

| Tool | Provider | Role in this project |
|------|----------|----------------------|
| **Grok Build CLI** | xAI | Interactive coding agent (file edits, tests, packaging, git, docs) |
| **Grok 4.5** | xAI | Language model backend used by Grok Build CLI |

## Scope of use

Assistance covered **software engineering and research-engineering** tasks, including (as applicable over the course of development):

- scaffolding and refining the Python package layout;
- implementing and debugging library modules and experiment drivers;
- writing and updating tests;
- packaging (`pyproject.toml`), CLI entry points, and reproducibility notes;
- documentation (README, citation metadata) and repository housekeeping;
- linking the codebase to the companion preprint (arXiv metadata, release tags).

## Out of scope / non-claims

- Generative AI systems are **not** co-authors of the paper or software.
- They are **not** claimed as independent sources of mathematical theorems or proofs.
- Any AI-suggested text, code, or bibliography was subject to human review before inclusion.

## Human responsibility

**Alessandro Linzi** designed the research program, owns all scientific claims, verified proofs and computational outputs against the intended mathematics and the literature, and takes full responsibility for the published paper, the software release (e.g. v0.1.0), and any errors.

## Where this is disclosed

| Venue | Location |
|-------|----------|
| Companion paper | Acknowledgements (see `paper/main.tex`; arXiv:2608.03625) |
| This repository | [README § AI use disclosure](README.md#ai-use-disclosure) and this file |

## Template for journal submission

If a journal asks for a dedicated “Use of AI” or “Declaration of generative AI” statement, the following short form is usually sufficient (adapt if the journal template differs):

> Development of the companion software and related repository materials was assisted by Grok Build CLI with the language model Grok 4.5 (xAI), used for programming tasks such as implementation, tests, packaging, documentation, and repository maintenance. No generative AI system is an author of this work. All mathematical statements, proofs, experimental designs, reported tables, and final text were reviewed and verified by the author, who takes full responsibility for the content and for the correctness of computational results.

## Updates

Update this file if the tools, model versions, or scope of AI assistance change in a material way (e.g. AI-assisted drafting of paper sections, figure generation, or automated theorem discovery).
