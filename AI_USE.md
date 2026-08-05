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
