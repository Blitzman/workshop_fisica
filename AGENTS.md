# Repository Agent Protocol

This repository contains two independent workstreams with different rules depending on the subdirectory:

- **If working in `axion_ns/` (or `report/`, `results/`):** read `axion_ns/AGENTS.md`, `CONVENTIONS.md`, and `LAB_NOTEBOOK.md` before writing code or running calculations.
- **If working in `ns_cooling/` (or `observations/`):** read `ns_cooling/AGENTS.md` and `ns_cooling/README.txt`. The `axion_ns/`-specific rules (SymPy from the action, zero free constants, `CONVENTIONS.md`) do **not** apply to `ns_cooling/` (decision D-026).

## Common Rules Across the Repository
- **Hypothesis and outcome:** before testing a non-trivial idea, record the hypothesis in the corresponding notebook/tracker (`LAB_NOTEBOOK.md` for axions, or the `ns_cooling` tracker); afterwards record whether it was confirmed or refuted and why.
- **External validation:** no new result counts as validated unless it first reproduces a published result by other authors.
- **Physical premise reviewed by the group:** before writing code for a new line of work, record the premise and review it with the group/José. Frame results as preliminary ("we obtain", "suggests"), not as discoveries.
- **HTML/LaTeX token economy:** when reading or editing long HTML reports (such as `docs/cooling_comparison.html` or `report/*.html`), use line ranges (`StartLine`/`EndLine`) without re-reading entire files.
