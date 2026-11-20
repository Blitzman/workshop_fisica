# Repository Agent Protocol

This repository contains two independent workstreams: axion physics and neutron-star cooling.

- **`axion_ns/` is the primary implementation for axion work.** Before writing code or running calculations, read `axion_ns/AGENTS.md`, `axion_ns/README.md`, the relevant module's `CONTRACT.md` where available, `CONVENTIONS.md`, and the relevant entries in `LAB_NOTEBOOK.md`. Use `axion_ns/MIGRATION.md` to map legacy steps to modules.
- **`archive/axion_ns/`, root `report/`, and root `results/` are the legacy axion workflow.** Use them for historical comparisons, reproducing legacy results, or work explicitly targeting the legacy implementation. Read `archive/axion_ns/AGENTS.md`, `CONVENTIONS.md`, and `LAB_NOTEBOOK.md` before writing code or running calculations there.
- **If working in `ns_cooling/` (or `observations/`):** read `ns_cooling/AGENTS.md` and `ns_cooling/README.txt`. The `axion_ns/`-specific rules (SymPy from the action, zero free constants, `CONVENTIONS.md`) do **not** apply to `ns_cooling/` (decision D-026).

## Common Rules Across the Repository
- **Hypothesis and outcome:** before testing a non-trivial idea, record the hypothesis in the corresponding notebook/tracker (`LAB_NOTEBOOK.md` for axions, or the `ns_cooling` tracker); afterwards record whether it was confirmed or refuted and why.
- **External validation:** no new result counts as validated unless it first reproduces a published result by other authors.
- **Physical premise reviewed by the group:** before writing code for a new line of work, record the premise and review it with the group/José. Frame results as preliminary ("we obtain", "suggests"), not as discoveries.
- **HTML/LaTeX token economy:** when reading or editing long HTML reports (such as `docs/cooling_comparison.html` or `report/*.html`), use line ranges (`StartLine`/`EndLine`) without re-reading entire files.
