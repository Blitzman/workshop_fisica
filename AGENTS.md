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

## Verification and Run Policy (all workstreams)
Run the smallest check that answers the question. Nothing is rerun "to be safe".

1. **Ask before anything expected to take more than 2 minutes** (wall time on the machine at hand): a test, a module run, a gate script or a sweep. Say what it is, why it is needed and the expected time. Estimate from the runtimes below or the module docs; if unknown, ask.
2. **Pushes, merges, copies and file transfers:** verify with git only (`git diff --stat A B`, `git status`, hashes). If the trees are identical, nothing is run.
3. **Documentation, notebook or comment changes:** no tests.
4. **Code change in one module:** run that module's own tests (`pytest tests/test_<module>.py` or the matching `make test-*`) and its witness or quick target. Not the whole suite.
5. **A module's contracted outputs change** (tables, JSON artifacts): run the quick targets of the modules that read them. Full consumer runs only when the quick runs show a difference that matters, and under rule 1.
6. **Full sweeps (`make all`, `run_all`), slow tests and on-request comparison scripts:** only on request or before a release, under rule 1.
7. **Do not repeat a check** whose inputs have not changed since it last passed in the session. Report the earlier result instead.

Known wall times (axion_ns, 1–2 cores): `pytest tests` 10–20 s; S1 25 s; S1-witness, N1-witness, N2/N3/N4-quick under 1 min; T about 2 min; S1-extra about 2 min; N1 about 4 min; N3 about 4 min; D2 about 8 min; N2 about 17 min of CPU; N4 about 27 min; `run_all` (legacy) about 75 min. Anything at or above 2 minutes needs rule 1.
