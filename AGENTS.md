# Repository Agent Protocol

The active workstream is the axion physics in `axion_ns/`. Before writing code or running calculations, read `axion_ns/AGENTS.md`, `axion_ns/README.md`, the relevant module's `CONTRACT.md`, `CONVENTIONS.md`, and the relevant entries in `LAB_NOTEBOOK.md`.

- **The current `axion_ns/` is the reference.** `archive/` (including `archive/axion_ns/` and the migration scripts), root `report/` and root `results/` are the legacy workflow, kept for history only. Do not compare with them or run them unless José asks; `axion_ns/MIGRATION.md` records what was compared during the migration.
- **`ns_cooling/` and `observations/` are obsolete** (to be replaced by a new project). Do not work on them, test them or update them unless José asks.

## Common Rules
- **Hypothesis and outcome:** before testing a non-trivial idea, record the hypothesis in `LAB_NOTEBOOK.md`; afterwards record whether it was confirmed or refuted and why.
- **External validation:** no new result counts as validated unless it first reproduces a published result by other authors.
- **Physical premise reviewed by the group:** before writing code for a new line of work, record the premise and review it with the group/José. Frame results as preliminary ("we obtain", "suggests"), not as discoveries.
- **New numerical methods are discussed first:** before trying a different numerical method (integrator, variables, shooting or matching strategy, root finder, discretisation, tolerances used as a method change), explain to José what will be tried and for what purpose, and wait for his agreement before running anything.
- **HTML/LaTeX token economy:** when reading or editing long HTML reports (such as `report/*.html`), use line ranges without re-reading entire files.

## Module Stages
Each module states its stage in the first lines of its `CONTRACT.md`. Only José changes a stage; record the change in `LAB_NOTEBOOK.md`.

1. **Initial.** A new line of work: no reference yet, code and equations still changing.
   - Checks: only that it runs and gives sensible numbers, plus a known limit when it is cheap. No test files, no stored reference values, no comparisons with earlier runs.
   - Nothing is claimed as a result. Moves to *development* when it reproduces one published result by other authors.
2. **Development.** A reference exists; the work improves, extends or speeds up the module.
   - Checks: the published benchmark it reproduces, plus a few fast tests of the pieces that have stopped changing.
   - Comparisons with the module's previous outputs only when the numerics are deliberately changed; record them once in the notebook, not as permanent tests. Moves to *production* when José decides the results go into a paper or report.
3. **Production.** Results used in papers or reports; the committed outputs are the reference.
   - Equations are frozen; compiled equations are cached and rebuilt only on request (`--recompile`).
   - A change to code or outputs needs a documented comparison with the committed outputs before it is merged.

## Verification and Run Policy
Run the smallest check that answers the question. Nothing is rerun "to be safe".

1. **Ask before anything expected to take more than 2 minutes** (wall time on the machine at hand): a test, a module run or a sweep. Say what it is, why it is needed and the expected time. Estimate from the times below; if a time is not listed, ask.
2. **Pushes, merges, copies and file transfers:** verify with git only (`git diff --stat A B`, `git status`, hashes). If the trees are identical, nothing is run.
3. **Documentation, notebook or comment changes:** no tests.
4. **Code change in one module:** run that module's own tests (`pytest tests/test_<module>.py`) and its witness or quick target, as its stage allows. Not the whole suite.
5. **A module's contracted outputs change** (tables, JSON artifacts): run the quick targets of the modules that read them. Full consumer runs only when the quick runs show a difference that matters, and under rule 1.
6. **Full runs (`make all`, a module's full run) and on-request comparison scripts:** only on request or before a release, under rule 1.
7. **Do not repeat a check** whose inputs have not changed since it last passed in the session. Report the earlier result instead.

Known wall times (axion_ns, 2 cores, compiled equations cached). Measured 2026-11-20/21: `pytest tests` 13 s; S1 24 s; N3 `--quick` 21 s; N1 about 4.5 min; N3 about 4 min; N4 `--part all` about 27 min. From the module docs: R 1.5 s; D0 4 s; D1 24 s; S1 `--extra` about 2 min; T about 2 min; N2 `--quick` about 2 min; D2 about 8 min; N2 about 17 min of CPU. Not measured (ask): N1 `--nu-a 6`, N4 `--part fast`. Anything at or above 2 minutes needs rule 1.
