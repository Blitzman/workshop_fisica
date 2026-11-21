# Agent Protocol for the modular axion workflow

Read the module's `CONTRACT.md`, the modular `README.md`, and the shared `CONVENTIONS.md` and `LAB_NOTEBOOK.md` before changing code or running calculations.

## Scientific verification

- Derive equations from the action in SymPy and verify agreement with available references by reducing differences to exact zero modulo the background equations. Keep derivations independent of published transcriptions in `R_paper_refs/`. For new symbolic checks, verify that the comparator rejects at least two deliberate mutations.
- Compile numerical equations from verified symbolic artifacts. The compiled functions of N2, N3 and N4 are cached in `<module>/compiled/` (the equations are settled); rebuild them only on request, with the module's `--recompile` flag. Treat the equation PDFs as generated documentation. Document any existing exceptions in `MIGRATION.md`.
- Establish an independent published benchmark in the no-axion limit before treating a new axion result as validated. Use quantities derived from the equations or cited published values for physical coefficients. Record scientific hypotheses and outcomes in the lab notebook.

## Modules and runs

- Keep modules independently runnable. Exchange data through their contracted artifacts and use `common/` for shared infrastructure. Generate report numbers from module outputs.
- Tests and runs follow the **Verification and Run Policy** in the root `AGENTS.md` (smallest check first; ask before anything over 2 minutes). Witness and quick targets: `make S1-witness`, `N1-witness`, `N2-quick`, `N3-quick`, `N4-quick`, `T-fast`. Keep symbolic diagnostics to residuals and affected coefficients. Use exact polynomial reduction for spherical harmonics.
- Use `MIGRATION.md` to locate the original calculation and known differences. The archived workflow and its instructions are under `archive/axion_ns/`.
