# Agent Protocol for the modular axion workflow

Read the module's `CONTRACT.md`, the modular `README.md`, and the shared `CONVENTIONS.md` and `LAB_NOTEBOOK.md` before changing code or running calculations.

## Scientific verification

- Derive equations from the action in SymPy and verify agreement with available references by reducing differences to exact zero modulo the background equations. Keep derivations independent of published transcriptions in `R_paper_refs/`. For new symbolic checks, verify that the comparator rejects at least two deliberate mutations.
- Compile numerical equations from verified symbolic artifacts. Treat the equation PDFs as generated documentation. Document any existing exceptions in `MIGRATION.md`.
- Establish an independent published benchmark in the no-axion limit before treating a new axion result as validated. Use quantities derived from the equations or cited published values for physical coefficients. Record scientific hypotheses and outcomes in the lab notebook.

## Modules and runs

- Keep modules independently runnable. Exchange data through their contracted artifacts and use `common/` for shared infrastructure. Generate report numbers from module outputs.
- Start with the relevant module's witness or quick target where available. Use the full sweep after the focused checks pass. Keep symbolic diagnostics to residuals and affected coefficients. Use exact polynomial reduction for spherical harmonics.
- Use `MIGRATION.md` to locate the original calculation and known differences. The archived workflow and its instructions are under `archive/axion_ns/`.
