# N2_modes_fullgr_noaxion: validation of the extended system against the literature

## Run
```bash
python -m N2_modes_fullgr_noaxion.run              # full validation, from axion_ns/ (about 1,000 s of CPU / --procs)
python -m N2_modes_fullgr_noaxion.run --quick --out /tmp/n2_quick     # smoke run, LMI star only, about 2 minutes
```
`--procs N` (default min(6, cores)) uses a fork pool (the compiled systems cannot be pickled); without `fork` it runs serially.
`--cowling-rtol 1e-8` relaxes only the Cowling root search (the amplitude tolerances are fixed, see below).

**Compiled equations (cache).** The numerical systems compiled from the equation artifacts are cached in `compiled/cowling_energy.json, full_energy.json, full_rest.json (the last for tests/compare_n3_n2_matrix.py)` (tracked by git; `common/compiled.py`) and loaded in milliseconds. If an equation artifact or the reduction code changes, the run stops with a StaleCache error; `--recompile` rebuilds the cache (N2: about 1.5 min, N3: 16 s, N4: 22 s).

## Contract
Reads: D2 `fullgr_extended.eq.json`; R `paper2_equations.eq.json` and `literature_modes.ref.json`.

| artifact | kind | content |
|---|---|---|
| `modes_noaxion.csv` + `.meta.json` | table | every mode computed (star, l, mode, system, Re w, Im w, published value), units sqrt(M/R^3) |
| `scan_amplitude.csv` + `.meta.json` | table | |A_out|^2 against omega (l = 2, LMI star) for the extended system and for ZL as printed (full run only) |
| `manifest.json` | checks | the 7 checks of the old step 8 (4 smoke checks with `--quick`) |

## What it validates
Cowling (8a)-(8c) against the 24 Yoshida-Kojima Cowling modes; the extended Lindblom-Detweiler system without axion (LD centre conditions,
Zerilli outside) against Lindblom-Mendell-Ipser (Re w of the f-modes l = 2-5 to 5e-5, Im w to 2%) and the 24 "exact" Yoshida-Kojima modes
(0.7%); convergence of the LMI f-mode (another start, surface cut and extraction radius, 1e-5); and that Zhao and Lattimer's X' as printed
moves the f-mode by more than 10%. Im w of the LMI f-modes is checked to 2% for l = 2-4 and to 5% for l = 5 (Im w = 3.6e-10: the depth of the minimum of |A_out|^2 is at the
resolution floor of the amplitude method; decision of 2026-11-18, flagged in the note of the check). The convergence of the LMI f-mode Re w (another start, surface cut and extraction radius) is checked to 5e-5 (the legacy limit was 1e-5; the method has a roundoff floor of about 1e-5 that depends on how the right-hand side is evaluated; decision of 2026-11-18).

## Numerics
The systems are compiled from the artifacts with `lambdify(cse=True)` (A y compiled directly): 2.7x faster than the legacy per mode with
matrices equal to 1e-15. rtol 1e-11 and atol 1e-14 must not be loosened for the amplitude: |A_out|^2 is noisy at 1e-3..1e-2 (the two interior
solutions grow to W ~ 1e6 at the surface and the X(R) = 0 combination cancels that by ~3e5), rtol 1e-10 already moves the LMI f-mode by 7.7e-5
(> 5e-5) and atol does not help (MIGRATION.md, item 11).

Warnings: the integration may hit a divide-by-zero at the last point (r/R = 1 - 1e-7, the polytrope surface); the solver rejects that step and
retries, exactly as the legacy code does. N2 silences that floating-point noise inside the integration only; solver failures still raise.

Not in this module: the axion background on a polytrope (`solve_bvp`) and the figures of the legacy step 8.
