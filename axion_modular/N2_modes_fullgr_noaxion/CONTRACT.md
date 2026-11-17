# N2_modes_fullgr_noaxion: validation of the extended system against the literature

## Run
```bash
python -m N2_modes_fullgr_noaxion.run              # full validation, from axion_modular/ (about 1,000 s of CPU / --procs)
python -m N2_modes_fullgr_noaxion.run --quick --out /tmp/n2_quick     # smoke run, LMI star only, about 2 minutes
```
`--procs N` (default min(6, cores)) uses a fork pool (the compiled systems cannot be pickled); without `fork` it runs serially.
`--cowling-rtol 1e-8` relaxes only the Cowling root search (the amplitude tolerances are fixed, see below).

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
moves the f-mode by more than 10%.

## Numerics
The systems are compiled from the artifacts with `lambdify(cse=True)` (A y compiled directly): 2.7x faster than the legacy per mode with
matrices equal to 1e-15. rtol 1e-11 and atol 1e-14 must not be loosened for the amplitude: |A_out|^2 is noisy at 1e-3..1e-2 (the two interior
solutions grow to W ~ 1e6 at the surface and the X(R) = 0 combination cancels that by ~3e5), rtol 1e-10 already moves the LMI f-mode by 7.7e-5
(> 5e-5) and atol does not help (MIGRATION.md, item 11).

Warnings: the integration may hit a divide-by-zero at the last point (r/R = 1 - 1e-7, the polytrope surface); the solver rejects that step and
retries, exactly as the legacy code does. N2 silences that floating-point noise inside the integration only; solver failures still raise.

Not in this module: the axion background on a polytrope (`solve_bvp`) and the figures of the legacy step 8.
