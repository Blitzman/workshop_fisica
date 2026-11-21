# N3_modes_fullgr_axion: modes with the axion and the dynamic metric (beyond Cowling)

## Run
```bash
python -m N3_modes_fullgr_axion.run                  # steps 10 and 25: 12 values of nu_a + convergence study (about 70 + 26 mode searches)
python -m N3_modes_fullgr_axion.run --quick          # only the 6 kHz star, about 1 minute
```
`--procs N` (default min(6, cores)): fork pool over stars (step 10) and tasks (step 25).

**Compiled equations (cache).** The numerical systems compiled from the equation artifacts are cached in `compiled/extended.json` (tracked by git; `common/compiled.py`) and loaded in milliseconds. If an equation artifact or the reduction code changes, the run stops with a StaleCache error; `--recompile` rebuilds the cache (N2: about 1.5 min, N3: 16 s, N4: 22 s).

## Contract
Reads: D2 `fullgr_extended.eq.json`; S1 `eos_bsk26.csv`, `star_noaxion.csv`, `star_nua<nu>.csv`, `manifest.json`; R `paper2_tables.ref.json`; N1 `modes_cowling.csv` (optional).

| artifact | kind | content |
|---|---|---|
| `modes_axion.csv` + `.meta.json` | table | for every nu_a and mode (f, p1, p2): background-only and coupled Re, Im (kHz), and the no-axion references |
| `exterior_convergence.json` | data | step 25: branch tests, the 6 delicate modes under 4 variants, negative controls |
| `manifest.json` | checks | the checks of steps 10 and 25 |

## Method
The mode is a complex zero of the INGOING gravitational amplitude (Muller). The numerical matrices come from D2's extended system (symbolic compile, about 25 s) as
y' = A y with A = A_-1/w^2 + A_0 + w^2 A_1; the three real coefficient matrices (times r^2, smooth at the centre) are evaluated EXACTLY at the nodes of S1's table and spline-interpolated
(`tables.py`), so a right-hand-side call costs about 10 us instead of about 520 us. Phi'' is analytic (not a finite difference). The matching radius is `R + 0.7 (r_end - R)`.

## Measured against the legacy (nu_a = 6 kHz, f-mode, same theta_c)
tables vs legacy exact path on the same star: 4.0e-9 (Re), 3.5e-8 (Im); against the legacy on the r-integrated star: 1.2e-8, 3.1e-6 (star definition: 9 m truncation, 10 m start);
f-mode root 4.3 s against 195 s (45x); convergence of the root: rtol 1e-10 -> 1e-12 moves Re 5.7e-10 and Im 4.9e-8; matching radius 50 -> 60 km moves 1.4e-8 and 3.4e-6.
The 12-star comparison is `tests/n3_gate.py` (see MIGRATION.md for the result).

Not in this module: the LMI polytrope end-to-end check (needs a polytrope EOS table in S1; N2 validates the same equations) and the comparison with N2's matrix
(`tests/compare_n3_n2_matrix.py`).
