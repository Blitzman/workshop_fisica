# N4_modes_radial: radial (l = 0) modes of the star with the axion condensate (full GR)

**Stage: development** (Paper I Tables 1-2 not yet within the pre-set criteria). Set by José; see the root `AGENTS.md`, Module Stages.

## Run
```bash
make S1-extra                                   # once: the extra stars of S1 (about 2 minutes, not tracked by git)
python -m N4_modes_radial.run                   # steps 22, 26 and 15 (about 10 minutes on one core)
python -m N4_modes_radial.run --part all        # plus 19, 20, 16 and the analysis parts (the contour counts: long; --procs divides it)
python -m N4_modes_radial.run --part analysis   # steps 18, 21, 23, 24 (about 100 s); they read the outputs of steps 15 and 20
```
**Compiled equations (cache).** The numerical systems compiled from the equation artifacts are cached in `compiled/radial.json` (tracked by git; `common/compiled.py`) and loaded in milliseconds. If an equation artifact or the reduction code changes, the run stops with a StaleCache error; `--recompile` rebuilds the cache (N2: about 1.5 min, N3: 16 s, N4: 22 s).

## Contract
Reads: D0 `background.eq.json`; D1 `radial_linear.eq.json`; S1 `eos_bsk26.csv`, `star_*.csv`, `manifest.json` and the untracked `outputs_extra/`; R `paper1_tables.ref.json`.

| artifact | step | content |
|---|---|---|
| `radial_contour.json` | 22 | the exterior condition against the complex ray and the starting point |
| `radial_paper1.json` | 26 | Tables 1-3 of Paper I against ours (28 + 28 + 14 modes) |
| `radial_star.json` | 15 | turning point, F, the gravitational-atom ladder, H1, the retired artefact, stability on the imaginary axis, robustness |
| `radial_crossing.json` | 19 | avoided crossing F - gravitational atom n = 1 |
| `radial_damped.json` | 20 | count in the damped half-plane (and the old-lapse positive control) |
| `radial_stability.json` | 16 | argument-principle count in the unstable half-plane |
| `radial_rydberg.json` | 18 | Rydberg law of the gravitational-atom ladder (reads step 15) |
| `radial_deep.json` | 21 | deep zeros of step 20 against the matching radius (reads step 20) |
| `radial_eos.json` | 23 | sensitivity to the thermodynamic inconsistency of the BSk26 fits (three c_s^2 definitions; reads step 15) |
| `radial_scattering.json` | 24 | width from the phase of the real-frequency scattering (reads step 15) |
| `manifest.json` | | the checks |

Same keys as the legacy `results/radial_*.json`. Step 17 (Font et al. 2002) is not part of this module by decision; step 15 and 19's star construction on the fly is replaced by S1's extra stars.

## Method
The radial system y = (xi, DP, dtheta, dtheta') is derived from D0's and D1's equation artifacts (`systems.py`, a port of the legacy `radial.py`) and compiled as flat lists. For this system
A(r, w) = A_0(r) + w^2 A_1(r) EXACTLY (residual 1e-17): `tables.py` evaluates both coefficient matrices exactly at the stored nodes of S1's table and interpolates r^2 A, so one right-hand side
costs about 7 us instead of about 200 us and one `det_entire` about 0.11 s instead of about 2.7 s. `modes.py` ports `interior`, `outward`, `det`, `det_entire`, `exterior_psi`, `fluid_fraction`
(the winding and Muller functions are the legacy ones verbatim; `exterior_L` is N1's). The matching radius is `R + 0.7 (r_end - R)`.
