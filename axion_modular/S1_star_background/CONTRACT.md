# S1_star_background: BSk26 star with axion condensate: profile tables

## Run
```bash
python -m S1_star_background.run            # from axion_modular/ ; `--help` lists the options
```
Outputs go to `S1_star_background/outputs/` (override with `--out`). Upstream artifacts are read from the sibling `outputs/` folders
(override with the `--*-dir` options), so this module runs on its own as soon as those files exist.

## Contract
Reads: nothing (BSk26 analytic fits in `bsk26.py`).

| artifact | kind | content |
|---|---|---|
| `eos_bsk26.csv + .meta.json` | table | n [fm^-3], eps and p [MeV fm^-3] (monotone-in-p subset used by the solver) |
| `star_noaxion.csv, star_nua<nu>.csv + .meta.json` | table | dense STATE table: r, ln p, m, Phi, theta, theta' on 2,499 nodes (2,000 inside, surface-clustered; 500 in the tail); geometric units, km; metadata: M, R, theta_c, Phi_shift, r_end, r_k, constants, integrator settings, interpolation error |
| `manifest.json `extra`` | data | `stars` (M, R, theta_c, rows, error, seconds), `condensate_exists`, `nu_a_list` |

Every run also writes `outputs/manifest.json` (checks, runtime, environment, sha256 of each artifact and of each input).

Consumers interpolate the five columns with a cubic spline (not-a-knot), as TWO splines split at the node `r_k` (the radius where the EOS is switched off; one spline across it rings). S1 measures exactly that interpolation against the integrator's dense solution at ~28,000 radii and fails its own check above 1e-5 (the achieved error is in the metadata). Everything else (energy density, c_s^2, couplings, Phi', lambda', m') is computed by the consumer from the state and `eos_bsk26.csv`. The integrator runs at rtol 1e-6 with theta_c located to 1e-6; the original precision is available through `--rtol 1e-11 --shoot-tol 0`.
