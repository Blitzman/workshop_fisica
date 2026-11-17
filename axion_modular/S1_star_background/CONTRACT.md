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
| `star_noaxion.csv, star_nua<nu>.csv + .meta.json` | table | r, y = (m, p, Phi, theta, theta') and y' at the nodes; geometric units, km; metadata: M, R, theta_c, Phi_shift, r_end, constants, Hermite error |
| `manifest.json `extra`` | data | `stars` (M, R, theta_c, rows, error, seconds), `condensate_exists`, `nu_a_list` |

Every run also writes `outputs/manifest.json` (checks, runtime, environment, sha256 of each artifact and of each input).

Consumers rebuild y(r) by cubic Hermite interpolation (N1/star_table.py). The grid is refined (up to 128 subdivisions of each integrator step) until the interpolant reproduces the integrator's dense solution within 1e-6 (the achieved error is in the metadata and S1 fails its own check if the target is not met). The integrator itself runs at rtol 1e-6 with theta_c located to 1e-6; the original precision is available through `--rtol 1e-11 --shoot-tol 0 --tol 1e-9 --k0 8`.
