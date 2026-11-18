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
| `star_noaxion.csv, star_nua<nu>.csv + .meta.json` | table | dense STATE table: r, ln p, m, Phi, theta, theta' at about 7,200 nodes (6,198 inside: uniform in ln p + uniform in r + geometric near the centre; 1,000 in the fluid-free tail, where ln p = -inf); geometric units, km; metadata: M, R (at p = 1e-10 p_c), theta_c, Phi_shift, r_end, constants, integrator settings, interpolation error |
| `manifest.json `extra`` | data | `stars` (M, R, theta_c, rows, error, seconds), `condensate_exists`, `nu_a_list` |

Every run also writes `outputs/manifest.json` (checks, runtime, environment, sha256 of each artifact and of each input).

The star is integrated against ln p (start at r = 10 m from the centre series, stop at p = 1e-10 p_c, then the vacuum axion tail in r; `star_logp.py`). Consumers interpolate the columns with cubic splines (not-a-knot), as TWO splines split at the node `R_stop` (`r_k_km` in the metadata; the tail has no ln p). S1 measures exactly that interpolation against the exact state (`StarLogP.state_at_r`) at ~28,000 radii and fails its own check above 1e-5 (the achieved error is in the metadata; m, Phi and theta' are judged against max(|value|, 1% of the column maximum) because they vanish at the centre). Everything else is computed by the consumer from the state and `eos_bsk26.csv`. Defaults: ODE rtol 1e-8, theta_c located to 1e-12 (the shooting tolerance sets the valid axion tail, which N3 needs); the original precision is available through `--rtol 1e-11 --shoot-tol 0`.
