# S1_star_background: BSk26 star with axion condensate: profile tables

**Stage: development** (shooting changed on 2026-11-20 to the decaying-branch match; reproduces Paper I M and R). Set by José; see the root `AGENTS.md`, Module Stages.

## Run
```bash
python -m S1_star_background.run            # from axion_ns/ ; `--help` lists the options
```
Outputs go to `S1_star_background/outputs/` (override with `--out`). Upstream artifacts are read from the sibling `outputs/` folders
(override with the `--*-dir` options), so this module runs on its own as soon as those files exist.

## Contract
Reads: nothing (BSk26 analytic fits in `bsk26.py`).

| artifact | kind | content |
|---|---|---|
| `eos_bsk26.csv + .meta.json` | table | n [fm^-3], eps and p [MeV fm^-3] (monotone-in-p subset used by the solver) |
| `star_noaxion.csv, star_nua<nu>.csv + .meta.json` | table | dense STATE table: r, ln p, m, Phi, theta, theta' at about 7,200 nodes (6,198 inside: uniform in ln p + uniform in r + geometric near the centre; 1,000 in the fluid-free tail, where ln p = -inf); geometric units, km; metadata: M, R (at p = 1e-10 p_c), theta_c, delta_c, Phi_shift, r_end, r_match, shooting diagnostics, constants, integrator settings, interpolation error |
| `manifest.json `extra`` | data | `stars` (M, R, theta_c, delta_c, r_match, theta_match, decay_mismatch, trials, rows, error, seconds), `condensate_exists`, `nu_a_list` |

Every run also writes `outputs/manifest.json` (checks, runtime, environment, sha256 of each artifact and of each input).

The star is integrated against ln p (start at r = 10 m from the centre series, stop at p = 1e-10 p_c; `star_logp.py`). With the condensate the scalar is carried in the interior as u = ln(pi - theta), so the shooting parameter u_c = ln(delta_c) stays resolvable when theta_c rounds to pi (about 20 kHz and above). The nonlinear vacuum tail (TOV + Klein-Gordon with backreaction) is integrated in r to the matching radius r_m = R + 5/m_a, where theta_c is fixed by requiring theta'/theta to equal the decaying-branch logarithmic derivative of the linear exterior equation on Schwarzschild (inward Riccati integration; Brent in u_c, xtol 1e-12). Beyond r_m the tail is continued on that branch (m and Phi from the full vacuum equations) to r_end = R + 15/m_a, so the whole tabulated tail is on the decaying solution and r_end no longer depends on the shooting tolerance. Metadata add `delta_c`, `r_match_km`, `theta_match`, `decay_mismatch` (residual theta' - S theta at r_m divided by m_a |theta|), `match_lambda` and `tail_lambda`.

Consumers interpolate the columns with cubic splines (not-a-knot), as TWO splines split at the node `R_stop` (`r_k_km` in the metadata; the tail has no ln p). S1 measures exactly that interpolation against the exact state (`StarLogP.state_at_r`) at ~28,000 radii and fails its own check above 1e-5 (the achieved error is in the metadata; m, Phi and theta' are judged against max(|value|, 1% of the column maximum) because they vanish at the centre). Everything else is computed by the consumer from the state and `eos_bsk26.csv`. Options: `--shoot-tol` (Brent xtol in u_c), `--match` and `--tail` (in reduced Compton wavelengths 1/m_a), `--rtol`. The legacy outward-bisection solver in r (`star.Star`) is kept as the physics base class and as a reference.

Note for N3/N4: their matching radius `R + 0.7 (r_end - R)` follows `r_end`, which is now R + 15/m_a for every star (previously 12-14 /m_a, set by where the outward shooting turned). At a fixed matching radius the new and previous tables give the same modes to ~1e-7; see the lab notebook (2026-11-20) for the effect of the moved radius.
