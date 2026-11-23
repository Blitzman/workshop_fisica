# D1_derive_radial: Radial linearization and simplified model (Paper I)

**Stage: production** (equations derived from the action and verified; frozen). Set by José; see the root `AGENTS.md`, Module Stages.

## Run
```bash
python -m D1_derive_radial.run            # from axion_ns/ ; `--help` lists the options
```
Outputs go to `D1_derive_radial/outputs/` (override with `--out`). Upstream artifacts are read from the sibling `outputs/` folders
(override with the `--*-dir` options), so this module runs on its own as soon as those files exist.

## Contract
Reads: D0 `field_equations`, `background`; R `paper1_equations`.

| artifact | kind | content |
|---|---|---|
| `radial_linear.eq.json/.pdf` | derived | Gtt, Gtr, Grr, KG, FLUID_r, BAR linearized (e^{-i w t} factored out) |
| `toy_model.eq.json/.pdf` | derived | Sec. IV: (4.2), (4.5), char. polynomial, (4.12), `f_4_13_tan` (matching determinant, symbols k_p, k_m, k_a, k_out, x_0, alpha, lambda2) |

Every run also writes `outputs/manifest.json` (checks, runtime, environment, sha256 of each artifact and of each input).

Checks of the old steps 1, 2, 3 (27) with the same names; consumed by T (bridge).
