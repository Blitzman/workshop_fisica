# D2_derive_nonradial: Cowling and beyond-Cowling equations (Paper II + new work)

## Run
```bash
python -m D2_derive_nonradial.run            # from axion_modular/ ; `--help` lists the options
```
Outputs go to `D2_derive_nonradial/outputs/` (override with `--out`). Upstream artifacts are read from the sibling `outputs/` folders
(override with the `--*-dir` options), so this module runs on its own as soon as those files exist.

## Contract
Reads: D0 `background`; R `paper2_equations`, `ld_sources`.

| artifact | kind | content |
|---|---|---|
| `cowling_linear.eq.json/.pdf` | derived | (8a)-(8d) linearized from the action in Cowling, plus `Trt_ax` (axion energy flux) |
| `fullgr_linear.eq.json/.pdf` | derived | 11 linearized equations with the Regge-Wheeler perturbed metric (Ett, Etr, Err, Etth, Erth, KG, FLUID_r, FLUID_th, BAR, Eang, Etl) |
| `fullgr_extended.eq.json/.pdf` | candidate (verified) | extended system LD + axion; written ONLY if every equation passes its check |

Every run also writes `outputs/manifest.json` (checks, runtime, environment, sha256 of each artifact and of each input).

Checks of the old steps 5 and 7 (52). Derivation files never import R; R is read by the checks and by the assembly of the candidate.
