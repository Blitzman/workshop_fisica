# D0_derive_action: Action -> field equations and static background

**Stage: production** (equations derived from the action and verified; frozen). Set by José; see the root `AGENTS.md`, Module Stages.

## Run
```bash
python -m D0_derive_action.run            # from axion_ns/ ; `--help` lists the options
```
Outputs go to `D0_derive_action/outputs/` (override with `--out`). Upstream artifacts are read from the sibling `outputs/` folders
(override with the `--*-dir` options), so this module runs on its own as soon as those files exist.

## Contract
Reads: nothing.

| artifact | kind | content |
|---|---|---|
| `field_equations.eq.json/.pdf` | derived | Gtt, Gtr, Grr, KG, FLUID_r, BAR, ident_ax for general fields Phi, lambda, a, n, xi of (t, r) |
| `background.eq.json/.pdf` | derived | static limit of the above (TOV + Klein-Gordon form) |

Every run also writes `outputs/manifest.json` (checks, runtime, environment, sha256 of each artifact and of each input).

Symbols and function names are those of `common/notation.py` (the contract).
