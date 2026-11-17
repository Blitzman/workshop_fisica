# N1_modes_cowling: Cowling modes and gravitational-atom levels

## Run
```bash
python -m N1_modes_cowling.run            # from axion_modular/ ; `--help` lists the options
```
Outputs go to `N1_modes_cowling/outputs/` (override with `--out`). Upstream artifacts are read from the sibling `outputs/` folders
(override with the `--*-dir` options), so this module runs on its own as soon as those files exist.

## Contract
Reads: S1 tables + manifest (part `tables`); R `paper2_tables.ref.json` (both parts).

| artifact | kind | content |
|---|---|---|
| `modes_cowling.csv + .meta.json` | table | per (nu_a, mode in f/p1/p2): no-axion, background-only and coupled (Re, Im) frequencies in kHz, M, R |
| `atom_levels.csv + .meta.json` | table | bound levels n = 3, 4 from the Schwarzschild exterior for 4 interior variants |
| `paper2_atom_levels.png, levels.json` | figure | Fig. 5-style plot and its data |

Every run also writes `outputs/manifest.json` (checks, runtime, environment, sha256 of each artifact and of each input).

Checks of the old steps 9 and 6, restricted to the stars present in S1's outputs.
