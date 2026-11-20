# T_toy_model: Toy-model numerics (Paper I Fig. 3 and radial branches)

## Run
```bash
python -m T_toy_model.run            # from axion_ns/ ; `--help` lists the options
```
Outputs go to `T_toy_model/outputs/` (override with `--out`). Upstream artifacts are read from the sibling `outputs/` folders
(override with the `--*-dir` options), so this module runs on its own as soon as those files exist.

## Contract
Reads: D1 `toy_model.eq.json` (bridge check only). The numerics need no input.

| artifact | kind | content |
|---|---|---|
| `toy_qnm_landscape.png, fig_toy.json` | figure | Fig. 3 of Paper I and the variant linearized at theta = pi |
| `radial_branches.json` | data | threshold, branches, avoided crossing, instability, winding numbers (old step 14) |

Every run also writes `outputs/manifest.json` (checks, runtime, environment, sha256 of each artifact and of each input).

Includes the numerical bridge: hand-coded Eq. (4.13) vs the equation derived in D1, at random complex points.
