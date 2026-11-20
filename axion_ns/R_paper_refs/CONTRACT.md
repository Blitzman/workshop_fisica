# R_paper_refs: Paper equations and tables as data

## Run
```bash
python -m R_paper_refs.run            # from axion_ns/ ; `--help` lists the options
```
Outputs go to `R_paper_refs/outputs/` (override with `--out`). Upstream artifacts are read from the sibling `outputs/` folders
(override with the `--*-dir` options), so this module runs on its own as soon as those files exist.

## Contract
Reads: nothing. The equations and numbers are typed in `paper1.py`, `paper2.py`, `ld_sources.py`, `tables.py` (moved verbatim from `axion_ns`).

| artifact | kind | content |
|---|---|---|
| `paper1_equations.eq.json/.pdf` | transcribed | Paper I (2510.16245): background (2.6a-d) keys `B(2.6a) p'`..., radial perturbations (3.3a-e) keys `P(3.3a) KG`... |
| `paper2_equations.eq.json/.pdf` | transcribed | Paper II (2610.37308): (8a)-(8d), keys as in the paper |
| `ld_sources.eq.json/.pdf` | transcribed | Lindblom-Detweiler as typed in 4 sources (`ZL/`, `BU/`, `KU/`, `ZH/`), minimal corrections (`ZL_FIX/`, `BU_FIX/`) and `METFIX/` (base of the extended system) |
| `paper2_tables.ref.json` | data | Paper II Tables I-IV (`nu_a_kHz_rows`, columns documented in the file) |
| `literature_modes.ref.json` | data | Lindblom-Mendell-Ipser 1997 (f-modes l = 2-5, Re and Im) and Yoshida-Kojima 1997 (f, p1, p2, l = 2-5, M/R = 0.1, 0.2; exact and Cowling), in units of sqrt(M/R^3); read by N2 |
| `paper1_tables.ref.json` | data | Tables 1-3 of Paper I (radial modes of the BSk26 star with the condensate): T1 f0..f3 (kHz, with the no-axion row), T2 (Re kHz, Im Hz), T3 (Re, Im kHz of f1a, f2a); read by N4 |

Every run also writes `outputs/manifest.json` (checks, runtime, environment, sha256 of each artifact and of each input).

Rule: nothing in a derivation module may import this module. Consumers read its artifacts.
