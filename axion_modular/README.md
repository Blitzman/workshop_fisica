# axion_modular

A reorganization of `axion_ns/` into modules that run **independently**: modules never import each other. They exchange
files (equations as JSON + PDF, tables as CSV + metadata) through `common/`, a small library that holds no physics except
the notation contract. `axion_ns/` is untouched and still works.

```
                R_paper_refs ──────────────┬──────────────┬─────────────┐
 (paper equations & tables as data)        ▼              ▼             ▼
D0_derive_action ───────────────► D1_derive_radial   D2_derive_nonradial   N1_modes_cowling ◄── S1_star_background
(action → field equations,             │  (Paper I)    (Paper II + new)         (modes, atom)      (BSk26 + condensate,
 static background)                    ▼                                                           profile tables)
                                  T_toy_model          report (reads the manifests of all)
```
Not part of this layout yet: N2 (full-GR modes without axion), N3 (with axion), N4 (radial modes), A (radial analysis).

## Run

```bash
pip install numpy scipy sympy mpmath matplotlib pytest      # + a LaTeX (pdflatex) for the equation PDFs
cd axion_modular
make R D0 D1 D2 S1-witness N1-witness T report test          # everything except the long runs
make S1 N1                                                   # all 12 stars of Paper II Table II (long, see below)
python -m D1_derive_radial.run --help                         # every module has its own entry point and options
```
A module needs only the files of its upstream modules (see its `CONTRACT.md`), not their code or a previous run in the same
process. Each run writes `<module>/outputs/manifest.json` with checks, runtime, hashes of its inputs and outputs;
`python -m report.build` turns the manifests into `report/index.html` and flags **STALE** modules (an input changed since the run).

## The two kinds of artifact

| | produced by | read by | notes |
|---|---|---|---|
| `*.eq.json` | D0, D1, D2, R | D1, D2, T | SymPy expressions as a de-duplicated node table + provenance + content hash. Tampering is detected on load. |
| `*.pdf` | same | people | Typeset from the JSON (`python -m common.make_pdf --all` rebuilds all of them from the JSON alone). **Output only**: numerical code must load the JSON. Retyping from a PDF is exactly the transcription error the repo's rules forbid. |
| `*.csv` + `.meta.json` | S1, N1 | N1 | 17 significant digits; metadata with units, parameters, producer, sha256. |

Kinds of equation artifact: `derived` (from the action), `transcribed` (as typed in a paper; R), `candidate` (hand-assembled,
written only after its checks pass: `fullgr_extended`).

## What is checked, and against what

* **Same results as the original.** `tests/test_golden_checks.py` compares the checks of D0+D1 (old steps 1-3), D2 (5, 7),
  T (4, 14) and N1 (6) with the committed `results/checks.json`: same names (they embed the computed numbers) and same outcome.
  `tests/test_data_vs_legacy.py`: `fig_toy.json`, `radial_branches.json`, `levels.json` equal the legacy files.
* **N1 against the old modes** at nu_a = 6 kHz: `tests/test_n1_witness.py` (vs `results/fig_axcow.json`, 4 digits: all 9 numbers agree) and,
  at full precision against the original code (`make test-slow`). Measured at nu_a = 6 kHz, relative differences from the original:

  | S1 setting | Re nu (worst) | Im nu of p1 (1.7055 mHz) | Im nu of p2 |
  |---|---|---|---|
  | original precision (rtol 1e-11, Hermite 1e-9) | 1.6e-8 | 2.2e-5 | 4e-8 |
  | **current defaults (all 1e-6)**, measured with a Hermite cap of 64 | 3.6e-8 | 5.5e-6 | 1.2e-6 |

  For all 12 stars, N1 on the 1e-6 tables agrees with the legacy `fig_axcow.json` (4 digits) within 0.005 percentage points and
  passes its 6 paper checks. There is **no full-precision comparison with the original for the other 11 stars** (only the 4-digit
  one), and none for the 2.8 kHz threshold star at base precision. See `MIGRATION.md`, "Deliberate differences" 5 and 9.
* **Reproducibility.** With `PYTHONHASHSEED=0`, all 10 equation artifacts have identical content hashes across independent runs
  (D1 and D2 were run twice). Inputs are recorded by the upstream *content* hash, so a rerun of an unchanged module does not
  make its consumers look stale; the Dummy symbols of `Subs` objects are serialized by order of appearance (SymPy's own index is
  random per process).
* **Independence** (`tests/test_independence.py`): no module imports another; derivation files never import the paper references.
* **Equation -> numerics bridge** in T: the hand-coded Eq. (4.13) equals the matching determinant derived in D1 (lambdify).
* **Artifact round trip** with awkward expressions, and tamper detection.

## Measured runtimes (this sandbox: 1 CPU core, Python 3.12)

R 1.5 s, D0 4 s, D1 24 s, D2 420-460 s, T 107 s, tests 1 s. For comparison, the legacy `checks.json` (6-core i5) lists
steps 1-3 = 63 s, 5+7 = 543 s, 4+14 = 238 s.

S1 at the default tolerances (all 1e-6): about 5 s for the witness (no-axion star + nu_a = 6 kHz) and about 47 s for all 12 stars
(measured in the experiment; the final defaults raise the Hermite cap to 128 and were not rerun). At the original precision
(`--rtol 1e-11 --shoot-tol 0 --tol 1e-9 --k0 8`) one star takes about 108 s. N1 then dominates: 986 s for the 12 stars on one
core (50-60 s per star, but 167 s at 2.8 kHz and 235 s at 15 kHz); N1 witness 68 s.

## Known limitations

* S1's TOV + Klein-Gordon integrator and N1's perturbation equations are hand-typed, as in the original; they are **not** compiled
  from the equation artifacts (`axion_ns/AGENTS.md` section 1). Only T has the bridge so far. See `MIGRATION.md`.
* The shipped S1/N1 outputs are the **witness** (nu_a = 6 kHz + the no-axion star): the Table I/II/IV comparisons of N1 therefore
  cover one star, and the Table IV mass check is descriptive. `make S1 N1` computes all 12 (S1 about 1 min, N1 about 16 min on
  one core). Each star table is roughly 2-9 MB at the default tolerances (estimated from the row counts, which range from 7,000 to
  39,000), so the tables of the 12 stars (about 55 MB) are not shipped.
* The shipped S1 tables were produced before the tolerances were consolidated (rtol 1e-11, Hermite 1e-9); regenerate them with
  `make S1-witness` (or `make S1`) and then `make N1-witness`. Until N1 is rerun, `report/index.html` shows it as STALE.
* Python 3.12 (the legacy `pyproject.toml` asks for >= 3.13); one CPU core.
* PDFs of the raw D0 equations contain terms too long to read; the readable forms are `background`, `radial_linear`, etc.
* `MIGRATION.md` maps every legacy file to its new place and lists the deliberate differences.
