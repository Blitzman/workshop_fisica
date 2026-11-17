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
  at full precision against the original code (`make test-slow`, measured): worst relative difference in Re nu **1.6e-8**
  (typically 1e-9 to 1e-10); Im nu of the weakly damped p1 mode (1.7055 mHz) differs by 2e-5 relative. This is the price of the
  CSV + Hermite hand-over (S1 bounds the interpolation error at 1e-9; the modes amplify it by about 10-15x, more for Im nu).
* **Reproducibility.** With `PYTHONHASHSEED=0`, all 10 equation artifacts have identical content hashes across independent runs
  (D1 and D2 were run twice). Inputs are recorded by the upstream *content* hash, so a rerun of an unchanged module does not
  make its consumers look stale; the Dummy symbols of `Subs` objects are serialized by order of appearance (SymPy's own index is
  random per process).
* **Independence** (`tests/test_independence.py`): no module imports another; derivation files never import the paper references.
* **Equation -> numerics bridge** in T: the hand-coded Eq. (4.13) equals the matching determinant derived in D1 (lambdify).
* **Artifact round trip** with awkward expressions, and tamper detection.

## Measured runtimes (this sandbox: 1 CPU core, Python 3.12)

R 1.5 s, D0 4 s, D1 24 s, D2 420-460 s, S1 witness 220 s (one star + no-axion; ~200 s per additional star), N1 witness 68 s,
T 107 s, tests 1 s. For comparison, the legacy `checks.json` (6-core i5) lists steps 1-3 = 63 s, 5+7 = 543 s, 4+14 = 238 s.

## Known limitations

* S1's TOV + Klein-Gordon integrator and N1's perturbation equations are hand-typed, as in the original; they are **not** compiled
  from the equation artifacts (`axion_ns/AGENTS.md` section 1). Only T has the bridge so far. See `MIGRATION.md`.
* The shipped S1/N1 outputs are the **witness** (nu_a = 6 kHz + the no-axion star): the Table I/II/IV comparisons of N1 therefore
  cover one star, and the Table IV mass check is descriptive. `make S1 N1` computes all 12 (S1 ~200 s per star on one core;
  each star table is 4-9 MB, so they are not shipped).
* Python 3.12 (the legacy `pyproject.toml` asks for >= 3.13); one CPU core.
* PDFs of the raw D0 equations contain terms too long to read; the readable forms are `background`, `radial_linear`, etc.
* `MIGRATION.md` maps every legacy file to its new place and lists the deliberate differences.
