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
                                  T_toy_model                │
                                                             ▼
                                                    N2_modes_fullgr_noaxion (validation of the extended system vs LMI / Yoshida-Kojima)
                                       report (reads the manifests of all)
```
Not part of this layout yet: N3 (full-GR modes with axion), N4 (radial modes), A (radial analysis).

## Run

```bash
pip install numpy scipy sympy mpmath matplotlib pytest      # + a LaTeX (pdflatex) for the equation PDFs
cd axion_modular
make R D0 D1 D2 S1 N1 N2 T report test      # everything; D2 (about 8 min) and N2 (about 17 min of CPU, divided by the cores) are the long ones
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
* **N1 against the old results.** All 12 stars, current outputs: vs the legacy `results/fig_axcow.json` (4 digits) the shifts agree
  within 0.005 percentage points (`tests/test_n1_witness.py` checks nu_a = 6 kHz); N1 passes its paper checks (19/19 with the atom
  part). At nu_a = 6 kHz against the ORIGINAL code at full precision (`make test-slow`), relative differences:

  | S1 / N1 setting | Re nu (worst) | Im nu of p1 (1.7055 Hz) | Im nu of p2 (12.77 Hz) |
  |---|---|---|---|
  | S1 at the original precision (rtol 1e-11), Hermite CSV tables | 1.6e-8 | 2.2e-5 | 4e-8 |
  | S1 at 1e-6, Hermite CSV tables (previous route) | 3.6e-8 | 5.5e-6 | 1.2e-6 |
  | **current: S1 at 1e-6, dense state table** | 3.1e-8 | 5.2e-5 | 1.6e-6 |

  The Im nu of the weakly damped p1 mode is the delicate quantity (the damping changes by about 0.09 mHz); the real parts are stable.
  Against the previous route, over all 12 stars: Re nu within 5e-8, Im nu within 1e-5 for the damped modes (typically 1e-6; 5e-5 for
  p1 at 6 kHz and 7e-5 for p2 at 8 kHz). There is **no full-precision comparison with the original for the other 11 stars**. See
  `MIGRATION.md`, "Deliberate differences" 5, 9 and 10.
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

Current S1 (all 1e-6, dense state table): 25 s for the no-axion star and all 12 values of nu_a (about 2 s per star, measured). At the
original precision (`--rtol 1e-11 --shoot-tol 0`) one star took about 108 s. Current N1 on those tables: mode search 203 s for
the 12 stars (10-16 s per star, but 37 s at 2.8 kHz and 43 s at 15 kHz; previous route 966 s), 214 s for the whole `make N1`
including the atom levels. N1 reads each star through one spline lookup (about 9 microseconds per background evaluation against
about 60 for the previous Hermite/EOS route).

## N2 status
N2 (validation of the extended system without axion against Lindblom-Mendell-Ipser 1997 and Yoshida-Kojima 1997) was built and
smoke-tested (`make N2-quick`, LMI star only) but **its full validation (`make N2`, 7 checks) has not been run by the author**: run it
before relying on it. What was measured: the compiled matrices equal the legacy ones to 1e-15, one Cowling evaluation agrees to 3e-14
(1.6x faster), the LMI f-mode agrees with LMI and the legacy code within the repo's tolerances (2.7x faster), and loosening rtol breaks
the amplitude method (`MIGRATION.md`, item 11).

## Known limitations

* S1's TOV + Klein-Gordon integrator and N1's perturbation equations are hand-typed, as in the original; they are **not** compiled
  from the equation artifacts (`axion_ns/AGENTS.md` section 1). Only T has the bridge so far. See `MIGRATION.md`.
* The shipped S1/N1 outputs are the full set (12 values of nu_a plus the no-axion star), produced with the current defaults. Each
  star table is about 0.3 MB (2,499 rows), 3.7 MB for the 13 tables.
* Python 3.12 (the legacy `pyproject.toml` asks for >= 3.13); one CPU core.
* PDFs of the raw D0 equations contain terms too long to read; the readable forms are `background`, `radial_linear`, etc.
* `MIGRATION.md` maps every legacy file to its new place and lists the deliberate differences.
