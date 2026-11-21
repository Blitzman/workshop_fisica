# axion_ns

A reorganization of the original axion implementation into modules that run **independently**: modules never import each other. They exchange
files (equations as JSON + PDF, tables as CSV + metadata) through `common/`, a small library that holds no physics except
the notation contract. The original package is preserved in `archive/axion_ns/`.

```
                R_paper_refs ──────────────┬──────────────┬─────────────┐
 (paper equations & tables as data)        ▼              ▼             ▼
D0_derive_action ───────────────► D1_derive_radial   D2_derive_nonradial   N1_modes_cowling ◄── S1_star_background
(action → field equations,             │  (Paper I)    (Paper II + new)         (modes, atom)      (BSk26 + condensate,
 static background)                    ▼                                                           profile tables)
                                  T_toy_model                │
                                                             ▼
                                                    N2_modes_fullgr_noaxion (validation of the extended system vs LMI / Yoshida-Kojima)
                                  N3_modes_fullgr_axion (beyond Cowling, with the axion: reads D2 + S1 + R tables, N1 optionally)
                                  N4_modes_radial (radial modes: reads D0 + D1 equations, S1 tables and extra stars, R's Paper I tables)
                                       report (reads the manifests of all)
```
Every legacy step is now ported except step 17 (the Font et al. comparison, dropped by decision) and the pieces listed as not ported in `MIGRATION.md`.

## Run

```bash
pip install numpy scipy sympy mpmath matplotlib pytest      # + a LaTeX (pdflatex) for the equation PDFs
cd axion_ns
make R D0 D1 D2 S1 S1-extra N1 N2 N3 N4 T report test  # everything; D2 (about 8 min) and N2 (about 17 min of CPU, divided by the cores) are the long ones
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

* **The current version is the reference (2026-11-21).** No test compares with the legacy package any more. The migration-era comparisons
  (checks, data files, N1/N3/N4 against the legacy results) and their numbers are recorded in `MIGRATION.md`; the comparison scripts are
  in `archive/axion_ns_migration_scripts/` and the legacy gates were removed.
* **Test suite** (`make test`, about 10 s): S1 shooting and table nodes, the shared reader `common/star_table.py` against S1's formulas, N2,
  N3 and N4 table pieces, artifacts and independence. Each star is built once per session (`tests/conftest.py`). The interpolation error of
  S1's tables is checked only by S1 when it writes them. `python -m tests.compare_n3_n2_matrix` (about 1 minute, on request) checks that
  N3's interior matrix equals N2's with the axion off.
* **Reproducibility.** With `PYTHONHASHSEED=0`, all 10 equation artifacts have identical content hashes across independent runs
  (D1 and D2 were run twice). Inputs are recorded by the upstream *content* hash, so a rerun of an unchanged module does not
  make its consumers look stale; the Dummy symbols of `Subs` objects are serialized by order of appearance (SymPy's own index is
  random per process).
* **Independence** (`tests/test_independence.py`): no module imports another (shared code lives in `common/`); derivation files never import the paper references.
* **Equation -> numerics bridge** in T: the hand-coded Eq. (4.13) equals the matching determinant derived in D1 (lambdify).
* **Artifact round trip** with awkward expressions, and tamper detection.

## Measured runtimes (this sandbox: 1 CPU core, Python 3.12)

R 1.5 s, D0 4 s, D1 24 s, D2 420-460 s, T 107 s, tests 1 s. For comparison, the legacy `checks.json` (6-core i5) lists
steps 1-3 = 63 s, 5+7 = 543 s, 4+14 = 238 s.

Current S1 (rtol 1e-8, Brent match onto the decaying branch, dense state table): 24 s for the no-axion star and all 12 values of nu_a
(1.2-3.3 s per star, 2 cores, Python 3.14); `make S1-extra` about 2 minutes. The legacy r-method at the original precision (rtol 1e-11,
theta_c to 1e-16) took about 108 s per star. Current N1 on those tables: mode search 203 s for
the 12 stars (10-16 s per star, but 37 s at 2.8 kHz and 43 s at 15 kHz; previous route 966 s), 214 s for the whole `make N1`
including the atom levels. N1 reads each star through one spline lookup (about 9 microseconds per background evaluation against
about 60 for the previous Hermite/EOS route).

## S1 / N1 / N3 status (the ln p star and the N3 port)
* **S1 integrates the structure against ln p** (start at r = 10 m from the centre series, stop at p = 1e-10 p_c, then the vacuum axion tail in r) and writes tables whose nodes
  follow the solution (uniform in ln p, uniform in r, geometric near the centre; 7,197 rows, 0.8 MB). R and M are defined at p = 1e-10 p_c, about 9 m before the old EOS floor
  (no-axion star: M = 1.5359 Msun, R = 11.7069 km). **Shooting (2026-11-20):** the scalar is carried inside the star as u = ln(pi - theta), and theta_c is fixed by
  matching the nonlinear vacuum tail at R + 5/m_a onto the decaying branch of the linear exterior (inward Riccati, Brent in u_c); the tail is continued on that branch to
  R + 15/m_a, so the whole tabulated tail is valid and r_end no longer depends on the shooting tolerance. About 1.2-3 s per star, 24 s for the 13 stars (previous outward
  bisection to 1e-12: 75-114 s); stars above ~20 kHz, where theta_c rounds to pi, are now solvable (checked to 50 kHz). Against the previous tables: theta_c, R and M agree
  to 2e-8, 1e-8 km and 1e-9 Msun; the worst interpolation error is 3.5e-8. The previous shooting history (1e-6 / 1e-6 shortened the tail to 41 km at 11 kHz and moved the
  p2 damping by 4.7e-3) is in MIGRATION.md item 12 and the lab notebook.
* **N1** reads the new tables with a two-line change (the tail has no ln p): 19 of 19 paper checks pass, the coupled shifts agree with the legacy `fig_axcow.json` to 0.0056
  percentage points (36 values); the run takes 249 s.
* **N3** (legacy steps 10 and 25) is table-based: one right-hand-side call costs about 10 us instead of 520 us; the whole run takes 256 s on one core (legacy: 5,064 s on 6
  workers). Against the committed `results/fig_axgr.json` (migration gate, now removed): **111 of 114 comparisons agree within 1e-3**, every real part within 4.4e-6, the
  background-only shifts within 3.9e-6. The three above 1e-3 are damping values: p1 at 5 kHz (1.4e-3), p2 at 10 kHz (2.0e-3) and the no-axion p2 (3.2e-3 relative but
  7e-7 Hz absolute, the noise floor of the search). The dampings of modes near the axion threshold are reproducible only to about 1e-3: stopping the star 10 m earlier or later
  moves them by about 1e-3 in either direction (p2 at 10 and 11 kHz), and the original-precision run differs from the legacy by 1.1e-3 at 11 kHz; for p1 at 5 kHz the cause is
  not established (the cut does not move it).
* **Step 25** reproduces the legacy numbers. Three of its checks (H33 on Re nu, negative controls 3 and 4) fail in the committed legacy `exterior_convergence.json` as well and fail
  here with the same magnitudes; the H31 "above threshold" check was relaxed from `== 0` to `< 1e-12` (roundoff 2e-17).
* Fast tests: 26 passed, 1 skipped. The N3 interior matrix with the axion off equals N2's to 1.9e-15 (`python -m tests.compare_n3_n2_matrix`, manual).
* Not ported: the end-to-end LMI polytrope check of the legacy step 10 (needs a polytrope EOS table in S1; N2 validates the same equations).

## N4 status (radial modes)
N4 ports the legacy steps 15, 16, 19, 20, 22 and 26 (step 17, the Font et al. comparison, is not ported by decision) and, as extra parts (`--part analysis`), the radial analysis of steps 18, 21, 23 and 24 (about 100 s on one core; formerly planned as module A). The radial matrices are derived from D0's and D1's artifacts and
satisfy A = A_0 + w^2 A_1 exactly, so two coefficient tables are enough; a `det_entire` evaluation takes 0.11 s instead of about 2.7 s and the whole module about 35 minutes on one core (legacy: 23,205 s
wall on 6 workers). It needs S1's extra stars (`make S1-extra`, about 2 minutes, untracked). Against the committed legacy `results/radial_*.json` (migration gate, now removed): Table 1 and 2 of
Paper I within 9.5e-5 (real parts) and 1.1e-3 (damping), the F and H1 modes within 2.2e-8 and 1.8e-8, the avoided-crossing gap 49.1 Hz (legacy 49.07), all counts equal except three at the noise floor or
tied to the matching radius (step 20: deep region 0 against the legacy's 2 with no root located; old-lapse control 5 against 8; old-lapse deep 2 against 0). N4's checks: 22: 4/4, 26: 3 OK / 2 failures / 2 descriptive
(the two failures are failing in the committed legacy result too, with the same numbers), 15: 13/13, 19: 2/2, 20: 4/4, 16: 7/7 (6/7 on the owner's second computer, where the negative control 1 depends on the platform: accepted). Item 13 of `MIGRATION.md` lists the adaptations; the caveats above were reviewed and accepted by the project owner on 2026-11-18 (same item).

## N2 status
N2 (validation of the extended system without axion against Lindblom-Mendell-Ipser 1997 and Yoshida-Kojima 1997) was run in full by the project owner on two computers with identical
results: 450 s on 4 processes and 360 s on 6 (my estimate for one core: 15-25 minutes). Of its 7 checks, 5 pass. The Im w check failed only for l = 5 (Im w = 3.8e-10, +4.7% against LMI; l = 2-4 are
within 0.4%); by decision of 2026-11-18 l = 5 has a 5% tolerance, flagged in the note of the check. The convergence check (Re w moves 3e-5 on both computers; the legacy limit was 1e-5 and the legacy value 2e-6) has a limit of 5e-5 since the decision of 2026-11-18; with it and the l = 5 tolerance all 7 checks are expected to pass, but the validation has not been rerun with the new limits.
The legacy passed both; the results do not depend on the computer, so this is not platform noise. The only difference in the integration is how the right-hand side is evaluated (A y compiled with CSE
against the matrix compiled without it); the roundoff difference, amplified by the amplitude-minimum method, is the likely cause, but that has not been tested. What was measured before: the compiled
matrices equal the legacy ones to 1e-15, one Cowling evaluation agrees to 3e-14 (1.6x faster), the LMI f-mode agrees with LMI and the legacy code within the repo's tolerances (2.7x faster), and
loosening rtol breaks the amplitude method (`MIGRATION.md`, item 11).

## Known limitations

* S1's TOV + Klein-Gordon integrator and N1's perturbation equations are hand-typed, as in the original; they are **not** compiled
  from the equation artifacts (`axion_ns/AGENTS.md` “Scientific verification”). Only T has the bridge so far. See `MIGRATION.md`.
* The shipped S1/N1 outputs are the full set (12 values of nu_a plus the no-axion star), produced with the current defaults. Each
  star table is about 0.3 MB (2,499 rows), 3.7 MB for the 13 tables.
* Python 3.12 (the legacy `pyproject.toml` asks for >= 3.13); one CPU core.
* PDFs of the raw D0 equations contain terms too long to read; the readable forms are `background`, `radial_linear`, etc.
* `MIGRATION.md` maps every legacy file to its new place and lists the deliberate differences.
