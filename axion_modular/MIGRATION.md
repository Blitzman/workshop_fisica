# Migration map: `axion_ns/` -> `axion_modular/`

`axion_ns/` is untouched and still runs (`uv run python -m axion_ns.run_all`). Every moved file keeps its logic; the
header of each new file says what it was moved from and what changed.

## Where each legacy file went

| legacy | new | what changed |
|---|---|---|
| `geometry.py` | `common/geometry.py` | copied unchanged |
| `compare.py` | `common/modulo.py` | only the import line |
| `model.py` (symbols, `eps_of`, `p_of`) | `common/notation.py` | definitions only; the notation contract |
| `model.py` (`derive_field_equations`, `background`) | `D0_derive_action/derive.py` | unchanged |
| `model.py` (`linearize`) | `D1_derive_radial/derive.py` | unchanged |
| `pipeline.py` | `common/checks.py` (`report`, `ok_field`) | the cached derivations are replaced by artifacts |
| `paper.py` | `R_paper_refs/paper1.py` | import line only |
| `nonradial.py` (`PAPER`) | `R_paper_refs/paper2.py` | verbatim |
| `nonradial.py` (derivation) | `D2_derive_nonradial/derive_cowling.py` | verbatim |
| `fullgr.py` (LD transcriptions) | `R_paper_refs/ld_sources.py` | `NR.x` -> `x` |
| `fullgr.py` (derivation) | `D2_derive_nonradial/derive_fullgr.py` | verbatim |
| `fullgr.py` (`NEW`, `EXT`) | `D2_derive_nonradial/ext_system.py` | `EXT` is built by `build_ext(LD, (8d))` from R's artifacts |
| `step1_background.py` | `D0/run.py` (2 checks) + `D1/checks_radial.py` (5 checks) | split by who owns the check |
| `step2_perturbations.py` | `D1/checks_radial.py` | verbatim logic |
| `step3_toy_model.py` | `D1/toy.py` | (3.3) from R's artifact; also returns the derived objects |
| `step5_nonradial.py` | `D2/checks_cowling.py` | `derived()` / `background_rules()` injected |
| `step7_fullgr.py` | `D2/checks_fullgr.py` | module imports replaced by `init(...)` |
| `bsk26.py` | `S1/bsk26.py` | unchanged |
| `axstar.py` (background half) | `S1/star.py` | unchanged |
| `axstar.py` (mode half) | `N1/modes.py` | `st.sol.sol` -> `st.y_at`, `g2` -> `st.g2` |
| `axstar.py` (`matter`, `rhs`) | `N1/star_table.py` | copies, validated against S1 in `tests/` |
| `step9_axion_cowling.py` | `N1/run_tables.py` (+ the 3 background-only checks in `S1/run.py`) | stars read from tables; works on any subset of nu_a |
| `step6_nonradial_numerics.py` | `N1/atom_levels.py` | tables from R |
| `step4_numerics.py` | `T/toy.py` | dps fixed at 30 (as in the full pipeline) |
| `step14_radial_branches.py` | `T/branches.py` | `f_toy` from `T/toy.py` |
| `figdata.py` | removed | each module writes its own JSON/CSV |
| `run_all.py` | `Makefile` | dependency order instead of one process |
| `report/build.py` (+ siblings) | `report/build.py` | generic, manifest-driven; the legacy builders stay in `report/` of the legacy tree |
| `step11-13`, `gmodes*`, `rotation*`, `magnetar*` | `archive_unvalidated/README.md` | not ported (retired hypotheses) |
| `step8`, `fmodes`, `step10`, `axgr`, `step25`, `radial`, `step15-17,19,20,22,26`, `step18,21,23,24` | not yet | modules N2, N3, N4, A |

## Deliberate differences

1. **Hidden dependencies became files.** `step10 <- step9`, `step18 <- step15`, `step21 <- step20`, `step25 <- step10` no
   longer exist inside the ported part: N1 reads S1's tables and R's tables, T reads D1's artifact, D1/D2 read D0's and R's.
2. **D0 also writes the background equations**, so D1 and D2 depend only on D0 and R, not on each other.
3. **`step1` was split** because two of its checks concern the derivation itself (D0) and the rest compare with the paper (D1).
4. **The paper's equations live in R** and derivation files cannot import them (enforced by `tests/test_independence.py`).
5. **S1 hands the star to N1 as a dense STATE table** (see item 10) instead of the integrator's dense output. S1 writes the
   interpolation error it achieved into the metadata and fails its own check if it exceeds 1e-5; the effect on the modes is measured
   by `tests/verify_vs_legacy_n1.py`. (An intermediate version used Hermite tables with derivatives, 2-9 MB per star.)
6. **N1 reports only the stars it was given.** Check texts state how many values they cover; a full run needs all 12 stars.
7. **Inputs are fingerprinted by content.** `input_record` stores the `content_sha256` of equation artifacts (file hashes change with
   every run because of timestamps), and Dummy symbols are serialized by order of appearance. Found by checking that two runs of D1
   produced different hashes for identical equations.
8. **`PYTHONHASHSEED=0`** is fixed by the Makefile so that artifact hashes are reproducible.
9. **S1 tolerances consolidated at 1e-6** (background ODE `rtol` and shooting precision). The original code used rtol 1e-11 and
   located theta_c to 1e-16 (about 66 integrations per star). `atol` (1e-30), `max_step` (0.25 km), the EOS table and the event
   thresholds are unchanged. The original precision is still available:
   `python -m S1_star_background.run --rtol 1e-11 --shoot-tol 0`. Single runs, one CPU core, Python 3.12:

   | | original precision | 1e-6 |
   |---|---|---|
   | S1, no-axion star + nu_a = 6 kHz | 109.6 s | 4.9 s |
   | star at 6 kHz: theta_c, M, R, Phi_shift (relative change) | - | 1.7e-9, 3.8e-8, 3.2e-8, 2.2e-8 |

   The tail ends earlier at 1e-6 (r_end 66.6 km against 119.5 km at 6 kHz: the decaying branch is trusted over a shorter
   distance); N1 uses `r_end` only as the radius where the tabulated background hands over to Schwarzschild with theta = 0. The 1 km
   step cap tried in the experiment (about 35x instead of 22x faster) is not adopted.
10. **Dense state table (S1) and one lookup per radius (N1).** S1 writes, per star, ln p, m, Phi, theta, theta' on 2,499 nodes
    (2,000 inside the star, clustered toward the surface with exponent 5; 500 uniform in the axion tail), about 0.3 MB. N1 evaluates
    the five columns with one cubic-spline lookup and computes everything else (energy density, c_s^2, couplings, Phi', lambda', m')
    from the state and the EOS table (itself one stacked polynomial). Decisions, each measured:
    * Tabulating the derived quantities against r was tried first and rejected: c_s^2 and friends inherit the EOS structure of the
      crust (pointwise error 1e-3 even with the best grid).
    * Two splines, split at the node r_k where the EOS is switched off: with one spline across r_k, ln p rings in the first tail
      intervals (p wrong by up to 20x where p ~ 1e-15). It was found because the validation floor hid it.
    * Surface clustering matters most: a uniform grid misses the right-hand side by 2e-2 within 2 m of R; clustering at both ends
      gives up to 1.8e-4; surface clustering with 1000 nodes gives <= 6e-6 and with 2000 nodes <= 2.3e-7 (stars at 3 and 10 kHz).
    * Cost does not depend on the number of nodes (binary search plus a polynomial), so 2,000 interior nodes are free.
    * Interpolation error written by S1 (state variables vs the integrator's dense solution, ~28,000 radii per star): at most 4.7e-6
      over the 12 stars (nu_a = 7), 2.1e-6 for the no-axion star.

    Measured effect (12 stars, one core): S1 25 s for everything; N1 mode search 203 s against 966 s with the Hermite tables
    (4.8x; 10-16 s per star, 37 s at 2.8 kHz and 43 s at 15 kHz against 167 s and 235 s). Against the previous route, Re nu agrees
    within 5e-8 and Im nu within 1e-5 for the damped modes (typically 1e-6). Against the ORIGINAL code at 6 kHz (full precision) Re nu
    agrees within 3.1e-8; Im nu of the weakly damped p1 mode (1.7055 Hz) is off by 5.2e-5 (the previous route: 5.5e-6; S1 at
    the original precision: 2.2e-5), i.e. 0.09 mHz. N1's paper checks pass (19/19), and against `fig_axcow.json` (4 digits) the
    shifts agree within 0.005 percentage points.

    What it does NOT show: no full-precision comparison with the original code for the other 11 stars (in particular the 2.8 kHz
    star next to the condensate threshold); the cause of the larger Im nu of p1 at 6 kHz is not identified (it may be the 1e-6 to
    5e-6 interpolation errors of Phi and theta', or the S1 tolerance; more nodes cost nothing and could be tried); and N1 keeps its
    copy of the `matter`/`rhs` formulas, checked against S1 in `tests/test_physics_crosscheck.py`.

## Known limitations (not fixed here)

* S1's `Star` and N1's `matter`/`rhs` and `modes.py` are hand-typed ODEs, as in the original. They are NOT compiled from
  D0/D2's artifacts (rule of `axion_ns/AGENTS.md` section 1, "direct numerical bridge"). The artifacts that would back them
  exist (`background`, `cowling_linear`); T has the pattern (`T/bridge.py`). Open item for the next iteration.
* The check-count regression (`tests/test_golden_checks.py`) compares against the committed `results/checks.json`.
* Python 3.12 was used; the legacy project declares >= 3.13.
