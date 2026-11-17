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
5. **S1/N1 use Hermite tables** instead of the integrator's dense output. S1 writes the interpolation error it achieved into
   the metadata and fails its own check if it exceeds the tolerance (default 1e-6); the effect on the modes is measured by
   `tests/verify_vs_legacy_n1.py` (see item 9).
6. **N1 reports only the stars it was given.** Check texts state how many values they cover; a full run needs all 12 stars.
7. **Inputs are fingerprinted by content.** `input_record` stores the `content_sha256` of equation artifacts (file hashes change with
   every run because of timestamps), and Dummy symbols are serialized by order of appearance. Found by checking that two runs of D1
   produced different hashes for identical equations.
8. **`PYTHONHASHSEED=0`** is fixed by the Makefile so that artifact hashes are reproducible.
9. **S1 tolerances consolidated at 1e-6** (background ODE `rtol`, shooting precision, Hermite table tolerance). The original code
   used rtol 1e-11 and located theta_c to 1e-16 (about 66 integrations per star). `atol` (1e-30), `max_step` (0.25 km), the EOS
   table and the event thresholds are unchanged. The original precision is still available:
   `python -m S1_star_background.run --rtol 1e-11 --shoot-tol 0 --tol 1e-9 --k0 8`.

   Evidence (single runs, one CPU core, Python 3.12):

   | | original precision | 1e-6 defaults |
   |---|---|---|
   | S1, witness (no-axion + nu_a = 6 kHz) | 109.6 s | 4.9 s |
   | S1, 11 other stars | not run (about 108 s each expected) | 42 s |
   | star at 6 kHz: theta_c, M, R, Phi_shift (relative change) | - | 1.7e-9, 3.8e-8, 3.2e-8, 2.2e-8 |
   | N1 modes at 6 kHz vs the original code | Re nu 1.6e-8, Im p1 2.2e-5 | Re nu 3.6e-8, Im p1 5.5e-6, Im p2 1.2e-6 |
   | N1, 12 stars, vs legacy `fig_axcow.json` (4 digits) | - | within 0.005 percentage points; 6/6 paper checks pass |

   What it does NOT show: a full-precision comparison with the original for the 11 other stars, in particular the 2.8 kHz star
   next to the condensate threshold, where theta_c = pi - delta is most sensitive; and the 1 km step cap (configuration B in the
   experiment) was only tried at 6 kHz and is not adopted. The Hermite cap was raised from 64 to 128 subdivisions because the
   7 kHz star reached only 1.6e-6 at 64 (above the 1e-6 target); **the raised cap was not run**. With the integrator at 1e-6 the
   table ends earlier (r_end 66.6 km against 119.5 km at 6 kHz, because the decaying branch is trusted over a shorter distance);
   N1 uses `r_end` only as the radius where the tabulated background hands over to Schwarzschild with theta = 0.
   N1's cost does not drop (it is dominated by background lookups inside its own right-hand side), so the end-to-end gain is about 2.7x.

## Known limitations (not fixed here)

* S1's `Star` and N1's `matter`/`rhs` and `modes.py` are hand-typed ODEs, as in the original. They are NOT compiled from
  D0/D2's artifacts (rule of `axion_ns/AGENTS.md` section 1, "direct numerical bridge"). The artifacts that would back them
  exist (`background`, `cowling_linear`); T has the pattern (`T/bridge.py`). Open item for the next iteration.
* The check-count regression (`tests/test_golden_checks.py`) compares against the committed `results/checks.json`.
* Python 3.12 was used; the legacy project declares >= 3.13.
