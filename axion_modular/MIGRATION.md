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
5. **S1/N1 use Hermite tables** instead of the integrator's dense output. Interpolation error is bounded by S1 (1e-9) and
   the effect on the modes is measured by `tests/verify_vs_legacy_n1.py`.
6. **N1 reports only the stars it was given.** Check texts state how many values they cover; a full run needs all 12 stars.
7. **Inputs are fingerprinted by content.** `input_record` stores the `content_sha256` of equation artifacts (file hashes change with
   every run because of timestamps), and Dummy symbols are serialized by order of appearance. Found by checking that two runs of D1
   produced different hashes for identical equations.
8. **`PYTHONHASHSEED=0`** is fixed by the Makefile so that artifact hashes are reproducible.

## Known limitations (not fixed here)

* S1's `Star` and N1's `matter`/`rhs` and `modes.py` are hand-typed ODEs, as in the original. They are NOT compiled from
  D0/D2's artifacts (rule of `axion_ns/AGENTS.md` section 1, "direct numerical bridge"). The artifacts that would back them
  exist (`background`, `cowling_linear`); T has the pattern (`T/bridge.py`). Open item for the next iteration.
* The check-count regression (`tests/test_golden_checks.py`) compares against the committed `results/checks.json`.
* Python 3.12 was used; the legacy project declares >= 3.13.
