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
| `axstar.py` (background half) | `S1/star.py` (physics) + `S1/star_logp.py` (integration against ln p) | `star.py` unchanged; the ln p integrator is new (item 12) |
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
| `fmodes.py` (EOS, `Star`) | `N2/polytrope.py` | `integrate` takes the compiled system and explicit rtol/atol, and counts right-hand-side calls |
| `fmodes.py` (systems) | `N2/systems.py` | read from D2's `fullgr_extended` and R's `paper2_equations`; `lambdify(cse=True)`; EXT and ZL from ONE symbolic solve |
| `fmodes.py` (modes) | `N2/modes.py` | same algorithms, explicit `rtol` and `atol` |
| `fmodes.py` (`YK97`), `step8_numerics.py` (`LMI97`) | `R_paper_refs/literature.py` | the published numbers as data |
| `step8_numerics.py` (compute, checks) | `N2/validation.py`, `N2/run.py` | the 7 checks with the same texts; fork pool with `--procs` |
| `step8_numerics.py` (axion background on a polytrope, m_eff, figures) | not ported | the schematic axion background and the figures |
| `axgr.py` (symbolic half) | `N3/systems.py` | `fullgr.EXT` from D2's artifact; the matrices are compiled as flat lists |
| `axgr.py` (`bgvals`) | `N3/background.py` | vectorised; Phi'' analytic instead of a finite difference |
| `axgr.py` (`basis`, `incoming`, `center`, `zerilli_amplitudes`) | `N3/modes.py` | right-hand sides from `N3/tables.py` (Laurent tables) |
| `axstar.py` (`exterior_L`, `muller`) | `N3/modes.py` | extracted verbatim from N1's copy |
| `step10_axion_gr.py` | `N3/validation.py`, `N3/run.py` | table-based; the LMI polytrope check and the comparison with N2's matrix are not ported (the latter is `tests/compare_n3_n2_matrix.py`) |
| `step25_exterior_convergence.py` | `N3/validation.py` (`step25`) | starting points from N3's own step 10 table |
| `radial.py` (derivation: `reduced_equations`, `first_order_system`, `interior_matrix`, `exterior_matrix`) | `N4/systems.py` | from D0's `background.eq.json` and D1's `radial_linear.eq.json` instead of `pipeline` |
| `radial.py` (`interior`, `outward`, `det`, `det_entire`, `exterior_psi`, `fluid_fraction`) | `N4/modes.py` | right-hand sides from `N4/tables.py` (A = A_0 + w^2 A_1 tabulated as r^2 A); `exterior_psi` takes the exterior flat function and beta |
| `radial.py` (`winding`, `winding_converged`, `muller`) | `N4/modes.py` | verbatim |
| `step15`, `step16`, `step19`, `step20`, `step22`, `step26` (`*_radial_star`, `*_radial_stability`, `*_avoided_crossing`, `*_damped_count`, `*_exterior_contour`, `*_paper1_tables`) | `N4/validation.py`, `N4/run.py` | stars from S1 (the 12 and the extra set); adaptations in item 13 |
| `step26_paper1_tables.py` (Tables 1-3 of Paper I) | `R_paper_refs/paper1_radial_tables.py` | the paper's numbers as data |
| `step17_radial_validation.py` | not ported | external validation against Font et al. 2002 (needs a polytrope star); dropped by decision |
| `step18,21,23,24` | not yet | module A |

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
9. **S1 tolerances consolidated at 1e-6** (superseded for the ln p star: the defaults are now rtol 1e-8 and shooting 1e-12, item 12) (background ODE `rtol` and shooting precision). The original code used rtol 1e-11 and
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
11. **N2: compiled with CSE, tolerances unchanged.** The numerical matrices are compiled from the artifacts with `lambdify(cse=True)`
    (A y compiled directly; EXT and ZL share one symbolic solve). Measured on the LMI star, l = 2, one core:
    * the extended system from D2's artifact equals the legacy `fullgr.EXT` (6/6 keys); the star is identical to 12 digits; the matrices
      agree to 1.1e-15 (LD) and 3.9e-16 (Cowling);
    * one Cowling evaluation: 3.2e-14 relative, 0.19 s -> 0.12 s; one amplitude evaluation: 1.66 s -> 0.55 s (3.0x), but the VALUES differ by
      1.3e-3 relative at the resonance (the amplitude is a near-cancellation: 5e-6 against a scale of 0.7);
    * the LMI f-mode through `full_mode`: legacy Re w = 1.201606 (+6.1e-6 from LMI), Im w = 3.943e-4 (-1.8e-3); new Re w = 1.201600
      (-9.6e-8), Im w = 3.946e-4 (-9.9e-4); 54 s -> 20 s (2.7x). Both pass the criteria (5e-5, 2%); they differ from each other by 6e-6 in Re w
      and 8e-4 in Im w: that is the noise of the method at rtol 1e-11;
    * **loosening rtol breaks the amplitude method**: rtol 1e-10 moves Re w by +7.7e-5 (> 5e-5), 1e-9 by +1.9e-3 (Im w +83%), 1e-8 gives Im w
      -72%. `atol` is not the lever: from 1e-12 to 1e-24 the amplitude wobbles by +-1% (resonance) / +-0.1% (off resonance) without converging,
      and the cost grows 170x (2 million rhs calls at 1e-24); at rtol 1e-10 no atol rescues it (value -23%..-34%).
    * Conditioning (measured): both interior solutions reach W ~ 1.6e6 and 9e5 at R(1 - 1e-7) and the X(R) = 0 combination cancels W by 5e5
      (2.5e5 at 1.45 u); this explains a W error of ~5e-6 and the resonance noise (4e-3 predicted, 1e-2 observed) but NOT the off-resonance noise
      (1e-5 predicted, 1e-3 observed): the rest is not understood. A better-conditioned matching (integrating the regular solution inward from
      the surface to a fitting point) would be a change of method that needs its own validation against LMI and Yoshida-Kojima; not done.
    * Cowling X(R) is well conditioned: rtol 1e-8 changes it by 1e-9 (2.2x faster); `--cowling-rtol 1e-8` is an option whose effect on the
      24 Yoshida-Kojima roots has not been checked.
    * The matrix A depends on w only through w^2 as A = A_-1/w^2 + A_0 + w^2 A_1 (fit residual 1e-15 at one radius; symbolic degrees agree);
      tabulating the three coefficient matrices against r would make one right-hand side call ~6 us instead of ~26 us (synthetic ceiling, with an
      accuracy that was not tested). Not implemented.

    The full validation (`make N2`) was run by the project owner on two computers with identical results (450 s on 4 processes, 360 s on 6): 5 of the 7 checks pass. The Im w check failed only for
    l = 5 (Im w = 3.8e-10, +4.7% against LMI; l = 2-4 within 0.4%); **decision (2026-11-18, project owner): l = 5 has a 5% tolerance**, flagged in the note of the check (the name of the check is
    unchanged so that the golden test still matches it). The convergence check still fails (Re w moves 3e-5 against a limit of 1e-5; legacy 2e-6): open. The legacy passed both checks and the results
    do not depend on the computer; the cause is probably the roundoff of the CSE-compiled A y against the legacy matrix product, amplified by the amplitude-minimum method (hypothesis, not tested).

12. **S1 integrates against ln p; the tables follow the solution; N3 tabulates the Laurent coefficients (measured, `tests/compare_s1_logp.py` and the N3 prototype).**
    * *Integrator.* `StarLogP` integrates (r, m, Phi, theta, theta') against s = ln p from r = 10 m (centre series, errors at 10 m: p 6.5e-13, theta 1e-13, theta' 6e-6 at 6 kHz) until
      p = 1e-10 p_c, then the vacuum tail in r with the same events and shooting as `Star`. At the original precision it reproduces the r-method star: theta_c to 8.4e-10 of delta, R and
      M at p = 1e-10 p_c to 3e-10 / 7e-11 (6 kHz); at the S1 defaults both agree with the reference to 1e-7..2e-6. No speed gain for axion stars (the cost is the tail, not the structure);
      the gain is in the tables. The structure is cut 9 m before the old EOS floor: it drops a layer of mass 3.5e-10 M.
    * *Nodes.* Uniform in ln p (r spacing 3.3 cm near the surface instead of 19 cm), uniform in r and 200 geometric nodes from 1 m to 1 km (the centre is a regular singular point), a
      clustered tail of 1,000 nodes. `StarLogP.state_at_r` returns the exact state at any r (centre series / Newton inversion on r^2 / tail): the first version inverted r instead of r^2,
      which converged slowly near the centre and left 1e-4 in theta' there (found by the validation, fixed). The validation floor for m, Phi and theta' (which vanish at the centre) is 1% of
      the column maximum; the kink where the series meets the ODE (10 m) is 3e-8 of the maximum in absolute terms.
    * *N3 tables.* A(r, w) = A_-1/w^2 + A_0 + w^2 A_1 exactly (residual 1.5e-16 at complex w); the three real coefficient matrices are evaluated exactly at the stored nodes and spline-interpolated;
      what is interpolated is r^2 A (plain A fails 0.56 at 1-10 m from the centre). Right-hand-side error against the exact path: <= 1e-11 in the last 50 m below the surface (the old grid: O(1) at 9 m),
      <= 1.4e-6 around the axion wall, <= 4.6e-15..1.1e-9 near the centre; one call 11 us against 519 us.
    * *Matching radius.* The legacy `r_match` (first |theta| < 1e-7, else the end of the table) lands at the end of the table, where the tail of theta_0 is corrupted by the shooting
      tolerance (1e-6). N3 uses R + 0.7 (r_end - R); for the 6 kHz f-mode the root moves 1.4e-8 (Re) and 3.4e-6 (Im) between 50 and 60 km, 2.3e-7 and 4.7e-5 at 30 km.
    * *Mode-level agreement (6 kHz f-mode, same theta_c).* tables vs the legacy exact path on the same star 4.0e-9 (Re) and 3.5e-8 (Im); vs the legacy on the r-integrated star 1.2e-8 and 3.1e-6;
      root in 4.3 s against 195 s; the root does not depend on rtol (1e-10 -> 1e-12 moves 5.7e-10 and 4.9e-8), unlike N2's amplitude minimum.
    * *S1 tolerances for N3 (ladder at nu_a = 11 kHz, p2 damping; legacy 5.032 Hz).* 1e-6 / 1e-6: 5.0557 (4.7e-3, valid tail 41 km); 1e-6 / 1e-10: 5.0387; 1e-8 / 1e-8: 5.0496;
      **1e-8 / 1e-12: 5.0332 (2.0e-4, tail 64 km)**; original 1e-11 / 1e-16: 5.0378 (1.1e-3). The shooting tolerance sets the valid tail and delta theta(R). The new defaults are rtol 1e-8
      and shoot 1e-12 (5-9 s per star, 75 s for the 13 stars).
    * *12-star gate* (`tests/n3_gate.py` against the committed `fig_axgr.json`, 6 digits): 111 of 114 comparisons within 1e-3; real parts within 4.4e-6, background-only shifts within 3.9e-6.
      **Above 1e-3: three damping values**: p1 at 5 kHz 0.327943 -> 0.328412 Hz (1.4e-3), p2 at 10 kHz 3.910590 -> 3.902849 (2.0e-3), no-axion p2 0.000210 Hz (3.2e-3 relative, 7e-7 Hz absolute: the Muller
      resolution is about 1e-5 Hz). *Cut of the star* (p_frac 1e-10 -> 1e-11, i.e. 10 m later): p2 at 10 kHz 3.9029 -> 3.9057 (towards the legacy), p2 at 11 kHz 5.0332 -> 5.0377 (away), p1 at 5 kHz unchanged
      (0.328412 -> 0.328405): the near-threshold dampings depend on the surface definition at the 1e-3 level in either direction, so they are not reproducible better than that; the cause of the p1 at
      5 kHz difference is not established (candidates not tested: the 10 m start, the tables near the wall, the shooting tail).
    * *Step 25* (`N3/validation.py`): reproduces the legacy values (control 3: 0.00151 Hz here, 0.00173 Hz legacy; control 4: 7.5e-7 against 8.5e-7 Hz; H33 Re: 2e-7 against 2e-6). Three of its checks are
      `False` in the committed legacy `exterior_convergence.json` too (H33 on Re nu against the 1e-7 criterion, negative controls 3 and 4: the n_far = 1 control does not move a mode whose damping is
      dominated by axion emission); they fail here identically. The H31 "above threshold" check asked for a difference of exactly 0 and now allows 1e-12 (roundoff 2e-17).
    * N1 on the new tables: 19 of 19 paper checks; coupled shifts within 0.0056 percentage points of the legacy `fig_axcow.json`.
    * N3 interior matrix with the axion off against N2's: 1.9e-15 (`tests/compare_n3_n2_matrix.py`).

13. **N4: the radial system from D0/D1 artifacts, with coefficient tables (assessment and runs, one core).**
    * *Structure.* The radial matrices are a 4 x 4 interior (12 non-zero entries, 3,639 operations, 986 after CSE) and a 2 x 2 exterior (3 entries, 68 -> 45), derived from D0's background
      equations and D1's linearised radial equations (20 s to compile). A(r, w) = A_0(r) + w^2 A_1(r) EXACTLY (fit residual 2.5e-16 interior, 1.1e-16 exterior; no 1/w^2 term, unlike N3), so
      only two coefficient matrices are tabulated (r^2 A, on the stored S1 nodes). Right-hand-side error against the exact arithmetic: <= 3.4e-11 near the centre, 1.1e-11 in the last 50 m,
      3.5e-13 in the tail, <= 7e-8 elsewhere. The entries of A grow to 1e13 towards the surface (as W in N2); the product A y keeps its accuracy.
    * *Cost.* In the legacy a `det_entire` evaluation made about 10,000 `bgvals` calls (83% of the time in `bgvals`, 11% in A): 203 us per call; with the tables 7.2 us per call and 0.11 s per
      evaluation instead of about 2.7 s. The legacy contour counts were small (step 20: 5,399 points; step 16: about 5,000 plus 3 x 459 axis points), so the whole N4 takes about 35 minutes on one
      core against the legacy 23,205 s wall on 6 workers.
    * *Modes (nu_a = 6 kHz, committed `radial_star.json`).* F = 3.4492123635 kHz (2.2e-8), H1 = 8.114832172 - 0.003649453i (Re 1.8e-8, damping 1.2e-5).
    * *Steps (time on one core; the results of the committed legacy in brackets).* 22: 27 s, 4/4. 26: 115 s, Table 1 within 9.5e-5, Table 2 real parts within 9.5e-5 and imaginary parts within
      1.1e-3; the two failing checks (Table 1 max deviation 1.09% at f0, 2.8 kHz; Table 2 dampings 11/22 within 25%) are `False` in the committed legacy result with the same numbers (known
      disagreements with the paper); Table 3, descriptive: Muller converges in 8/14 seeds (legacy 6/14). 15: 83 s, 13/13 (turning point p_c = 904.9, legacy 904.96; ladder n_eff
      2.012 ... 5.967 as the legacy). 19: about 300 s, 2/2, minimum gap 49.1 Hz at nu_a = 3.509 kHz (legacy 49.07 Hz at 3.5086). 20: 595 s, 4/4; index 2 in the main region and 0 in the thin strip as
      the legacy, fluid modes 8.1148 - 3.65i and 11.67 - 21.7i (4.3e-5); the counts that differ are in settings at the noise floor or tied to the matching radius: deep region 0 (legacy 2, with no
      root located: its H23 check was `False`, here `True`), old-lapse control 5 (legacy 8) and old-lapse deep 2 (legacy 0); a cavity resonance moves with r_m (80 km here, 117 km in the legacy).
      16: 871 s, 7/7 (index 1 beyond the maximum mass and 0 before it without axion; the axion control at p_c = 1000 has index 1 and a sign change at 0.58 kHz; 6 and 3.4 kHz have index 0).
    * *Extra stars.* S1 supplies, in an untracked `outputs_extra/` (about 20 MB, `make S1-extra`, 2 minutes): nu_a = 3.4, 3.5 and the 15 values of step 19, the nu_a = 6 star at p_c = 1000, and
      no-axion stars for p_c = 700..1100 in steps of 10 (58 stars, worst interpolation error 6.9e-7). The legacy built them on the fly.
    * *Adaptations.* (a) The turning point of step 15 is read from the fixed grid with cubic splines (the legacy used a bounded minimiser and a root finder on fresh stars). (b) Step 19 uses the 15
      values of nu_a the legacy ended up with after its run-time refinement. (c) The matching radius is R + 0.7 (r_end - R) (see item 12). (d) Robustness of Im nu: the start-radius variant
      r0 = 1e-4 km is r0 = 3e-3 km (the tables start at 1e-3 km) and the matching-radius variants are 0.8 and 0.6 of r_m. (e) The window of the retired "axion resonance" (step 15) is "above the
      threshold and damped": the legacy window was tied to its r_m (the resonance sits at 6.290 - 0.802i for r_m = 80 km and moves to 8.261 - 1.742i at 64 km). (f) `real_root` rejects poles of
      `det`: the exterior wave vanishes at r_m for some real frequencies, and with a shorter r_m than the legacy's a pole fell inside the bracket of the n = 2 gravitational-atom level.
      (g) `exterior_psi` uses g''(0) = -beta/4 instead of `axstar.g2(0.0)`.
    * *Not ported.* Step 17 (Font et al. 2002 and the Garching table): it needs a polytrope star; dropped by decision.
    * *Decision (2026-11-18, project owner).* The N4 caveats reported with the package were reviewed and accepted as they are: (1) the four differences from the committed legacy results (Table 2,
      nu_a = 8 kHz, f3: damping 1.1e-3 relative; step 20: deep region 0 against 2, old-lapse control 5 against 8, old-lapse deep 2 against 0; their explanation in terms of the noise floor and the
      matching radius was not proven); (2) the two checks of step 26 that fail here with the same numbers as in the committed legacy result (Table 1 maximum deviation 1.09% and Table 2 dampings
      11/22 within 25%); (3) the adaptations (a)-(g) above; (4) the multi-process path of N4 was not exercised when the package was built (single-process runs, in two pieces with a merged manifest); it was
      exercised afterwards on a second computer (next bullet); (5) the step 16 negative control 1 failing on that computer (next bullet): accepted as irrelevant.
    * *Second computer (project owner, 2026-11-18; 12 cores, Python 3.14.10, numpy 2.5.4, scipy 1.18.1, sympy 1.14.0).* S1-extra: 58 stars, worst interpolation error 6.8e-7, about
      10.5 s per axion star. N4 `--part fast` (6 processes) 286 s (compile 31 s; steps 22: 49 s, 26: 67 s, 15: 138 s); `--part all --procs 4`: step 16 451 s, 22 53 s, 26 70 s, 19 159 s,
      step 15 862 s (an unexplained slowdown of this serial step, 6x its time in the previous run; not reproduced; ignored by decision), whole run 1,917 s. Results identical to the
      single-process runs (same numbers and flags) for steps 22, 26, 15, 19 and 20; step 16 gives 6/7: the negative control 1 (the broken variant without the e^(i k r_m) factor and with coarse
      sampling) returned index 0 instead of a non-zero one. The index of that broken variant is an aliasing artefact that depends on the platform (-1.0 in the single-process run, +3.0 in the
      committed legacy); accepted as irrelevant. The adaptive sampling also differed slightly (nu_a = 3.4: 1,924 points against 1,546, one tile needing a third refinement), with the same
      index 0. `tests/n4_gate.py` on that computer: identical to the single-process run, 143 comparisons with the same 4 above tolerance (Table 2 nu_a = 8 kHz f3 damping 1.1e-3; step 20 deep
      region 0 against 2, old-lapse control 5 against 8, old-lapse deep 2 against 0); worst differences by file: Table 2 1.1e-3, turning point 6.5e-5, avoided-crossing gap 1.1e-6, fluid-mode
      damping 4.4e-5; the check outcomes that differ from the legacy are H23 (False -> True) and the step 16 control (True -> False).

## Known limitations (not fixed here)

* S1's `Star` and N1's `matter`/`rhs` and `modes.py` are hand-typed ODEs, as in the original. They are NOT compiled from
  D0/D2's artifacts (rule of `axion_ns/AGENTS.md` section 1, "direct numerical bridge"). The artifacts that would back them
  exist (`background`, `cowling_linear`); T has the pattern (`T/bridge.py`). Open item for the next iteration.
* The check-count regression (`tests/test_golden_checks.py`) compares against the committed `results/checks.json`.
* Python 3.12 was used; the legacy project declares >= 3.13.
