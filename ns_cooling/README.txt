CURRENT RATE SELECTION
----------------------
The driver defaults to Friman-Maxwell modified Urca, using the existing
modified_urca_emissivity implementation and its superfluid factors.
Only MURCA changes: nucleon bremsstrahlung remains Bottaro-Caputo-Fiorillo.
No factor-of-20 experiment is enabled.
--murca bottaro selects the alternative MURCA calculation.
--old-rates is retained for compatibility and selects ALL standard rates,
including bremsstrahlung; it cannot be combined with --murca bottaro.
The startup message reports both selections. Each run retains its callable
registry snapshot for consistent energy diagnostics and process monitoring.
Direct calls to run_cooling still use the currently registered processes;
these driver defaults do not silently override a caller's custom registry.

DEFAULT MODEL: BSk26, 1.4 Msun, corrected Gudmundsson envelope.

VERSION 16 — late time stepping and neutrino dominance
-----------------------------------------------------
Default stopping condition is now L_gamma at infinity < 1e31 erg/s,
checked at every accepted state. The first state below threshold is retained.
No final-age or minimum-temperature limit is enabled by default.
--final-years and --min-temperature remain optional early-stop limits;
--min-luminosity sets the luminosity threshold.

Before 1e5 yr the existing startup, adaptive controller and 1000 yr cap remain.
After each accepted step reaching 1e5 yr, the next proposed dt is 1.05 times
that accepted dt, without the old cap. Trials exceeding the 5% maximum cell
relative temperature change are still rejected and retried with smaller dt.
This is an adaptive accuracy safeguard, not a mathematical stability proof.

At t=0 and every accepted step we integrate each registered neutrino process:
L_nu,infinity = sum Q_local * proper_cell_volume * exp(2 Phi).
Only accepted states enter the monitor. The dominant process is the largest
integrated luminosity, not the locally dominant emissivity or total cooling sink.
The saved registry snapshot defines process order. All-zero rates give 'none'.
Modified Urca and other registered groups retain their current grouping.

neutrino_processes.csv stores times, dominant process and all process luminosities.
cooling_neutrino_dominance.png is always saved, including with --no-animation.
It has one panel per star, avoiding ambiguous overlapping backgrounds. Shading
transitions use logarithmic midpoints of the accepted samples (no extra physical
precision implied). Red shading replaces neutrino shading wherever L_gamma,infinity exceeds the
sum of all neutrino luminosities. The CSV still records the leading neutrino process. The existing animation is retained alongside this final figure.

Validation of v16:
- Full BSk26 1.0, 1.4, 1.8 Msun runs stop at the first accepted luminosity below
  1e31 erg/s; final ages 322282, 227926, 163575 yr, respectively.
- All adaptive accepted steps satisfy the 5% maximum relative change.
- Late proposals grow by 1.05 except when rejected and reduced.
- 1.4 Msun history up to the old 1e5 yr stopping point is bitwise unchanged.
- Sum of monitored L_nu agrees with an independent total-emissivity evaluation
  at all accepted states for all three masses (relative tolerance 1e-12).
- Halving temperature-change controls changes late L_gamma by at most 0.30%
  for 1.4 Msun over the common time interval after 1e5 yr.
- The full pytest suite was not run: pytest is unavailable in this environment.

Older documentation follows; stopping defaults above supersede older defaults.

COOL1D - neutron-star cooling in general relativity, Python / NumPy, B = 0

V15 REVIEW FIXES (2026-11-05)
- Every CoolingRun records the ordered neutrino-process and heat-capacity callable
  snapshots present at the start of the run. energy_budget evaluates those saved
  functions directly, without consulting or changing the current registries.
  A run can therefore be diagnosed after leaving updated_rates(), or after a
  comparison has selected other rates. Runs without snapshots are rejected.
- updated_rates() and compare() restore the exact previous neutrino registry,
  including custom processes and order, also after exceptions or nested contexts.
  --old-rates explicitly selects the standard rates even in a reused interpreter.
- The fallback transport effective masses use the actual neutron and proton
  densities (nb*Yn, nb*Yp), including muonic compositions. The EOS provider still
  supplies its own Landau masses. The BSk26 default cooling curves are unchanged.
- Snapshots preserve callable selection, not mutable globals used inside custom
  functions. Keep physical constants and captured mutable parameters unchanged
  between a run and its diagnostic. Registries must not change during integration.
- The muon transport approximation and time-step settings are unchanged.
- Verification: the three default BSk26 runs to 1e5 yr reproduce the v15
  luminosity curves, temperature profiles and energy ledgers exactly. The
  updated-rate diagnostic before/after leaving its context now agrees exactly.
- Regression checks cover diagnostics after changing Q and Cv registries, nested
  contexts and exceptions, comparison restoration, and muonic fallback masses.


RUNNING
Open cooling_driver.py in VS Code and use "Run Python File in Terminal".
Do not run it cell by cell in the Jupyter interactive window.
From a terminal:
    python -m pip install -r requirements.txt
    python cooling_driver.py
or:
    make cooling PYTHON=python3
By default this runs ONE BSk26 star of 1.4 solar masses with the corrected
Gudmundsson envelope (input_bsk26). Other inputs remain available:
    python cooling_driver.py --input input_bsk22
The driver defaults to Friman-Maxwell modified Urca and Bottaro-Caputo-Fiorillo
nucleon bremsstrahlung. To change MURCA only:
    python cooling_driver.py --murca bottaro
To select all standard rates, including bremsstrahlung:
    python cooling_driver.py --old-rates
The run prints which rates are used. Scripts that call run_cooling directly use the standard processes
of the registry unless they call emissivity_updated.use_updated_rates().
Default paths are relative to the driver; no PATH changes are needed.

ANIMATION
The evolution of all the stars is computed first; then a window opens with two
synchronized panels, in which every star has its own colour (legend in the left panel):
- L at infinity vs time, logarithmic scale: each curve grows up to the time of the frame.
- log10 T at infinity vs log10 density (temperature in K, density in CGS) at the
  time of the frame.
Pause/Resume and Restart buttons. The window stays open with the final state.
The frames are 240 times equally spaced in log t, from the first step to the end of
the longest run, so all the stars are shown at the same physical time (the profiles
are interpolated between the accepted steps, linearly in log t and log T). A star
whose run ended earlier (minimum temperature reached) keeps its final state.
Default 40 ms per frame (about 10 s).
    python cooling_driver.py --interval-ms 20
A graphical Matplotlib backend is required; do not force Agg. In VS Code run
the file in the terminal with your usual interpreter.

OUTPUT
With one star the files are written directly in results_cooling/; with several,
in one subfolder per star (named after the mass, for example results_cooling/1.40_Msun/),
plus results_cooling/cooling_curves.csv with the luminosity curves of all the stars
(columns star, time_years, L_infinity_erg_s).
results_cooling/cooling_curve.csv: time and L at infinity.
results_cooling/temperature_profiles.csv: density and Tinf for all steps.
results_cooling/profile_times.csv: correspondence between step and time.
results_cooling/time_steps.csv: every accepted step, its dt and the maximum
relative temperature change.
results_cooling/energy_budget.csv: with --check-energy (see below).
No static figures or videos are generated automatically.
To compute without opening a window: --no-animation.
--out changes the results folder; --input changes the input file;
--background / --save-background read / write the background table (see
BACKGROUND; with several stars --save-background writes one file per star, with the
name of the star added; --background reads a single star).

MODEL AND EVOLUTION
Cooled domain: from the centre to the base density (base_density, default 1e10 g/cm^3); the layers
above it are represented by the envelope (Tb-Ts) relation.
Initial Tinf is set in the input file (T_8 = 1e2, i.e. 1e10 K).
First 20 steps: dt = 0.01 yr, fixed, no adaptive control (backward Euler).
Afterwards: maximum relative T change between cells <= 5%. If it is exceeded,
the step is rejected and repeated from the previous state. If it is below 1%,
the next step is increased (maximum factor 2) aiming at a 2.5% change. Inside
the band dt is kept until 1e5 yr, with upper limit dt = 1000 yr. After 1e5 yr
the proposed step grows by 5%, with the same rejection control and no upper cap.
Safety lower limit 1e-12 yr.
If a finite final time is requested, the run ends at the first step reaching it; the last step is not shortened,
so the final time can exceed the requested one by up to one step.
Optional early stopping at a final time or a minimum local temperature
(both disabled by default; default stop is L_gamma,infinity < 1e31 erg/s):
    python cooling_driver.py --final-years 100000 --min-temperature 10000000
All these numbers are fields of CoolingConfig (config.py).
Line 5 of the input file is not used.
The luminosity curve starts after 0.01 yr; the initial profile is shown at t = 0.

STRUCTURE
Data flow: provider -> Background -> evolution (never the other way round).
background.py: Background (everything that does not evolve in time), Points,
  face-table reading/writing, derived quantities and the composition check.
input_model.py: reads the input file and provides the Background, built from an
  equation of state. superfluid.py: pairing gap models and critical temperatures.
eos_table.py / tov.py / mesh.py / stellar_structure.py: the star from an equation
  of state (see STELLAR MODEL FROM AN EOS). constants.py: physical constants.
state.py: LocalState, the temperature plus the background at the points of
  interest; shared by all the microphysics.
registry.py: the Registry class used to plug in new physics (see below).
heat_capacity / emissivity / emissivity_factors / conductivity /
  conductivity_core / conductivity_crust / effective_mass: microphysics
  (see MICROPHYSICS). microphysics.py: evaluates all the microphysics of a state.
emissivity_updated.py: updated modified Urca and nucleon bremsstrahlung (see UPDATED
  NEUTRINO EMISSIVITIES). compare_emissivities.py: cooling curves with both sets of rates.
config.py: CoolingConfig, the settings of the time integration.
thermal_geometry / thermal_matrix / thermal_solver: grid, coefficients and
  tridiagonal solver; one time step.
thermal_boundary / thermal_fluxes: surface condition; Tb-Ts fits and their
  analytic derivatives (three envelopes).
cooling.py: evolution loop, adaptive steps, luminosity.
energy_balance.py: energy-balance check (see ENERGY-BALANCE CHECK).
cooling_driver.py: command line, execution and CSV output.
cooling_animation.py: animated window, independent of the solver.
numerical_validation.py: input protection.
tests/: test suite (see TESTS).

BACKGROUND (everything that does not evolve in time)
The evolution sees one object, background.Background. It holds, at the N faces
(conductivity, geometry) and at the N-1 cells (heat capacity, emissivity):
density, Ye, Yn, Yp, Ymu, Xh (mass fraction in nuclei), A, Z, exp(Phi), the critical
temperatures of crust neutrons, core neutrons and core protons, the Fermi
momenta kF_n, kF_p, kF_e and kF_mu and two sets of effective-mass ratios (one for
the heat capacity and the emissivities, one for the conductivity); plus the radii,
exp(Lambda), the surface gravity g14 and the envelope model. All arrays are
validated and read-only; the fields and units are listed in the docstring of
background.py. Neither the solver nor the microphysics nor the diagnostics
read files: they only see the Background.
The composition is per baryon: Ye, Ymu, Yn, Yp are particle numbers per baryon, so
n_i = Y_i rho / 1.66e-24 (LocalState also provides n_n, n_p, n_e, n_mu). In the core
charge neutrality is Yp = Ye + Ymu and baryon-number conservation is Yn + Yp = 1; every
Fermi momentum comes from its own number density, so kF_p^3 = kF_e^3 + kF_mu^3.
In the crust the nuclei carry the charge, Ye = Xh Z/A, and Yn + Xh = 1.
background.check_composition verifies these relations (the EOS builder calls it):
a violation in the core (1e-5) is an error; in the crust a charge mismatch above 25%
is a warning, because EOS tables satisfy Ye = Xh Z/A to rounding except in a
boundary layer next to the crust-core interface (up to 8% in BSk22 and 15% in BSk26,
down to about half of the interface density), and to ~1e-4 in the inner crust
(interpolation of the product between rows).
The included provider reads the input file (input_model.load_backgrounds), builds the star
from the EOS package (see STELLAR MODEL FROM AN EOS) and evaluates the pairing gap models there.
The complete primary description of a Background is a FACE TABLE, 13 columns
    radius_km, rho, ye, yn, yp, xh, aa, zz, ephi, elam, tc_crust, tc_n, tc_p
and an optional 14th, ymu (zeros if absent), plus the scalars g14 and
envelope_model (CSV header lines "# g14 = ...").
    python cooling_driver.py --save-background background.csv    # write it
    python cooling_driver.py --background background.csv         # run from it
The second form ignores the structure and superfluidity lines of the input file; the
initial temperature is still taken from it (or from --initial-temperature).
Cell values are face averages; the Fermi momenta and effective masses are
derived from the table by background.derived_fermi_quantities.
Adapting a different source of background information, in increasing order of
work, without touching the evolution:
- a table from another code or equation-of-state routine: write the columns
  in the CSV format above, or build the dictionary and call
  Background.from_face_table(table, g14=..., envelope_model=...).
- a routine that returns effective masses or Fermi momenta itself: pass it as
  derived= (it receives rho, ye, yn, yp and returns kf_n, kf_e, mstar_n,
  mstar_p, transport_mstar_n, transport_mstar_p; kf_p is computed from yp if the
  routine does not return it).
- everything from a routine, point by point: build the two Points objects and
  the Background directly.

STELLAR MODEL FROM AN EOS
The star is built from a cold EOS in
the CompOSE format. Add to the input file (input_bsk26 is the default input of the driver, one 1.4 Msun star):
    eos = BSk26              folder eos/BSk26/ with the files of the package
    mass = 1.0, 1.4, 1.8     target masses [Msun], one star each (or: central_pressure = ... [dyn/cm^2])
    base_density = 1e10      baryon mass density of the base of the cooled domain [g/cm^3]
    mesh = C                 preset A, B, C or D; optional overrides: dlnp_core,
                             dlnp_crust, blend, dr_core_max_m, dr_crust_max_m, refine
    python cooling_driver.py --check-energy
The structure line of the input file is then ignored. Several values in mass or
central_pressure define several stars (same mesh, pairing models and envelope for all).
The EOS packages BSk22 and BSk26 (Brussels-Montreal unified EOS, see their READMEs)
are included. input_bsk26 gives the central pressures of 1.0, 1.4 and 1.8 Msun found
with the TOV solver (6.0766217926446417e34, 1.2541308424984512e35 and
2.7326703375705250e35 dyn/cm^2); input_bsk22 those for BSk22 (3.964648164e34,
7.940366565e34 and 1.627427981e35). Writing mass = ... instead solves for the central
pressure on the stable branch of M(Pc), between the minimum mass (0.084 Msun for BSk22)
and the maximum (2.26 Msun); a target outside it is refused with those limits
(8-12 s per star).
Each EOS is one folder with the CompOSE files: eos.nb, eos.thermo and eos.compo are
required; eos.micro (Landau effective masses and 1S0 gaps, core only) and eos.mr
(mass-radius curve, used by the tests) are optional. Only one-dimensional T = 0 tables
are supported, and the reader refuses inconsistent ones (grid, e + p = n_b mu_b, charge
neutrality, phase order). A pressure that decreases by less than 0.1% between two
rows (BSk26 has one, 0.037% at data row 276) is accepted: that row is omitted from all
interpolated quantities, with a warning, and the files are not modified; larger
reversals are rejected.
Pipeline (stellar_structure.py)
- The structure is integrated once (tov.py) from the center to the outer edge of
  the table, in x = ln(P/Pc); exp(Phi) is matched to the Schwarzschild exterior
  there. The edge is the first row of the table (7.8e5 g/cm^3 for BSk22), not a
  computed photosphere; the layers above it (a few meters) are neglected.
- The faces are placed from the center to the base (mesh.py) with exact faces on
  the crust-core interface and on the base. The interface is the first core row
  of the table (the table brackets it between two rows).
- The EOS and the metric are evaluated at the faces and, independently, at the
  middle (in radius) of each cell, so no cell mixes crust and core values.
- rho is 1.66e-24 g * n_b, the convention of the microphysics (n = rho/1.66e-24),
  a baryon REST-MASS density (not the mass-energy density e/c^2 of the TOV equations, which
  differs from it by about 1 %), so the Fermi momenta reproduce the tabulated baryon density. The composition
  (including the muons) and its checks are described in BACKGROUND.
- Effective masses: the Landau masses of the EOS in the core (used by the heat
  capacity, the emissivities and the conductivity); in the crust, where the EOS
  has none, the tables of effective_mass. The critical temperatures come from
  the gap models of the input file.
- The surface: g14 is the gravity at the outer edge. The luminosity at infinity is
  conserved across the envelope, so the luminosity and the heat flux leaving the base
  carry the factor (R_s e^Phi_s / R_b e^Phi_b)^2 (1.036 for BSk22 at 1.4 Msun);
  Background.surface_factor, 1 for a table without surface values.
- --save-background writes the face table of such a star for inspection; it cannot
  be read back (its cells and masses come from the EOS, not from the faces).
Meshes (faces in ln P; Pc = 7.94e34 dyn/cm^2, 1.4 Msun BSk22)
    A  uniform, step 0.101                                            161 faces
    B  crust step 0.05                                                274 faces
    C  B + core cells <= 300 m                                        290 faces
    D  twice as fine as C                                             581 faces
(refine = 2 halves every spacing and cap, for convergence studies.)
THE DIRECT URCA EDGE MUST BE RESOLVED. The process is strong and switches on
sharply at the threshold, so its volume depends on the cells that contain the edge.
For the 1.4 Msun BSk22 star (edge at 5.3 km) L_inf relative to mesh D is 16-17% too low
on mesh A from 100 yr, and within 0.1% on mesh C. With Direct Urca switched off all
meshes agree within 0.5%. The builder warns when the edge lies in cells wider than 350 m
(a heuristic calibrated on this star). The warning is not a guarantee of convergence:
for a star with Direct Urca use mesh C or finer and compare meshes (refine = 2).
Validation (tests): an n = 1 polytrope reproduces the analytic radius and mass up to the
post-Newtonian correction (proportional to the compactness); the mass, hydrostatic and
metric equations hold on the profile; the mass-radius curve matches eos.mr within 0.13%
in radius and 0.03% in M_max; the Direct Urca threshold star has 1.148 Msun for BSk22
(the description of the EOS quotes 1.151 Msun).
Limits: no layers above the first row of the table; one nuclear species per row;
no 3P2 neutron gap or crust effective masses in the EOS; the interface is known
only to the spacing of the table.

MICROPHYSICS
Heat capacity (heat_capacity.py). Registered terms, summed:
- electrons and muons: degenerate relativistic free lepton gases,
  C = k_B^2 T E_F k_F / (3 (hbar c)^2), E_F = sqrt((hbar c k_F)^2 + m^2 c^4),
  with k_F from the density of each species (kf_e, kf_mu) and m_mu c^2 = 105.66 MeV.
  Muons exist only in the core (zero where Ymu = 0). The electron term is a fit with
  rounded constants of the same expression (agreement 0.3%). At the center of a
  1.4 Msun star the muons contribute 85% (BSk22) and 79% (BSk26) of the electron term.
- neutrons and protons: Landau quasiparticles, C = k_B^2 T m* k_F / (3 hbar^2), with
  the Landau effective mass m* = (m*/m) m_N of the background and k_F from the density
  of each species. The Landau mass is k_F / v_F, so it already contains the
  relativistic and the interaction corrections and no further relativistic factor
  appears. Reduced by superfluidity: singlet crust neutrons, triplet core neutrons,
  singlet core protons (heat_capacity_pairing_factor).
- ions: gas, liquid and Debye-solid regimes of the crust nuclei.
Neutrino emissivity (emissivity.py). Registered processes: modified Urca (neutron and
proton branches), neutron-neutron, neutron-proton, proton-proton and electron-proton
bremsstrahlung, Cooper pair breaking and formation (core and crust), direct Urca,
electron-ion bremsstrahlung, plasmon decay and pair annihilation. Direct Urca has an electron
and a muon channel, each allowed where its own threshold below is met. The standard modified Urca
treats the core leptons as electrons; the driver replaces it and the nucleon bremsstrahlung by the
updated rates (UPDATED NEUTRINO EMISSIVITIES below), which have the muon branches.
Conductivity (conductivity.py). Core electrons, muons and neutrons (conductivity_core.py) and
crust electrons (conductivity_crust.py). Composition per baryon in the core: Yp = Ye + Ymu,
Yn = 1 - Yp.
- Leptons (electrons and muons): kappa_l = C T / nu_l (n_l / n_0)^(2/3) v_l, n_0 = 0.16 fm^-3, with
  nu_l = nu_ll + nu_ln + nu_lp the lepton-lepton, lepton-neutron and lepton-proton collision
  frequencies (they share a screening wave number q, built from n_l, n_n, n_p and the nucleon
  effective masses) and v_l = p_F / mu_l the velocity factor of the finite mass: v = 1 for the
  electrons; the muons use the same expression with their own density and
  v_mu = p_F / sqrt(p_F^2 + m_mu^2). For ultra-relativistic leptons n_l / mu_l = n_l / p_F, which is
  the form of the electron expression. The nucleon superfluidity enters through the factors
  R (collisions), Z (screening) and R_C (neutron conduction) of conductivity_core.py.
- Approximations of the muon channel: each lepton is treated as the only lepton species in its
  collision frequencies and screening, so electron-muon collisions are not included and the
  channels are added as independent (the coupled electron-muon transport of Shternin and Yakovlev
  would give a smaller total); the collision frequencies of the muon are those of the electron
  expression for its own density; flux-tube scattering is not included.
- Neutrons: kappa_n from neutron-neutron and neutron-proton collisions, proportional to n_n.
- In the 1.4 Msun stars the muon conductivity is 0.5 (BSk26) and 0.6 (BSk22) of the electron one
  at the center and vanishes at the onset of the muons; the cooling curves change by less than 1%
  at ages of 1 yr to 1e5 yr (up to a few per cent where the curve is steepest).

UPDATED NEUTRINO EMISSIVITIES
emissivity_updated.py implements the modified Urca and the nucleon bremsstrahlung of
Bottaro, Caputo and Fiorillo (arXiv:2406.18640, JCAP 2024), a re-derivation of the
Friman-Maxwell framework: one-pion plus rho-meson exchange (m_pi = 138 MeV, m_rho = 769 MeV,
C_rho = 1.4, f_pi = 1) and Landau parameters (F0' = 0.7, G0 = 1.1, with f = pi^2 F / (2 m_N p_F(n))).
Closed forms used (natural units, one neutrino flavor unless stated):
- nn and pp bremsstrahlung, Eq. 3.11 with J of Eqs. 3.7-3.9 (p_F of the species), 3 flavors;
- np bremsstrahlung, Eq. 3.15 with the leading J^np of Eq. 3.16 (p_p << p_n), 3 flavors;
- modified Urca, neutron branch, Eq. 3.28 (triangular approximation): proportional to
  p_p and to the sum over the leptons of v_l = p_l / mu_l, so the muons contribute with their
  velocity and not as relativistic particles;
- modified Urca, proton branch, Eq. 3.39: (p_l + 3 p_p - p_n)^2 / mu_l for each lepton, open
  where p_l + 3 p_p > p_n (each lepton has its own threshold).
Compared with the standard rates: the rho-meson suppression, the exchange and interference
terms, the Landau-parameter structure of the amplitude, the proton Fermi momentum (not the
electron one) in the modified-Urca phase space, and the electron and muon channels together.
The short-range suppression factors of Friman and Maxwell (0.56 nn, 0.66 np, 0.70 pp, 0.68
modified Urca) are kept, as in the paper's comparisons; the nucleon effective masses enter as
overall powers, the Landau parameters use the bare mass. The pairing reduction factors are those
of the standard processes. The rates apply in the core; the crust keeps the standard neutron
bremsstrahlung, and electron-proton bremsstrahlung, Direct Urca, pair breaking and the other
processes are unchanged. Not included: the non-triangular modified-Urca integral (Eqs. 3.34
and 3.37, a four-dimensional numerical integral per point, which the paper finds 20-30% larger
for the neutron branch) and in-medium effects. A muon branch of Direct Urca is not part of
the paper.
Use (Friman-Maxwell MURCA and updated bremsstrahlung are the driver defaults):
    python cooling_driver.py --old-rates              # the standard rates instead
    from emissivity_updated import updated_rates      # context manager; or use_updated_rates()
    with updated_rates():                             # use_standard_rates() restores the others
        run = run_cooling(background, config)
Comparison of cooling curves of one star (standard rates, updated rates and variants that isolate
the muon channels, the modified Urca, the bremsstrahlung, the suppression factors and the size of
the non-triangular correction):
    python compare_emissivities.py --input input_bsk26 --star "1.40 Msun"
writes results_emissivities/emissivity_comparison.csv and .png and prints the luminosity ratios.
Validation (tests): the closed forms equal numerical integrals of the amplitudes of the paper
(Eqs. 3.3 and 3.13) to better than 1e-5; the bremsstrahlung reproduces the values of the paper's
Figs. 1 and 3 (for p_n = 1.75 fm^-1, T = 1e9 K, bare masses, 3 flavors: nn 4.3e19 erg cm^-3 s^-1, 10.5e19
without the rho meson; np 5.6e19, 6.7e19 with Landau parameters and no rho meson); the modified Urca
neutron branch relative to the standard rate in npe matter is 0.62 at p_n = 1.5 fm^-1 and 0.33 at
2.2 fm^-1; the lepton velocities, the thresholds, the T^8 and mass scalings and the registry
replacement are checked.

DIRECT URCA THRESHOLD
There is an electron channel (n -> p e nu) and a muon channel (n -> p mu nu). Each is allowed where
momentum conservation permits it in homogeneous matter: kf_n <= kf_p + kf_e for the electrons and
kf_n <= kf_p + kf_mu for the muons, with kf_p from the proton density. In beta equilibrium
mu_mu = mu_e, so kf_mu < kf_e and the muon channel opens at a higher density than the electron one.
Emissivity of each lepton: Q_l = 4.0e27 m*_n m*_p T9^6 (mu_l / E_ref) R_du, E_ref = 331.4 MeV, where
mu_l is the Fermi energy (sqrt((hbar c kf)^2 + m^2 c^4)) and mu_e / E_ref = (Ye rho / 2.656e14)^(1/3) is
the electron factor; where both channels are open and mu_mu = mu_e the emission doubles. The lepton
enters only through mu_l and its threshold. This was verified (tests/test_direct_urca_muons.py): the
spin-summed matrix element computed with explicit Dirac spinors has no mass term besides the lepton
energy, and the phase-space integral with exact kinematics gives Q = (457 pi / 10080) G^2 cos^2(theta_C)
(1 + 3 g_A^2) m*_n m*_p mu_l (k_B T)^6 for the electron and the muon alike (0.999 of the coefficient at
T = 0.1 MeV), and Q_mu / Q_e = mu_mu / mu_e. The sharp threshold is the usual approximation; within a
few MeV of it the thermal smearing of the Fermi surfaces matters (as for the electron channel).
In npe matter the electron threshold is Yp >= 1/9;
with muons the threshold depends on the composition: for BSk22 the threshold star is
1.148 Msun, against 0.883 Msun for Yp = 1/9. In BSk22 the threshold is at a central
pressure of 5.16e34 dyn/cm^2 (rho = 5.5e14 g/cm^3): the 1.0 Msun star does not reach it
(Pc = 0.77 of the threshold); in the 1.4 Msun star it is active out to 5.5 km (42%
of the radius, 16% of the mass) and in the 1.8 Msun star out to 8.2 km (63%, 46%).
In the default BSk26 stars both channels are closed everywhere: at the center
kf_n = 2.23, 2.40 and 2.61 fm^-1 against kf_p + kf_e = 1.72, 1.86 and 2.03 fm^-1 and kf_p + kf_mu =
1.52, 1.67 and 1.87 fm^-1 for 1.0, 1.4 and 1.8 Msun. In the BSk22 stars the muon channel is closed
at 1.0 and 1.4 Msun (at the center of the 1.4 Msun star kf_p + kf_mu = 2.094 against kf_n = 2.119 fm^-1) and
open at 1.8 Msun, in 7 of the 11 core cells where the electron channel is open.

UNITS AND ASSUMPTIONS
Reference density: every emissivity normalised by a reference density uses n0 = 0.16 fm^-3 of
Yakovlev et al. (2001), i.e. 2.656e14 g/cm^3 as a rest-mass density (constants.N0_FM3, constants.RHO0);
the pair-breaking rates use the actual Fermi momenta of the background, p_F/(m_N c) with m_N the neutron mass.
Not covered: the electron-proton bremsstrahlung (2.82e14 g/cm^3) and the fit constant 2.8e14 of the
electron-ion bremsstrahlung, which keep the values of their sources.
Microphysics: local temperature in K. Evolution: Tinf_8 = exp(Phi)*Tlocal/1e8.
Conductivity: CGS, converted to the internal coefficient with kappa*1e13*1e-40.
Time unit: 1 yr = 3.16e7 s (SECONDS_PER_YEAR in thermal_matrix.py), about
0.13% longer than the Julian year.
Crust: strongly degenerate electrons, point ions, Qimp = 0.1;
ni = rho*xh/(m_u*A). No thermal averaging, ion conductivity, nuclear form
factor or Umklapp freeze-out. The EXPINT function needed is built into
the same module. Source of the fit: Potekhin's electron-conductivity fit
(2019 version), COULAN 1999, TAUEESY 2016. Some constants are rounded to
single precision.

NUMERICAL CHOICES
Crank-Nicolson with a predictor-corrector (two microphysics evaluations and
two linear solves per step); backward Euler is available (see TIME INTEGRATION
SCHEME).
Emissivity derivative at fixed T/Tc, a deliberate approximation.
Volumes use the outer elam; the last cell is the base of the envelope.
Envelope exponent 1/1.82; the output Ts is computed after updating T.
Linf = 4*pi*R_cm^2*sigma*Ts_local^4*exp(2*Phi).
T/Tc is computed with the local temperature and Tc in Kelvin.
Masses, particle densities and Fermi momenta are precomputed.
Emissivities and core transport are vectorized; the crust conductivity is
evaluated cell by cell.

ENVELOPES
Last line of the input file:
    envelope_model = gudmundsson
Available values (no quotes): gudmundsson, fe, accreted.
The selector applies to the boundary condition of the matrix and to the Ts and
L computed after each step. If missing, gudmundsson is used; a wrong value
gives an explicit error. No other file needs to be modified.

Gudmundsson: Ts = 1e6*g14^0.25*(Tb8/1.288)^(1/1.82).
Fe and accreted: fits of Potekhin et al. (2003), ApJ 594, 404. Accreted is
that paper's accreted-envelope fit, not a pure-H model with configurable mass.
The three models use the same implicit linearization:
    x = Tb/1e8; F = sigma*Ts^4*1e-30
    cfluxb = dF/dx
    sfluxb = F - x*cfluxb
Each fit has its analytic derivative. For u = Tinf/1e8 and ephi = exp(Phi):
    offset = ephi^2*sfluxb; slope = ephi*cfluxb
so that offset + slope*u reconstructs ephi^2*F at the linearization point.
The matrix is evaluated at the old temperature; the output Ts and Linf are
evaluated on the updated profile. There is a single solve per step.
Tb and Ts are local inside the fits. Fe/accreted require zeta > 0; an invalid
input is rejected.

TIME INTEGRATION SCHEME
Default: Crank-Nicolson with a predictor-corrector (CoolingConfig.theta = 0.5,
predictor_corrector = True, fixed_steps_backward_euler = True).
Each step: (1) a backward-Euler predictor gives a first estimate of the new
state, using the coefficients (heat capacity, source linearization,
conductivity, boundary fit) at the old state; (2) the coefficients are
re-evaluated at the MIDPOINT of the old and predicted states and the
Crank-Nicolson corrector is solved from the old state. Cost: two microphysics
evaluations per step (1.96x the CPU time per step of backward Euler). The first
n_fixed_steps steps (the 20 fixed steps of 0.01 yr) use plain backward Euler:
it damps completely the fast outer-crust modes that this unresolved start would
otherwise seed and that Crank-Nicolson does not remove.
Plain backward Euler:
    python cooling_driver.py --theta 1 --no-predictor-corrector
or CoolingConfig(theta=1.0, predictor_corrector=False). Single-stage
Crank-Nicolson (theta = 0.5, predictor-corrector off) is also possible but
is first order.

Properties of the scheme, from the energy balance dE/dt = -(L_nu + L_gamma) with
E = sum V e^Phi int(Cv dT). "Resolved start" means dt_initial_years=1e-8 and
n_fixed_steps=0; references are the finest runs.
- Frozen coefficients limit Crank-Nicolson to first order for the nonlinear
  problem; midpoint coefficients give second order. With fixed steps and
  smooth physics (superfluidity off) the predictor-corrector converges at order
  1.97-2.07 while backward Euler and single-stage CN are first order.
- Full physics, resolved start: the energy-closure error converges at about
  second order. With fixed steps the observed order lies between 1 and 2, because
  the microphysics is not smooth (jump of the superfluid heat-capacity factor at Tc,
  validity thresholds of some neutrino processes).
- Without the backward-Euler start-up the predictor-corrector rings in the outer
  crust (rho 1e10-1e14 g/cm^3): upward temperature steps appear and grow after 2e4 yr.
- The first ~0.2 yr is under-resolved with the default start; this has no
  observable consequence (with a resolved start the results at ages >= 1 yr
  change by <= 0.5%).

ADDING NEW PHYSICS
Neutrino processes, conductivity channels and heat-capacity terms are kept in
three registries (registry.py). Adding or replacing an ingredient means writing
one function and registering it; no existing file has to be edited. The module
with your function must be imported before the run (for example from your own
script that then calls run_cooling). All three kinds have the same signature:
f(state) -> array of the shape of state.T, with state a state.LocalState.

    from emissivity import NEUTRINO_PROCESSES            # erg cm^-3 s^-1
    from conductivity import CONDUCTIVITY_CHANNELS       # erg cm^-1 s^-1 K^-1
    from heat_capacity import HEAT_CAPACITY_TERMS        # erg cm^-3 K^-1

    @NEUTRINO_PROCESSES.register("my_process")
    def my_process(state):
        return 1e20 * (state.T * 1e-9) ** 8 * (state.rho > 2.8e14)

    @CONDUCTIVITY_CHANNELS.register("my_channel")
    def my_channel(state):
        return np.where(state.crust, 1e18 * state.T, 0.0)

    @HEAT_CAPACITY_TERMS.register("my_term")
    def my_term(state):
        return 1e18 * state.T * (state.xh == 0)

The state has: T (local temperature, K); rho, ye, ymu, yn, yp, xh, aa, zz, ephi,
tc_crust, tc_n, tc_p, kf_n, kf_p, kf_e, kf_mu, mstar_n, mstar_p, transport_mstar_n,
transport_mstar_p (the background fields); nb (baryon number density) and the
number densities n_n, n_p, n_e, n_mu; the masks crust (nuclei present) and core;
direct_urca_allowed; tau_crust, tau_n, tau_p (= T/Tc, 10 without superfluidity);
superfluid. emissivity.process_factors(state) gives the pairing-reduction factors.
Replace an existing ingredient: register with the same name and replace=True (the
position in the sum is kept). Switch one off: REGISTRY.unregister(name). Registered
functions are summed in registration order; removing and registering again moves a
function to the end, which changes only round-off. The default entries are
listed by REGISTRY.names(). The energy-balance check (below) includes any new
ingredient automatically.

ENERGY-BALANCE CHECK
    python cooling_driver.py --check-energy
checks dE/dt = -(L_nu + L_gamma), with E = sum V e^Phi int(Cv dT), the neutrino
luminosity sum V e^(2 Phi) Q and the photon luminosity of the envelope, all
computed from the registered microphysics and the profiles of every accepted
step. It prints, for each decade of time, the closure error: the fraction of
the energy lost in the window that the luminosities do not account for, and
writes energy_budget.csv (E, L_nu, L_gamma, radiated energies, residual). The
closure of a window integrates the luminosities within it. Options:
--energy-start-years (default 1: the unresolved first ~0.2 yr never closes) and
--energy-tolerance (default 0.02); both must be positive and finite. A positive
error means that the star loses less energy than the luminosity integral says
(time-step lag). For the three default BSk26 stars (1.0, 1.4 and 1.8 Msun) the largest window errors
are 1.4%, 3.0% and 2.0% with the updated rates and 0.29%, 0.92% and 1.5% with the standard ones. The
3.0% of the 1.4 Msun star, above the default tolerance of 2%, is in the 0.1-1 kyr window, where the
delayed thermal relaxation of the updated rates makes L_inf drop steeply. It is a time-step error: halving
max_change, grow_below_change and target_change (CoolingConfig) gives 0.99% and quartering them 0.36%,
and the age at which L_inf falls to a given value then moves by about 2% in that window (the
luminosity itself by 7% at 400 yr, 0.3% at 300 yr and 0.5% at 1 kyr; below 0.05% at the other ages checked).
Use the check after adding physics or changing tolerances. The functions are in
energy_balance.py (energy_budget, energy_closure).

REFERENCES AND FORMULAS
The formulas of the model and their sources. References marked (believed) were not recorded in the
original implementation and are given from memory of the literature: check them against the papers
before quoting them. The paper of Bottaro, Caputo and Fiorillo and the expressions derived for this
version are cited exactly, and the places where the implementation differs from a source are stated.

Structure of the star
- Tolman-Oppenheimer-Volkoff equations in x = ln(P/Pc) (tov.py): Tolman (1939), Phys. Rev. 55, 364;
  Oppenheimer and Volkoff (1939), Phys. Rev. 55, 374. exp(Phi) is matched to the Schwarzschild exterior.
- Stable branch of M(Pc), between the minimum and the maximum mass: the turning-point criterion
  (Harrison, Thorne, Wakano and Wheeler 1965) (believed).
- EOS tables (eos_table.py): CompOSE format (compose.obspm.fr; Typel et al. 2015, 2022) (believed).
  BSk22 and BSk26 unified EOS: Pearson et al. (2018), MNRAS 481, 2994 and the Skyrme functionals of
  Goriely, Chamel and Pearson (2013), Phys. Rev. C 88, 024308 (both believed); see the README of each
  folder of eos/.

Thermal evolution
- Heat transport and energy balance in general relativity with the redshifted temperature
  T_inf = exp(Phi) T: Thorne (1977), MNRAS 178, 501 (believed); L_inf = 4 pi R^2 sigma Ts^4 exp(2 Phi).
- Time integration: Crank and Nicolson (1947), Proc. Camb. Phil. Soc. 43, 50; the predictor-corrector
  with midpoint coefficients is a standard second-order scheme for nonlinear problems.
- Envelopes: Gudmundsson-type relation Ts = 1e6 g14^0.25 (Tb8 / 1.288)^(1/1.82) (Gudmundsson, Pethick and
  Epstein 1983, ApJ 272, 286, believed); Fe and accreted fits of Potekhin et al. (2003), ApJ 594, 404.

Heat capacity (heat_capacity.py)
- Electrons and muons: degenerate relativistic gas, C = k_B^2 T E_F k_F / (3 (hbar c)^2), E_F = sqrt((hbar c
  k_F)^2 + m^2 c^4) (standard: Landau and Lifshitz, Statistical Physics; Shapiro and Teukolsky 1983, believed).
  The electron term keeps a fit with rounded constants of the same expression (agreement 0.3%).
- Nucleons: Landau quasiparticles, C = k_B^2 T m* k_F / (3 hbar^2) with the Landau effective mass m* = k_F / v_F
  (Landau 1957; Baym and Pethick 1991, Landau Fermi-Liquid Theory, believed). The tests check the lepton and
  nucleon expressions against a numerical Fermi-Dirac gas.
- Superfluid reduction factors and the dimensionless gap fits: Levenfish and Yakovlev (1994), Astron. Lett. 20,
  43 and Yakovlev, Levenfish and Shibanov (1999), Physics-Uspekhi 42, 737 (both believed).
- Ions (gas, liquid and Debye-solid fits): source not recorded in the original implementation.

Neutrino emissivities (emissivity.py, emissivity_factors.py)
- Modified Urca and nucleon bremsstrahlung, standard rates: Friman and Maxwell (1979), ApJ 232, 541, in the
  parametrization of Yakovlev, Kaminker, Gnedin and Haensel (2001), Phys. Rep. 354, 1; superfluid reduction:
  Yakovlev and Levenfish (1995), A&A 297, 717 and Levenfish and Yakovlev (1994).
- Direct Urca: Lattimer, Pethick, Prakash and Haensel (1991), Phys. Rev. Lett. 66, 2701; emissivity prefactor
  Q = (457 pi / 10080) G^2 cos^2(theta_C) (1 + 3 g_A^2) m*_n m*_p mu_l (k_B T)^6 as in Yakovlev et al. (2001).
  The threshold k_n <= k_p + k_l uses the proton Fermi momentum from the proton density. The muon channel
  (lepton factor mu_mu, threshold k_n <= k_p + k_mu) is the same expression for the muon, verified here
  numerically as described in DIRECT URCA THRESHOLD; it is not taken from a published formula.
- Cooper pair breaking and formation: Yakovlev, Kaminker and Levenfish (1999), A&A 343, 650 (believed) and
  Yakovlev et al. (2001). Plasmon decay and pair annihilation: Itoh et al. (1996), ApJS 102, 411 (believed).
  Electron-ion bremsstrahlung: Haensel, Kaminker and Yakovlev (1996), A&A 314, 328 (believed). Electron-proton
  bremsstrahlung: Kaminker and Haensel (1999), Acta Phys. Pol. B 30, 1125 (believed).
- Updated modified Urca and nucleon bremsstrahlung (emissivity_updated.py): Bottaro, Caputo and Fiorillo (2024),
  JCAP, arXiv:2406.18640: Eqs. 3.3, 3.7-3.9 and 3.11 (nn, pp), 3.13 and 3.15-3.17 (np), 3.22 and 3.28 (modified
  Urca, neutron branch, triangular), 3.39 (proton branch). Parameters of the paper: m_pi = 138 MeV, m_rho =
  769 MeV, C_rho = 1.4 (Bonn potential, Machleidt, Holinde and Elster 1987, Phys. Rep. 149, 1), Landau parameters
  F0' = 0.7 and G0 = 1.1 (Sjoberg 1976; Anastasio and Brown 1977, as cited there), G = G_F cos(theta_C),
  g_A = 1.27. Differences with the paper: the nucleon effective masses enter as overall powers, the pairing
  reduction factors are those of the standard processes, the Friman-Maxwell suppression factors beta (0.56 nn,
  0.66 np, 0.70 pp, 0.68 modified Urca) are kept, and the rates apply in the core only. Not implemented: the
  non-triangular integrals (Eqs. 3.34 and 3.37) and in-medium effects.

Superfluidity (superfluid.py)
- Gap models 1-16 (crust neutron 1S0, proton 1S0, core neutron 3P2) with the parameters of GAP_MODELS: the
  numbering follows the NSCool code (Page) (believed). Dimensionless gap v(T/Tc): Levenfish and Yakovlev (1994)
  and Yakovlev, Levenfish and Shibanov (1999) (believed).

Thermal conductivity (conductivity_core.py, conductivity_crust.py)
- Core electrons and neutrons, the collision frequencies and the factors R, Z, R_C: Shternin and Yakovlev (2007),
  Phys. Rev. D 75, 103004 (electron and muon conduction by the exchange of transverse plasmons) and Baiko,
  Haensel and Yakovlev (2001), A&A 374, 151 (neutrons), with the superfluid factors of Shternin and Yakovlev
  (2008) (all believed; the original implementation does not name them).
- Muon conduction: the adaptation of the electron expression described in MICROPHYSICS (own approximation). The
  coupled electron-muon treatment of Shternin and Yakovlev (2007) is not implemented.
- Crust electrons: Potekhin's electron-conductivity fit (2019 version), COULAN (1999) (Potekhin, Baiko, Haensel and
  Yakovlev 1999, A&A 346, 345, believed) and TAUEESY (2016), as named in the original implementation; Qimp = 0.1.

TESTS
    python -m pip install -r requirements-dev.txt     (pytest)
    make test                                         (python -m pytest tests, ~2 min)
Covered: background validation, composition (charge neutrality, baryon number, muon
fractions and Fermi momenta) and table round trip; the local state; the registries
(adding, replacing and removing physics); the heat capacity (the lepton and Landau
formulas against a numerical Fermi-Dirac gas, their limits, the muon term and the sum of the
registered terms); microphysics regimes (superfluid heat-capacity jump); the core conductivity with muons (the npe
results reproduced exactly, the lepton expression and its velocity factor, the composition, the
response to pairing and the registered channel); envelope slopes against
finite differences and the Fe fit; the order of the theta-method against an exact
solution (1 for backward Euler, 2 for Crank-Nicolson); monotonic cooling with the default start; the energy balance (the default
scheme must close, plain backward Euler and a deliberately inconsistent solver must
be flagged, new physics passes automatically); input parsing; the command line.
For the EOS star: the reader (BSk22 and BSk26) and its refusal of corrupted tables, the
structure equations and their analytic limit, the mass solver and the stable branch,
the mesh, the background conventions, the effective masses, the Direct Urca threshold
against the published value, the surface factor (the heat removed at the base equals
the reported luminosity), the resolution warning and the mesh convergence. For the updated
emissivities: the closed forms against numerical integrals, the values of the paper's figures, the
muon channels and thresholds, the scalings, the registry and the comparison
script. For the direct Urca with muons: the matrix element with explicit spinors, the phase-space
integral with exact kinematics, the thresholds and the doubling where both channels are open.

Validation of the current default: the 1.4 Msun Gudmundsson run with FM MURCA
completed in 414 accepted steps (3 rejected), stopping at 741495.54 yr
on the photon luminosity threshold. No factor-of-20 scaling is applied.
