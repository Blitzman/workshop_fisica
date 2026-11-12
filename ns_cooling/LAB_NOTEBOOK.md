# Lab notebook — neutron-star cooling

Rules: `../AGENTS.md` applies here unchanged (hypothesis before testing, external validation first, no free constants, every reported number generated from `results/*.json`, physical premise validated with José before code).

Line references `paper:N` point to `paper/marino2024.txt` (text extracted from the PDF José sent; both are kept locally and not committed — the paper is arXiv:2404.05371).

## 0. Status

- 2026-11-12. Folder created. Paper read in full from the PDF. Premise and claim list pending José's validation. Started anyway (user decision) with steps that only reproduce published numbers.
- 2026-11-12. H1 confirmed: BSk24 structure and direct-Urca threshold reproduce the paper.
- 2026-11-12. H2: SLy4 structure reproduces Douchin & Haensel 2001; no direct Urca in any stable star. Next: the cooling solver, validated against published curves before any sweep.
- 2026-11-12. H3 confirmed: integrated heat capacity and nn-bremsstrahlung luminosity of the isothermal core reproduce the fits of Ofengeim et al. 2017. Next: step 4, the cooling curve with a published envelope relation.
- 2026-11-12. H4: standard slow cooling reproduces Yakovlev et al. 2011 up to 1.6 M☉; 9–11 % hotter at 1.8–2.0 M☉ (open, I1). Next: envelope relation and photon stage, then the BSk24 1.4 M☉ curve of Potekhin & Chabrier 2018.

## 1. The brief (José, 2026-11-11)

A 1D cooling code, no magnetic field, with up-to-date microphysics, to confirm or refute statements in Marino, Dehman, Kovlakas, Rea, Pons & Viganò 2024 (Nat. Astron. 8, 1020) that he considers not rigorous. In his words: a full update of things not updated in years, "much simpler, 1D code, but with real microphysics".

## 2. What the paper does (verified against the PDF)

- **Sources.** PSR J0205+6449 (841 yr), PSR B2334+61 (7700 yr), CXOU J0852−4617 (2500–5000 yr) (paper:155–157). Thermal luminosities, in 10³² erg/s: 1.44 (+0.6/−0.5), 0.78 (+0.30/−0.10), 2.0 ± 0.6, 90 % confidence (Supplementary Table 6, paper:2541).
- **Simulations.** 81 runs of the 2D magneto-thermal code: 3 EoS (SLy4, BSk24, GM1A) × 3 masses (1.4, 1.6, 1.8 M☉) × 9 dipole fields from 1e12 to 7e13 G, no toroidal field (paper:194–213).
- **Thermal treatment.** Heat diffusion solved only in the crust, from an initial 1e10 K; core taken as isothermal; outer boundary through a Ts–Tb relation for an iron envelope; blackbody surface (paper:1352–1379).
- **Superfluidity.** One gap model, "the superfluid gap proposed by" Ho et al. 2015 (PRC 91, 015806), without saying which of that paper's gaps (paper:1314). The hyperon runs use a different gap model (paper:1401).
- **Other microphysics.** Not listed; the text refers to the Potekhin, Pons & Page 2015 review and Potekhin's public codes (paper:1326–1330).
- **Fast-cooling thresholds quoted.** BSk24: direct Urca above 8.25e14 g/cm³, active "when the mass is larger than 1.6 M☉". GM1A: enhanced cooling at all masses from 1.4 M☉, density 5.949e14 g/cm³. SLy4: none, even at 2 M☉ (paper:237, 1422–1430, 1489).
- **Conclusions.** SLy4 ("minimal cooling") cannot reach the three luminosities for any mass or field; classifier probability < 5 % (paper:228–232, 284–288). "Only EoSs allowing a fast-cooling process in the first few kyr" fit. About 75 % of EoS in a simplified nucleonic meta-model (Margueron et al. 2018) lack the proton fraction for fast cooling and are "potentially excluded" (paper:315–322).

## 3. Statements to test (draft for José)

- **S1. "Minimal cooling cannot reach these luminosities."** Tested with one non-dUrca EoS and one gap model. The paper itself says gaps "significantly impact the cooling timescales" and defers the exploration (paper:1288–1312). Open: neutron triplet gap and pair-breaking emission, in-medium enhancement of modified Urca, updated conductivities and heat capacity.
- **S2. "Unless extremely unusual choices are made for the gap models, cooling is typically accelerated"** (paper:1306). Stated without a calculation. Directly testable by a gap sweep.
- **S3. "Iron envelope is the most extreme (coldest) case"** (paper:1470–1477). True as a rule of thumb at young ages; to be checked with current envelope fits, including partially accreted ones.
- **S4. Isothermal core from the start.** For the 841-yr source the crust–core relaxation is not obviously over; a full 1D solve (core + crust) checks whether this biases the young end of the curves.
- **S5. Direct-Urca thresholds** for BSk24 and GM1A (density and mass). Checkable against the EoS papers (Pearson et al. 2018; Gusakov et al. 2014).
- **S6. "About 75 % of EoS excluded."** Not from the simulations. Recomputable as a dUrca-threshold census over a stated EoS set; also depends on S1.
- **S7. Statistics.** A classifier whose classes are the 81 runs can only rank those runs; "< 5 %" for SLy4 is relative to that set. Not a physics check; noted, lower priority.
- Out of scope for now: the X-ray spectral analysis, distances and ages.

## 4. Premise and order of magnitude (to be confirmed, not our result)

The sources sit near 1e32 erg/s at 1e3–1e4 yr. Published minimal-cooling curves at those ages are usually quoted around 1e33–1e34 erg/s. The question is whether updated microphysics without direct Urca closes a gap of one to two orders of magnitude. This figure must be confirmed when we reproduce published curves; it is from memory of the literature, not from a calculation here.

## 5. Proposed plan (after José's answer)

1. **External validation first.** Build stars from tabulated EoS (TOV) and reproduce published numbers: mass–radius and dUrca thresholds from the EoS papers.
2. **Reference cooling curves.** Reproduce published curves from other authors with their ingredients (Page et al. 2004 minimal cooling; the public NSCool code) before adding anything new.
3. **Reproduce the paper's B = 0 limit** for SLy4 and BSk24 (their lowest field, 1e12 G, should be close).
4. **Then the sweeps** for S1–S4, one ingredient at a time, each with its hypothesis written here first.

## 6. Questions for José

1. Which statements does he consider not rigorous? Do they match S1–S7?
2. Which ingredients does he consider outdated (gaps, in-medium modified Urca, conductivities, envelopes, heat capacity)?
3. Which gap from Ho et al. 2015 did the paper use, and which did the hyperon runs use?
4. Start from an existing 1D code of his, or write a new one and validate it against NSCool?
5. Which EoS set for the comparison, and from where (CompOSE tables)?
6. Can we have the 81 curves (or the SLy4 and BSk24 ones) as data, to compare number by number?

## 7. Hypotheses

(none yet)

### H1 · BSk24 star structure and direct-Urca threshold (2026-11-12)

Started before José's answer, at the user's request; this step only reproduces published numbers, so it does not depend on his choices.

- **Source of the EoS.** Potekhin's `bskfit18.f` (Pearson et al. 2018), downloaded from http://www.ioffe.ru/astro/NSG/BSk/ to `data/src/`. Coefficients are parsed from that file, not typed. Formulas ported: `BSKfit`, `BSkEofN`, `CORE18`.
- **Hypothesis.** With these fits and a TOV integration, BSk24 gives the numbers quoted in the paper (paper:1422–1425, 237): maximum mass 2.279 M☉ (within 0.5 %), central density there 2.26e15 g/cm³ (within 3 %), direct-Urca threshold density 8.25e14 g/cm³ (within 3 %), and a threshold mass of about 1.6 M☉ (between 1.55 and 1.65).
- **Direct-Urca condition.** Fermi momenta must close a triangle: p_Fn ≤ p_Fp + p_Fe, with Y_p = Y_e + Y_μ.
- **Negative controls.** (1) Ignoring muons in the triangle (Y_p = Y_e) must move the threshold density by more than 3 %. (2) BSk26 coefficients must miss the maximum mass by more than 1 %.
- **Status: CONFIRMED** (`step1_structure.py`, 6/6 checks, a few seconds; numbers below inserted from `results/step1_structure.json`).
  - Maximum mass 2.279 M☉ (+0.01 % vs the paper); central density there 2.282e+15 g/cm³ (+0.96 %).
  - Direct Urca opens at n_b = 0.4527 fm⁻³, 8.251e+14 g/cm³ (+0.01 %), proton fraction 0.1364; the star with that central density has 1.595 M☉ and R = 12.58 km. R(1.4 M☉) = 12.59 km.
  - Controls: without muons the threshold density moves +35.7 %; BSk26 gives 2.169 M☉.
  - *Reading:* the paper's BSk24 numbers follow from the public fits, so our star and theirs are the same object. Consequence for S5: the paper's three masses put only 1.8 M☉ clearly above the threshold; 1.6 M☉ sits 0.005 M☉ above it, so its direct-Urca core is tiny and its curve should be very sensitive to the exact mass. Still to compare with Pearson et al. 2018's own table (not in hand).

### H2 · SLy4 star structure and absence of direct Urca (2026-11-12)

- **Sources.** EoS table `sly4.dat` (Douchin & Haensel's own table, from http://www.ioffe.ru/astro/NSG/NSEOS/). Reference numbers and core composition from Douchin & Haensel 2001 (A&A 380, 151; arXiv:astro-ph/0111092), Tables 4 and 6 and Sect. 7.6, read from the PDF (`paper/dh2001/`, not committed). Both are parsed, not typed.
- **Hypothesis.** TOV on the table reproduces their Table 6: maximum mass 2.05 M☉, R = 9.99 km (each within 0.5 %), central baryon density 1.21 fm⁻³ and central pressure 1.38e36 dyn/cm² (each within 3 %). Their Table 4 composition gives the electron direct-Urca threshold at 1.35 fm⁻³ and the muon one at 1.44 fm⁻³ (each within 0.02 fm⁻³). Hence the paper's statement S5 for SLy4: no direct Urca in any stable star, including 2 M☉ (paper:1489).
- **Note.** Table 6 prints the central density as 2.86 in units of 1e14 g/cm³; their Table 5 (2.83e15 at 1.2 fm⁻³) shows the unit is 1e15. We compare with 2.86e15.
- **Negative controls.** (1) Dropping the relativistic corrections in TOV must change the maximum mass by more than 10 %. (2) Ignoring muons (x_e = x_p) must move the electron threshold by more than 0.05 fm⁻³.
- **Status: structure CONFIRMED; quoted thresholds REFUTED; the conclusion (no direct Urca) CONFIRMED** (`step2_sly4.py`, 7/9 checks; numbers from `results/step2_sly4.json`). Criteria not changed after the fact.
  - Maximum mass 2.050 M☉ (+0.01 %), R = 9.99 km (+0.02 %), n_c = 1.206 fm⁻³ (-0.32 %), P_c = 1.371e+36 dyn/cm² (-0.63 %). R(1.4 M☉) = 11.71 km, central density 0.538 fm⁻³.
  - Thresholds from their Table 4 with the triangle condition: electron 1.414 fm⁻³ (they quote 1.35), muon 1.496 fm⁻³ (they quote 1.44). Both come out 0.06–0.06 fm⁻³ higher than quoted: 2 FAIL. The table's last rows do not satisfy the condition at the quoted density, so the quoted figures and the table are not consistent with each other at this level; cause unknown (rounding of the table, or a different interpolation).
  - Either way the threshold is above the largest central density any stable SLy4 star reaches (1.206 fm⁻³; a 2 M☉ star has 0.956 fm⁻³). The paper's statement that SLy4 never activates direct Urca holds.
  - Controls: Newtonian structure and the muon-free triangle both break the result (electron threshold 1.147 fm⁻³ without muons).

## 8. Decisions

- 2026-11-12 (user). We write our own 1D code; we do not wait for a code from José. Standing permission to download arXiv papers into `paper/refs/` (not committed). Downloaded so far: astro-ph/0012122 (Yakovlev et al. 2001, neutrino emission), 1612.04672 (Ofengeim et al. 2017, analytic L_ν and C), 1711.07662 (Potekhin & Chabrier 2018, BSk24 cooling), 1507.06186 (Potekhin, Pons & Page 2015 review), astro-ph/9706148 (Potekhin, Chabrier & Yakovlev 1997, envelopes), 1010.1154 (Yakovlev et al. 2011, cooling rates).
- Build order: (3) integrated neutrino luminosity and heat capacity of an isothermal star, checked against published fits; (4) cooling curve of the isothermal star with a published envelope relation; (5) real microphysics (emissivities, effective masses, superfluid factors), each from its source; (6) full heat diffusion (needed for S4); then the sweeps.
- Side note for S5: Potekhin & Chabrier 2018 give the BSk24 direct-Urca threshold mass as 1.596 M☉ and remark that 1.6 M☉ exceeds it by 0.004 M☉ (refs/1711.07662, Sect. on Fig. 4). Our H1 value agrees. This is the independent source H1 was missing.

### H3 · Integrated neutrino luminosity and heat capacity of the isothermal core (2026-11-12)

- **Reference.** Ofengeim, Fortin, Haensel, Yakovlev & Zdunik 2017 (arXiv:1612.04672): Eqs. (2)–(4), (6)–(11) define the integrands and integrals; Eqs. (19)–(20) and Table 2 give analytic fits calibrated on eight EoS, among them SLy (Douchin & Haensel), for M from 1 M☉ to the maximum mass at T̃ = 1e9 K. Quoted accuracy: heat capacity rms 5 %, maximum 12 %; nn-bremsstrahlung rms 14 %, maximum 54 %; modified Urca rms 25 %, maximum 66 %.
- **What we compute.** With our SLy4 star (H2) and its Table 4 composition: C_core (non-superfluid: n, p with m* = 0.7 m; e, μ relativistic, their Eq. 6), L_nn (their Eq. 4) and L_MU (their Eq. 2). The integrals carry the redshift exp(2Φ) and the proper volume (their Eqs. 8–9), with T = T̃ exp(−Φ).
- **Hypothesis.** For M = 1.0, 1.2, …, 2.0 M☉: (a) C_core within 12 % of the fit at every mass; (b) L_nn within 54 % at every mass, rms at most 20 %; (c) their constant c0 = 1.12e20 follows from Eq. (6) with m* = 0.7 m_n at n0 within 1 %.
- **Not a pass/fail item.** L_MU depends on the factor Ω, which the paper only describes as "∼ 1"; we report the ratio with Ω = 1.
- **Universality test (prediction, not validation).** BSk24 is not in their calibration set; we expect the same tolerances to hold for it.
- **Negative controls.** (1) Dropping the redshift (Φ = 0 in T and in exp(2Φ)) must push L_nn outside 54 % for at least one mass. (2) Integrating over the whole star volume as if it were core at the central density (uniform-density star) must push C_core outside 12 %.
- **Status: CONFIRMED** (`step3_integrals.py`, 5/5 checks, about half a minute; numbers from `results/step3_integrals.json`).
  - SLy4, six masses: heat capacity within 3.9% of the fit (rms 2.3%); nn-bremsstrahlung within 13.0% (rms 8.4%). c0 recomputed: 1.125e+20 erg cm⁻³ K⁻¹.
  - BSk24 (not in their calibration): heat capacity within 9.6%, nn-bremsstrahlung within 15.5%. The universality they claim extends to it.
  - Modified Urca with Ω = 1: our integral is 0.27–0.41 of the fit, falling with mass. So their Ω is about 2.4–3.8 on average, not ∼1: it must contain the other branches of the process. To be settled in step 5 with the full expression of Yakovlev et al. 2001; prediction: agreement within their 66 %.
  - Controls: without redshift and with a uniform-density star the agreement breaks.
  - *Reading:* the structure integrals (redshift, proper volume, core boundary, composition) are right. This is the skeleton every later ingredient plugs into.

### H4 · Standard slow cooling: modified Urca and the interior temperature at 330 yr (2026-11-12)

- **Emissivities.** Yakovlev, Kaminker, Gnedin & Haensel 2001 (arXiv:astro-ph/0012122): modified Urca, neutron branch Eq. (140) and proton branch Eq. (142), with α_n = 1.13, β_n = 0.68; muon branches by the rule after Eq. (142) (p_Fe → p_Fμ, extra factor (n_μ/n_e)^(1/3)); nucleon bremsstrahlung Eqs. (165)–(167) with their α, β and three neutrino flavours. Read from the PDF pages 70–71 and 79. Check on the transcription: Eq. (165) with m* = 0.7 m gives the prefactor 1.77e19 that Ofengeim et al. quote.
- **Reference.** Yakovlev, Ho, Shternin, Heinke & Potekhin 2011 (arXiv:1010.1154), Eq. (14): the redshifted interior temperature of a thermally relaxed "standard candle" (non-superfluid, modified Urca, m* = 0.7 m_n) is T̃ = 3.45e8 K (1 − x)[1 + 0.12 (R/10 km)²](330 yr/t)^(1/6), x = 2GM/(Rc²). Fitted to about 170 stars from 17 EoS, SLy among them; rms 3 %, maximum 7 % for M < 2.4 M☉.
- **What we compute.** Since L_ν ∝ T̃⁸ and C ∝ T̃, the neutrino-stage solution is exact: T̃(t) = (6 q t)^(−1/6) with q = L_ν/(C T̃⁷) from our integrals (their Eq. 10).
- **Hypothesis.** (a) SLy4, M = 1.0–2.0 M☉: T̃(330 yr) within 7 % of Eq. (14) at every mass. (b) Prediction left open in H3: the full modified-Urca luminosity is within 66 % of the Ofengeim fit at every mass.
- **Note on sensitivity.** T̃ ∝ L^(−1/6): a 7 % band in T̃ is a 50 % band in luminosity. This test is coarse; it catches factor-of-two errors, not 10 % ones.
- **Negative controls.** (1) Bremsstrahlung only (no modified Urca) must leave the 7 % band. (2) No redshift must leave the 7 % band.
- **Status: (a) REFUTED at the two heaviest stars, confirmed below; (b) CONFIRMED** (`step4_standard_candle.py`, 4/5 checks; numbers from `results/step4_standard_candle.json`). Criteria not changed.
  - (a) T̃(330 yr) minus Eq. (14), for M = 1.0, 1.2, 1.4, 1.6, 1.8, 2.0 M☉: -1.9%, +0.3%, +2.7%, +5.6%, +8.6%, +10.7%. Inside 7 % up to 1.6 M☉; outside at 1.8 and 2.0 M☉. We are hotter there, i.e. our neutrino luminosity is lower than theirs by up to a factor 1.84.
  - (b) Full modified Urca over the Ofengeim fit: 0.73, 0.71, 0.68, 0.64, 0.60, 0.57. Inside their 66 %, but always below and falling with mass. Modified Urca is 96%–97% of the total; the muon branches carry 37%–44%.
  - *Reading.* Both references say the same thing: our modified-Urca luminosity is low at high mass. Three things known: (1) Yakovlev et al. 2011 state that the largest error of Eq. (14) below 2.4 M☉, about 7 %, is precisely the 1.9 M☉ SLy star, so part of the excess may be the fit and not us (sign not stated); (2) even for nn-bremsstrahlung, where the formula is unambiguous, our integral over the Ofengeim fit falls with mass (0.96 to 0.87), so the fits themselves drift at high mass; (3) Ofengeim et al. quote a modified-Urca prefactor of 1.75e+21, while Eq. (140) with m* = 0.7 m gives 1.494e+21: a factor 1.17 not explained by the text.
  - **Open (I1).** Whether the high-mass deficit is ours or the fits'. Not tuned. To settle with an independent curve for a named star (Potekhin & Chabrier 2018, BSk24 1.4 M☉) in the next step.
