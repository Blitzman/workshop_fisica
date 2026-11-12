# Lab notebook — neutron-star cooling

Rules: `../AGENTS.md` applies here unchanged (hypothesis before testing, external validation first, no free constants, every reported number generated from `results/*.json`, physical premise validated with José before code).

Line references `paper:N` point to `paper/marino2024.txt` (text extracted from the PDF José sent; both are kept locally and not committed — the paper is arXiv:2404.05371).

## 0. Status

- 2026-11-12. Folder created. Paper read in full from the PDF. **No code yet: premise and claim list pending José's validation.**

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
