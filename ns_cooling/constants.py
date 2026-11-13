"""Physical constants and unit conversions (CGS): one source for the whole code."""

G = 6.67430e-8  # gravitational constant [cm^3 g^-1 s^-2]
C = 2.99792458e10  # speed of light [cm s^-1]
M_SUN = 1.98841e33  # solar mass [g]: nominal GM_sun / G (decision D-007)
M_U = 1.66053906660e-24  # atomic mass unit [g]
M_NEUTRON = 1.67492749804e-24  # [g]
M_PROTON = 1.67262192369e-24  # [g]
K_B = 1.380649e-16  # Boltzmann constant [erg K^-1]
HBAR = 1.054571817e-27  # reduced Planck constant [erg s]
HBARC_MEV_FM = 197.3269804  # hbar c [MeV fm]
MEV_TO_ERG = 1.602176634e-6
MUON_MASS_MEV = 105.6583755  # m_mu c^2 [MeV]
ELECTRON_MASS_MEV = 0.51099895  # m_e c^2 [MeV]
NEUTRON_MASS_MEV = 939.56542  # m_n c^2 [MeV]
PROTON_MASS_MEV = 938.27209  # m_p c^2 [MeV]
FERMI_COUPLING_MEV = 1.1663788e-11  # G_F / (hbar c)^3 [MeV^-2]
MEV_TO_G = 1.78266192e-27  # 1 MeV/c^2 in g
MEV_FM3_TO_DYN_CM2 = 1.602176634e33  # 1 MeV fm^-3 in dyn cm^-2
FM3_PER_CM3 = 1e39  # (cm^-3) / (fm^-3)

RHO_TO_GEOM = G / C**2  # g cm^-3 -> cm^-2
PRESSURE_TO_GEOM = G / C**4  # dyn cm^-2 -> cm^-2

# The microphysics derives number densities as n = rho / 1.66e-24 (see background.py).
# A background built from an equation of state sets rho = 1.66e-24 g * n_b, so the Fermi
# momenta of the microphysics reproduce the tabulated baryon density exactly.
# This rho is a baryon REST-MASS density (n_b times a fixed mass unit); it is not the
# mass-energy density e/c^2 of the TOV equations, which differs from it by about 1 %.
MICROPHYSICS_MASS_UNIT = 1.66e-24  # [g]

# Reference baryon density n0 of the Yakovlev et al. (2001) emissivity formulas (decision D-014),
# and the same density as a rest-mass density in the convention above (2.656e14 g cm^-3).
N0_FM3 = 0.16  # [fm^-3]
RHO0 = N0_FM3 * FM3_PER_CM3 * MICROPHYSICS_MASS_UNIT  # [g cm^-3]
