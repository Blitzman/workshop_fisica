"""Tests of the heat capacity and conduction area (decision D-028). Run them with `make test-heat-conduction` when this part of the code changes, or on request.

Code covered: heat_capacity.py, conductivity_core.py, conductivity_crust.py, effective_mass.py, superfluid.py.
Every test states a hypothesis, a reference and a criterion, in the style of the validation steps of the GitHub code.
"""

import numpy as np
import pytest
from conductivity_core import core_conductivity, lepton_conductivity, muon_velocity
from constants import HBAR, HBARC_MEV_FM, K_B, MEV_TO_ERG, MUON_MASS_MEV, M_NEUTRON
from heat_capacity import degenerate_lepton_heat_capacity, heat_capacity, landau_heat_capacity
from scipy.integrate import quad
from scipy.optimize import brentq
from state import LocalState


ELECTRON_MASS_MEV = 0.51099895


# ------------------------------------------------------------------ the formulas against a numerical Fermi-Dirac gas
def numerical_heat_capacity(dispersion, k_f, T):
    """C = du/dT at fixed density of spin-1/2 fermions (natural units, MeV): the independent reference.

    u = (1/pi^2) int k^2 E(k) f dk, n = k_f^3 / (3 pi^2); the chemical potential is solved at each T.
    """

    def fermi(x):
        return 0.5 * (1.0 - np.tanh(0.5 * x))

    def integral(weight, mu, T):
        value, _ = quad(
            lambda k: k * k * weight(k) * fermi((dispersion(k) - mu) / T),
            0.0,
            2.0 * k_f,
            points=[k_f],
            limit=500,
            epsabs=0.0,
            epsrel=1e-13,
        )
        return value / np.pi**2

    n0 = k_f**3 / (3.0 * np.pi**2)

    def energy_density(T):
        e_f = dispersion(k_f)
        mu = brentq(
            lambda m: integral(lambda k: 1.0, m, T) - n0, e_f - 30 * T, e_f + 30 * T, xtol=1e-13
        )
        return integral(dispersion, mu, T)

    h = 0.05 * T
    return (energy_density(T + h) - energy_density(T - h)) / (2.0 * h)


def to_cgs(c_natural_mev3):
    """C [MeV^3] in natural units -> erg cm^-3 K^-1."""
    return c_natural_mev3 / HBARC_MEV_FM**3 * 1e39 * K_B


@pytest.mark.parametrize("kf_fm, mass", [(1.1, ELECTRON_MASS_MEV)])
def test_the_lepton_formula_matches_a_numerical_fermi_dirac_gas(kf_fm, mass):
    """Hypothesis: the lepton heat capacity is that of a degenerate relativistic gas.
    Reference: a numerical Fermi-Dirac integral at k_B T = 1e-3 E_F.
    Criterion: relative 1e-4 (electron, k_F = 1.1 fm^-1).
    """
    k_f = kf_fm * HBARC_MEV_FM  # MeV
    T_mev = 1e-3 * np.sqrt(k_f**2 + mass**2)  # k_B T = 1e-3 E_F: deeply degenerate
    numeric = numerical_heat_capacity(lambda k: np.sqrt(k * k + mass**2), k_f, T_mev)
    T_kelvin = T_mev * MEV_TO_ERG / K_B
    formula = degenerate_lepton_heat_capacity(kf_fm, mass, T_kelvin)
    assert formula == pytest.approx(to_cgs(numeric), rel=1e-4)


@pytest.mark.parametrize("kf_fm, ratio, mass_g", [(2.12, 0.51, M_NEUTRON)])
def test_the_nucleon_formula_matches_a_numerical_gas_of_landau_quasiparticles(kf_fm, ratio, mass_g):
    """Hypothesis: the nucleon heat capacity is that of Landau quasiparticles with the Landau mass.
    Reference: a numerical integral for a gas with e = k^2 / 2m* at k_B T = 1e-3 E_F.
    Criterion: relative 1e-4 (neutron, k_F = 2.12 fm^-1, m*/m = 0.51).
    """
    mass_mev = ratio * mass_g * 2.99792458e10**2 / MEV_TO_ERG  # m* c^2 [MeV]
    k_f = kf_fm * HBARC_MEV_FM
    T_mev = 1e-3 * k_f**2 / (2.0 * mass_mev)  # k_B T = 1e-3 times the Fermi energy
    numeric = numerical_heat_capacity(lambda k: k * k / (2.0 * mass_mev), k_f, T_mev)
    formula = landau_heat_capacity(kf_fm, ratio, mass_g, T_mev * MEV_TO_ERG / K_B)
    assert formula == pytest.approx(to_cgs(numeric), rel=1e-4)


def test_limits_of_the_lepton_formula():
    """Hypothesis: the lepton heat capacity has the right ultra-relativistic and non-relativistic
        limits, and vanishes at k_F = 0.
    Reference: C = k_B^2 T k_F^2 / (3 hbar c) and C = k_B^2 T m k_F / (3 hbar^2).
    Criterion: relative 1e-3, and exactly zero at k_F = 0.
    """
    T = 1e9
    k = 20.0  # fm^-1: ultra-relativistic for a muon
    ultra = K_B**2 * T * (k * 1e13) ** 2 / (3.0 * HBAR * 2.99792458e10)
    assert degenerate_lepton_heat_capacity(k, MUON_MASS_MEV, T) == pytest.approx(ultra, rel=1e-3)
    k = 0.02  # fm^-1: non-relativistic (k_F c = 4 MeV against 106 MeV)
    muon_g = MUON_MASS_MEV * MEV_TO_ERG / 2.99792458e10**2
    non_relativistic = K_B**2 * T * muon_g * (k * 1e13) / (3.0 * HBAR**2)
    assert degenerate_lepton_heat_capacity(k, MUON_MASS_MEV, T) == pytest.approx(
        non_relativistic, rel=1e-3
    )
    assert degenerate_lepton_heat_capacity(0.0, MUON_MASS_MEV, T) == 0.0


MASSES = (0.6, 0.55)


def conductivity(Y_mu=0.0, Y_e=0.06, T_K=3e8, rho=5e14, **kwargs):
    return core_conductivity(T_K, rho, Y_e, effective_masses=MASSES, Y_mu=Y_mu, **kwargs)


# ------------------------------------------------------------------ the muon channel
def test_the_muon_conductivity_is_the_electron_expression_with_the_muon_density_and_velocity():
    """Hypothesis: the muon conductivity is the electron expression with the muon density and velocity.
        This is our own approximation (README, MICROPHYSICS): no independent reference exists.
    Reference: the function lepton_conductivity applied to the muons and to the electrons.
    Criterion: relative 1e-12 for the conductivities and the collision frequencies.
    """
    Y_e, Y_mu, rho = 0.05, 0.03, 5e14
    r = conductivity(Y_mu=Y_mu, Y_e=Y_e, rho=rho)
    nb = rho / 1.66e15
    n_mu, n_p = nb * Y_mu, nb * (Y_e + Y_mu)
    n_n = nb * (1.0 - Y_e - Y_mu)
    kappa, nu_ll, nu_ln, nu_lp = lepton_conductivity(
        3.0, n_mu, n_n, n_p, *MASSES, 1.0, 1.0, 1.0, 1.0, velocity=muon_velocity(n_mu)
    )
    assert float(r["muons"]) == pytest.approx(kappa, rel=1e-12)
    assert float(r["nu_mm"]) == pytest.approx(nu_ll, rel=1e-12)
    assert float(r["nu_mn"]) == pytest.approx(nu_ln, rel=1e-12)
    assert float(r["nu_mp"]) == pytest.approx(nu_lp, rel=1e-12)
    # and the electrons are the same function with the electron density and velocity 1
    k_e = lepton_conductivity(3.0, nb * Y_e, n_n, n_p, *MASSES, 1.0, 1.0, 1.0, 1.0)[0]
    assert float(r["electrons"]) == pytest.approx(k_e, rel=1e-12)


def test_superfluid_heat_capacity_regimes(background):
    """Hypothesis: superfluidity leaves the heat capacity unchanged above T_c, reduces it strongly deep
        inside the gap, and keeps the jump just below T_c.
    Reference: the reduction factors of Yakovlev, Levenfish and Shibanov (1999) (README, believed).
    Criterion: a ratio of exactly 1 at 1e10 K, below 0.5 where T_c > 3 T, and above 1 where 1 <
        T_c/T < 1.5.
    """

    def ratio(T):
        T = np.full(background.n_cells, T)
        on = heat_capacity(LocalState(T, background.cells, superfluid=True))["total"]
        off = heat_capacity(LocalState(T, background.cells, superfluid=False))["total"]
        return on / off

    c = background.cells
    assert np.all(ratio(1e10) == 1.0)  # hotter than every Tc: identical to the normal state
    hot = ratio(2e8)
    deep = (c.xh == 0) & (c.tc_n > 3 * 2e8)
    assert deep.any() and np.all(hot[deep] < 0.5)
    near = (c.xh == 0) & (c.tc_n > 2e8) & (c.tc_n < 1.5 * 2e8)
    assert near.any() and np.all(hot[near] > 1.0)  # the heat-capacity jump at Tc is kept
