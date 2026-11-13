"""Heat capacity: the lepton and Landau-quasiparticle formulas against a numerical Fermi-Dirac gas, and the registered terms."""


import numpy as np
import pytest
from scipy.integrate import quad
from scipy.optimize import brentq

from constants import HBAR, HBARC_MEV_FM, K_B, M_NEUTRON, M_PROTON, MEV_TO_ERG, MUON_MASS_MEV
from heat_capacity import (
    degenerate_lepton_heat_capacity,
    electron_heat_capacity,
    heat_capacity,
    landau_heat_capacity,
    muon_heat_capacity,
)
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


@pytest.mark.parametrize(
    "kf_fm, mass", [(0.5, MUON_MASS_MEV), (0.9, MUON_MASS_MEV), (1.1, ELECTRON_MASS_MEV)]
)
def test_the_lepton_formula_matches_a_numerical_fermi_dirac_gas(kf_fm, mass):
    k_f = kf_fm * HBARC_MEV_FM  # MeV
    T_mev = 1e-3 * np.sqrt(k_f**2 + mass**2)  # k_B T = 1e-3 E_F: deeply degenerate
    numeric = numerical_heat_capacity(lambda k: np.sqrt(k * k + mass**2), k_f, T_mev)
    T_kelvin = T_mev * MEV_TO_ERG / K_B
    formula = degenerate_lepton_heat_capacity(kf_fm, mass, T_kelvin)
    assert formula == pytest.approx(to_cgs(numeric), rel=1e-4)


@pytest.mark.parametrize(
    "kf_fm, ratio, mass_g", [(2.12, 0.51, M_NEUTRON), (1.2, 0.50, M_PROTON), (0.99, 0.34, M_PROTON)]
)
def test_the_nucleon_formula_matches_a_numerical_gas_of_landau_quasiparticles(kf_fm, ratio, mass_g):
    mass_mev = ratio * mass_g * 2.99792458e10**2 / MEV_TO_ERG  # m* c^2 [MeV]
    k_f = kf_fm * HBARC_MEV_FM
    T_mev = 1e-3 * k_f**2 / (2.0 * mass_mev)  # k_B T = 1e-3 times the Fermi energy
    numeric = numerical_heat_capacity(lambda k: k * k / (2.0 * mass_mev), k_f, T_mev)
    formula = landau_heat_capacity(kf_fm, ratio, mass_g, T_mev * MEV_TO_ERG / K_B)
    assert formula == pytest.approx(to_cgs(numeric), rel=1e-4)


def test_limits_of_the_lepton_formula():
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


def test_the_nucleon_heat_capacity_is_proportional_to_the_landau_mass():
    c1 = landau_heat_capacity(2.0, 0.5, M_NEUTRON, 1e9)
    assert landau_heat_capacity(2.0, 1.0, M_NEUTRON, 1e9) == pytest.approx(2.0 * c1, rel=1e-14)
    assert landau_heat_capacity(1.0, 0.5, M_NEUTRON, 1e9) == pytest.approx(0.5 * c1, rel=1e-14)


# ------------------------------------------------------------------ the muon term in the star
@pytest.fixture(scope="module")
def core_state(bsk22_model):
    return LocalState(1e9, bsk22_model.background.cells, superfluid=False)


def test_the_muons_exist_where_the_background_has_them(core_state):
    cv = muon_heat_capacity(core_state)
    assert np.all(cv[core_state.crust] == 0.0) and np.all(cv[core_state.ymu == 0] == 0.0)
    muons = core_state.ymu > 0
    assert muons.any() and np.all(cv[muons] > 0)
    np.testing.assert_allclose(cv, heat_capacity(core_state)["muons"], rtol=0)


def test_the_lepton_formula_agrees_with_the_electron_term(core_state):
    core = core_state.core
    generic = degenerate_lepton_heat_capacity(core_state.kf_e, ELECTRON_MASS_MEV, core_state.T)
    np.testing.assert_allclose(generic[core], electron_heat_capacity(core_state)[core], rtol=5e-3)


def test_the_total_is_the_sum_of_the_terms(core_state):
    h = heat_capacity(core_state)
    assert set(h) == {"electrons", "muons", "neutrons", "protons", "ions", "total"}
    parts = sum(v for k, v in h.items() if k != "total")
    np.testing.assert_allclose(h["total"], parts, rtol=1e-14)


def test_the_centre_has_muons_comparable_to_the_electrons(core_state):
    h = heat_capacity(core_state)
    assert (
        0.5 < h["muons"][0] / h["electrons"][0] < 1.0
    )  # kf_mu < kf_e but E_F(mu) ~ 2 E_F(e) / ...
    assert 0.03 < h["muons"][0] / h["total"][0] < 0.15


