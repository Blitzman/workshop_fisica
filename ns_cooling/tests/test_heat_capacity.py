"""Heat capacity: the muon term and the Landau-quasiparticle nucleons, tested separately and together."""

from contextlib import contextmanager
from dataclasses import replace

import numpy as np
import pytest
from scipy.integrate import quad
from scipy.optimize import brentq

from config import CoolingConfig
from constants import HBAR, HBARC_MEV_FM, K_B, M_NEUTRON, M_PROTON, MEV_TO_ERG, MUON_MASS_MEV
from cooling import run_cooling
from heat_capacity import (
    HEAT_CAPACITY_TERMS,
    degenerate_lepton_heat_capacity,
    electron_heat_capacity,
    heat_capacity,
    heat_capacity_pairing_factor,
    landau_heat_capacity,
    muon_heat_capacity,
    neutron_heat_capacity,
    proton_heat_capacity,
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


# ------------------------------------------------------------------ the two changes, separately
def old_nucleon_formula(state, fraction, mass):
    """Relativistic gas of mass m*: the formula that the Landau-mass expression replaces."""
    x = (6.497229e-14 / mass) * (state.nb * fraction) ** (1.0 / 3.0)
    return 4.5507e11 * mass**2 * x * np.hypot(x, 1.0) * state.T


def old_neutrons(state):
    cv = old_nucleon_formula(state, state.yn, state.mstar_n)
    if state.superfluid:
        cv *= heat_capacity_pairing_factor(state.T, np.where(state.crust, state.tc_crust, 0.0), 1)
        cv *= heat_capacity_pairing_factor(state.T, np.where(state.core, state.tc_n, 0.0), 2)
    return cv


def old_protons(state):
    cv = old_nucleon_formula(state, state.yp, state.mstar_p)
    if state.superfluid:
        cv *= heat_capacity_pairing_factor(state.T, state.tc_p, 1)
    return cv


@contextmanager
def heat_capacity_terms(*, muons, landau):
    """Temporarily select the muon term and the nucleon formula; the registry is restored afterwards."""
    saved = list(HEAT_CAPACITY_TERMS.items())
    try:
        for name, _ in saved:
            HEAT_CAPACITY_TERMS.unregister(name)
        for name, term in saved:
            if name == "muons" and not muons:
                continue
            if not landau and name in ("neutrons", "protons"):
                term = {"neutrons": old_neutrons, "protons": old_protons}[name]
            HEAT_CAPACITY_TERMS.register(name)(term)
        yield
    finally:
        for name in list(HEAT_CAPACITY_TERMS.names()):
            HEAT_CAPACITY_TERMS.unregister(name)
        for name, term in saved:
            HEAT_CAPACITY_TERMS.register(name)(term)


def total(state, **options):
    with heat_capacity_terms(**options):
        return heat_capacity(state)["total"]


def test_the_registry_is_restored(core_state):
    names = list(HEAT_CAPACITY_TERMS.names())
    with heat_capacity_terms(muons=False, landau=False):
        assert "muons" not in HEAT_CAPACITY_TERMS.names()
    assert list(HEAT_CAPACITY_TERMS.names()) == names


def test_the_correction_is_the_removal_of_the_relativistic_factor(core_state):
    """With the same constants the two nucleon formulas differ exactly by sqrt(1 + x^2)."""
    s = core_state
    core = s.core & (s.yn > 0)
    x = (6.497229e-14 / s.mstar_n) * (s.nb * s.yn) ** (1.0 / 3.0)
    ratio = neutron_heat_capacity(s)[core] / old_neutrons(s)[core]
    np.testing.assert_allclose(
        ratio * np.hypot(x, 1.0)[core], 1.0, atol=3e-3
    )  # constants agree to ~0.1%
    assert (
        np.all(ratio < 1.0) and ratio[0] < 0.8
    )  # overestimated throughout the core, by 30% at the centre
    x_p = (6.497229e-14 / s.mstar_p) * (s.nb * s.yp) ** (1.0 / 3.0)
    ratio_p = proton_heat_capacity(s)[core] / old_protons(s)[core]
    np.testing.assert_allclose(ratio_p * np.hypot(x_p, 1.0)[core], 1.0, atol=3e-3)


def test_each_change_is_independent_and_additive(core_state):
    base = total(core_state, muons=False, landau=False)
    muons_only = total(core_state, muons=True, landau=False)
    landau_only = total(core_state, muons=False, landau=True)
    both = total(core_state, muons=True, landau=True)
    np.testing.assert_allclose(
        muons_only - base, muon_heat_capacity(core_state), rtol=1e-12, atol=0
    )
    np.testing.assert_allclose(both, muons_only + (landau_only - base), rtol=1e-12)
    np.testing.assert_allclose(
        both, heat_capacity(core_state)["total"], rtol=1e-14
    )  # the default is both


def test_the_size_of_each_change_at_the_centre(core_state):
    base = total(core_state, muons=False, landau=False)[0]
    muons_only = total(core_state, muons=True, landau=False)[0]
    landau_only = total(core_state, muons=False, landau=True)[0]
    both = total(core_state, muons=True, landau=True)[0]
    assert 1.05 < muons_only / base < 1.12  # muons add 8%
    assert 0.78 < landau_only / base < 0.86  # the Landau formula removes 18%
    assert 0.86 < both / base < 0.94  # -10% together


def test_the_crust_changes_only_through_the_free_neutrons(core_state):
    s = core_state
    base = total(s, muons=False, landau=False)
    both = total(s, muons=True, landau=True)
    outer = s.crust & (s.rho < 3e11)  # no free neutrons above the neutron drip: nothing changes
    np.testing.assert_array_equal(both[outer], base[outer])
    inner = s.crust & ~outer  # the dense neutron gas of the inner crust loses up to ~5%
    change = both[inner] / base[inner] - 1.0
    assert -0.06 < change.min() < -0.02 and change.max() < 0.003


# ------------------------------------------------------------------ the effect on the cooling curve
def luminosity_at(result, t):
    return float(
        np.exp(np.interp(np.log(t), np.log(result.times_years), np.log(result.luminosity_erg_s)))
    )


def test_the_cooling_effect_of_each_change(bsk22_model):
    """More heat capacity cools more slowly: the muons raise L, the Landau correction lowers it."""
    config = replace(CoolingConfig(), final_years=1e3)
    runs = {}
    for key, options in dict(
        base=(False, False), muons=(True, False), landau=(False, True), both=(True, True)
    ).items():
        with heat_capacity_terms(muons=options[0], landau=options[1]):
            runs[key] = run_cooling(bsk22_model.background, config)
    effect = {
        key: luminosity_at(runs[key], 1e3) / luminosity_at(runs["base"], 1e3) - 1.0 for key in runs
    }
    assert 0.003 < effect["muons"] < 0.03  # about +1%
    assert -0.05 < effect["landau"] < -0.01  # about -2.3%
    assert effect["landau"] < effect["both"] < effect["muons"]
    assert effect["both"] == pytest.approx(
        effect["muons"] + effect["landau"], abs=0.005
    )  # nearly additive
