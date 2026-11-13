"""Core conductivity with muons: the lepton expression adapted to the muons, and the composition."""

from contextlib import contextmanager
from dataclasses import replace

import numpy as np
import pytest

from conductivity import CONDUCTIVITY_CHANNELS, core_muon_channel, thermal_transport
from conductivity_core import core_conductivity, lepton_conductivity, muon_velocity
from config import CoolingConfig
from constants import HBARC_MEV_FM, MUON_MASS_MEV
from cooling import run_cooling
from state import LocalState

# Results of the npe expressions (no muons) before the muon channel existed.
NPE_REFERENCE = [
    (dict(T_K=1e8, rho=3e14, Y_e=0.05, effective_masses=(0.7, 0.65)),
     1.1931100575530629e24, 4.9652496405565764e23, 1.6896350216087205e24, 12034831661216.97),
    (dict(T_K=1e9, rho=6e14, Y_e=0.08, effective_masses=(0.55, 0.5)),
     4.2626439711037473e23, 2.5728033832737348e23, 6.835447354377482e23, 2334714903976443.0),
    (dict(T_K=3e8, rho=5e14, Y_e=0.07, tau_n=0.5, tau_p=0.3, effective_masses=(0.6, 0.55)),
     5.743049502146786e23, 1.3230468419907491e23, 7.066096344137536e23, 832175007652651.5),
    (dict(T_K=2e8, rho=4e14, Y_e=0.06),
     1.0307872636351306e24, 2.9724936168063277e23, 1.3280366253157633e24, 54005083419468.68),
]  # fmt: skip
MASSES = (0.6, 0.55)


def conductivity(Y_mu=0.0, Y_e=0.06, T_K=3e8, rho=5e14, **kwargs):
    return core_conductivity(T_K, rho, Y_e, effective_masses=MASSES, Y_mu=Y_mu, **kwargs)


# ------------------------------------------------------------------ npe matter is unchanged
@pytest.mark.parametrize("case", NPE_REFERENCE)
def test_without_muons_the_npe_results_are_reproduced(case):
    arguments, electrons, neutrons, total, nu_ee = case
    for extra in ({}, dict(Y_mu=0.0)):
        r = core_conductivity(**arguments, **extra)
        assert float(r["electrons"]) == pytest.approx(electrons, rel=1e-13)
        assert float(r["neutrons"]) == pytest.approx(neutrons, rel=1e-13)
        assert float(r["total"]) == pytest.approx(total, rel=1e-13)
        assert float(r["nu_ee"]) == pytest.approx(nu_ee, rel=1e-13)
        assert float(r["muons"]) == 0.0


# ------------------------------------------------------------------ the muon channel
def test_the_muon_conductivity_is_the_electron_expression_with_the_muon_density_and_velocity():
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


def test_the_finite_mass_enters_as_the_velocity_factor():
    args = (3.0, 0.01, 0.15, 0.02, *MASSES, 1.0, 1.0, 1.0, 1.0)
    ultra = lepton_conductivity(*args)[0]
    v = muon_velocity(0.01)
    assert lepton_conductivity(*args, velocity=v)[0] == pytest.approx(v * ultra, rel=1e-14)
    p = HBARC_MEV_FM * (3 * np.pi**2 * 0.01) ** (1 / 3)
    assert v == pytest.approx(p / np.hypot(p, MUON_MASS_MEV), rel=1e-14) and 0.5 < v < 0.9
    assert muon_velocity(10.0) > 0.99  # ultra-relativistic at very high density
    assert muon_velocity(1e-5) < 0.2  # non-relativistic near the onset


def test_muons_conduct_only_where_there_are_muons():
    mixed = core_conductivity(
        np.full(3, 3e8), np.full(3, 5e14), np.full(3, 0.05), effective_masses=(0.6, 0.55),
        Y_mu=np.array([0.0, 0.01, 0.03]),
    )  # fmt: skip
    assert mixed["muons"][0] == 0.0 and np.all(mixed["muons"][1:] > 0)
    assert mixed["muons"][1] < mixed["muons"][2]  # more muons, more conduction
    assert np.all(np.isfinite(mixed["total"]))
    np.testing.assert_allclose(
        mixed["total"], mixed["electrons"] + mixed["muons"] + mixed["neutrons"], rtol=1e-14
    )


def test_the_muon_conductivity_vanishes_continuously_at_the_onset():
    values = [float(conductivity(Y_mu=y)["muons"]) for y in (1e-8, 1e-6, 1e-4, 1e-2)]
    assert np.all(np.diff(values) > 0) and values[0] / float(conductivity()["electrons"]) < 1e-3


def test_the_pairing_acts_on_the_muons_as_on_the_electrons():
    """Proton pairing removes the lepton-proton collisions and raises both; neutron pairing lowers
    the screening (the factor Z) and the collision frequencies grow, which lowers both."""
    normal = conductivity(Y_mu=0.03)
    for kwargs, direction in (
        (dict(tau_p=0.3), 1.0),
        (dict(tau_n=0.5), -1.0),
        (dict(tau_n=0.5, tau_p=0.3), -1.0),
    ):
        paired = conductivity(Y_mu=0.03, **kwargs)
        change_e = float(paired["electrons"] / normal["electrons"]) - 1.0
        change_mu = float(paired["muons"] / normal["muons"]) - 1.0
        assert change_e * direction > 0 and change_mu * direction > 0, kwargs
        assert change_mu == pytest.approx(
            change_e, abs=0.1
        ), kwargs  # the same response, within 10 points
    assert float(conductivity(Y_mu=0.03, tau_p=0.3)["nu_mp"]) < 0.05 * float(
        normal["nu_mp"]
    )  # pairing suppresses it


# ------------------------------------------------------------------ the composition
def test_the_neutron_density_follows_the_proton_fraction_with_muons():
    base, with_mu = conductivity(), conductivity(Y_mu=0.03)
    assert float(with_mu["neutrons"] / base["neutrons"]) == pytest.approx(
        (1 - 0.06 - 0.03) / (1 - 0.06), rel=1e-13
    )


def test_the_proton_density_and_the_neutron_density_enter_the_lepton_screening():
    """With muons, n_p = n_e + n_mu is larger and n_n smaller than in npe matter at the same Y_e."""
    base, with_mu = conductivity(), conductivity(Y_mu=0.03)
    assert float(with_mu["electrons"]) != float(base["electrons"])


@pytest.mark.parametrize(
    "bad", [dict(Y_mu=-0.01), dict(Y_mu=0.95), dict(Y_mu=np.nan), dict(Y_e=0.5, Y_mu=0.5)]
)
def test_invalid_compositions_are_refused(bad):
    with pytest.raises(ValueError, match="Y_mu >= 0 and Y_e \\+ Y_mu < 1"):
        conductivity(**bad)


# ------------------------------------------------------------------ the registry and the star
def test_the_muon_channel_is_registered_between_the_electrons_and_the_neutrons():
    assert CONDUCTIVITY_CHANNELS.names()[:3] == ("core_electrons", "core_muons", "core_neutrons")


def test_in_the_star_the_muon_channel_is_zero_outside_the_muon_region(bsk22_model):
    cells = bsk22_model.background.cells
    state = LocalState(5e8, cells)
    kappa = thermal_transport(state)
    muons = (cells.ymu > 0) & (cells.xh == 0)
    assert muons.any() and np.all(kappa["kappa_core_muons"][muons] > 0)
    assert np.all(kappa["kappa_core_muons"][~muons] == 0.0)
    total = sum(v for k, v in kappa.items() if k.startswith("kappa_") and k != "kappa_total")
    np.testing.assert_allclose(kappa["kappa_total"], total, rtol=1e-14)
    ratio = kappa["kappa_core_muons"][0] / kappa["kappa_core_electrons"][0]
    assert 0.4 < ratio < 0.9  # about 0.63 at the centre of BSk22 1.4 Msun
    assert (
        kappa["kappa_core_muons"][muons][-1] < 0.02 * kappa["kappa_core_electrons"][muons][-1]
    )  # onset


@contextmanager
def without_the_muon_channel():
    saved = list(CONDUCTIVITY_CHANNELS.items())
    try:
        CONDUCTIVITY_CHANNELS.unregister("core_muons")
        yield
    finally:
        for name in list(CONDUCTIVITY_CHANNELS.names()):
            CONDUCTIVITY_CHANNELS.unregister(name)
        for name, channel in saved:
            CONDUCTIVITY_CHANNELS.register(name)(channel)


def test_the_cooling_curve_hardly_changes_with_the_muon_conductivity(bsk22_model):
    config = replace(CoolingConfig(), final_years=300.0)
    with_muons = run_cooling(bsk22_model.background, config)
    with without_the_muon_channel():
        assert "core_muons" not in CONDUCTIVITY_CHANNELS.names()
        without = run_cooling(bsk22_model.background, config)
    assert CONDUCTIVITY_CHANNELS.names()[1] == "core_muons" and core_muon_channel is not None

    def luminosity(run, t):
        return np.exp(np.interp(np.log(t), np.log(run.times_years), np.log(run.luminosity_erg_s)))

    for t in (10.0, 100.0, 300.0):  # measured: below 1% (-0.9% at 30 yr, -0.2% at 100 yr)
        assert abs(luminosity(with_muons, t) / luminosity(without, t) - 1.0) < 0.02
