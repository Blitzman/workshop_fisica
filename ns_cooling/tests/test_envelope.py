"""Tb-Ts envelope relations: analytic flux slopes and the Fe fit."""

import numpy as np
import pytest

from thermal_fluxes import MODELS, thermal_fluxes

G14 = 1.7392597708694928


@pytest.mark.parametrize("model", MODELS)
def test_analytic_slope_matches_finite_difference(model):
    for Tb in np.geomspace(2e7, 2e9, 25):
        h = 1e-5 * Tb
        up = thermal_fluxes(Tb + h, G14, model).flux_internal
        down = thermal_fluxes(Tb - h, G14, model).flux_internal
        finite_difference = (up - down) / (2 * h / 1e8)  # dF/dTb8
        assert thermal_fluxes(Tb, G14, model).flux_slope == pytest.approx(
            finite_difference, rel=1e-8
        )


def test_flux_offset_and_slope_reconstruct_the_flux_at_the_expansion_point():
    for model in MODELS:
        r = thermal_fluxes(3e8, G14, model)
        assert r.flux_offset + r.flux_slope * 3.0 == pytest.approx(r.flux_internal, rel=1e-12)


def test_fe_fit_uses_zeta_over_three():
    """Potekhin et al. (2003): Ts^4 = g14 [(7 zeta)^2.25 + (zeta/3)^1.25] (units of 1e6 K)."""
    for Tb in (1e8, 5e8, 2e9):
        t = Tb / 1e9
        zeta = t - 1e-3 * G14**0.25 * np.sqrt(7 * t)
        ts = 1e6 * (G14 * ((7 * zeta) ** 2.25 + (zeta / 3) ** 1.25)) ** 0.25
        assert thermal_fluxes(Tb, G14, "fe").surface_temperature_K == pytest.approx(ts, rel=1e-13)


def test_gudmundsson_scaling():
    ts = thermal_fluxes(1.288e8, 1.0, "gudmundsson").surface_temperature_K
    assert ts == pytest.approx(1e6, rel=1e-12)  # Tb = 1.288e8 K gives 1e6 K for g14 = 1


def test_invalid_inputs_and_models_are_rejected():
    with pytest.raises(ValueError):
        thermal_fluxes(-1.0, G14, "fe")
    with pytest.raises(ValueError):
        thermal_fluxes(1e8, G14, "unknown")
    with pytest.raises(ValueError, match="zeta"):
        thermal_fluxes(1e3, G14, "fe")
