"""Tests of the envelope and boundary area (decision D-028). Run them with `make test-envelope` when this part of the code changes, or on request.

Code covered: thermal_fluxes.py, thermal_boundary.py.
Every test states a hypothesis, a reference and a criterion, in the style of the validation steps of the GitHub code.
"""

import numpy as np
import pytest
from cooling import surface_luminosity
from thermal_fluxes import thermal_fluxes
from thermal_geometry import thermal_geometry
from thermal_matrix import evaluate_coefficients


G14 = 1.7392597708694928


@pytest.mark.parametrize("model", ["accreted"])
def test_analytic_slope_matches_finite_difference(model):
    """Hypothesis: the analytic slope of the envelope flux is its derivative with respect to T_b.
    Reference: centred finite differences.
    Criterion: relative 1e-8 for T_b between 2e7 and 2e9 K (accreted envelope).
    """
    for Tb in np.geomspace(2e7, 2e9, 25):
        h = 1e-5 * Tb
        up = thermal_fluxes(Tb + h, G14, model).flux_internal
        down = thermal_fluxes(Tb - h, G14, model).flux_internal
        finite_difference = (up - down) / (2 * h / 1e8)  # dF/dTb8
        assert thermal_fluxes(Tb, G14, model).flux_slope == pytest.approx(
            finite_difference, rel=1e-8
        )


def test_fe_fit_uses_zeta_over_three():
    """Hypothesis: the iron envelope follows the fit of Potekhin, Chabrier and Yakovlev (1997), Appendix A.
    Reference: Ts^4 = g14 [(7 zeta)^2.25 + (zeta/3)^1.25] in units of 1e6 K, with zeta = t - 1e-3
        g14^0.25 sqrt(7 t).
    Criterion: relative 1e-13 at T_b = 1e8, 5e8 and 2e9 K.
    """
    for Tb in (1e8, 5e8, 2e9):
        t = Tb / 1e9
        zeta = t - 1e-3 * G14**0.25 * np.sqrt(7 * t)
        ts = 1e6 * (G14 * ((7 * zeta) ** 2.25 + (zeta / 3) ** 1.25)) ** 0.25
        assert thermal_fluxes(Tb, G14, "fe").surface_temperature_K == pytest.approx(ts, rel=1e-13)


def test_gudmundsson_scaling():
    """Hypothesis: T_b = 1.288e8 K gives T_s = 1e6 K for g14 = 1.
    Reference: the Gudmundsson, Pethick and Epstein (1983) relation (README, to be checked).
    Criterion: relative 1e-12.
    """
    ts = thermal_fluxes(1.288e8, 1.0, "gudmundsson").surface_temperature_K
    assert ts == pytest.approx(1e6, rel=1e-12)  # Tb = 1.288e8 K gives 1e6 K for g14 = 1


def test_the_flux_leaving_the_base_is_the_luminosity_that_is_reported(background):
    """Hypothesis: the heat removed through the boundary equals the reported luminosity: energy is
        conserved across the envelope.
    Reference: L_inf = 4 pi R^2 sigma T_s^4 exp(2 Phi), on the default star at uniform T_inf = 5e8
        K.
    Criterion: relative 1e-9.
    """
    b = background
    u = np.full(b.n_cells, 5.0)  # T_inf / 1e8 K
    geometry = thermal_geometry(b)
    boundary = evaluate_coefficients(u, b, geometry=geometry).boundary
    removed = (
        geometry.areas[-1] * (boundary.flux_offset + boundary.flux_slope * u[-1]) * 1e40
    )  # erg/s
    assert removed == pytest.approx(surface_luminosity(u[-1], b, geometry), rel=1e-9)
