"""Direct Urca threshold from Fermi momenta, effective-mass plumbing and surface plumbing."""

from dataclasses import replace

import numpy as np
import pytest

import tov
from constants import RHO0
from cooling import surface_luminosity
from emissivity import direct_urca_emissivity, neutrino_emissivity
from state import LocalState
from synthetic import points_from_table
from thermal_boundary import surface_boundary
from thermal_geometry import thermal_geometry
from thermal_matrix import evaluate_coefficients


def npe_points(ye, rho=1e15):
    """Core-only npe points (Yp = Ye, Yn = 1 - Ye) with the given electron fractions."""
    ye = np.asarray(ye, dtype=float)
    n = len(ye)
    zeros = np.zeros(n)
    table = dict(
        radius_km=np.arange(n, dtype=float),
        rho=np.full(n, rho),
        ye=ye,
        yn=1 - ye,
        yp=ye,
        xh=zeros,
        aa=np.ones(n),
        zz=ye,
        ephi=np.full(n, 0.8),
        elam=np.full(n, 1.2),
        tc_crust=zeros,
        tc_n=zeros,
        tc_p=zeros,
    )
    return points_from_table(table)


# ------------------------------------------------------------------ Direct Urca threshold
YE_SCAN = np.array([0.10, 0.1105, 0.11110, 0.11112, 0.1112, 0.12, 0.2])  # 1/9 = 0.111111...


def test_the_threshold_does_not_depend_on_the_density():
    for rho in (3e14, 1e15, 2e15):
        state = LocalState(1e9, npe_points(YE_SCAN, rho=rho))
        assert list(state.direct_urca_allowed) == [False, False, False, True, True, True, True]


def test_the_proton_fermi_momentum():
    npe = LocalState(1e9, npe_points(YE_SCAN))
    np.testing.assert_array_equal(npe.kf_p, npe.kf_e)  # n_p = n_e without muons


def test_direct_urca_emissivity_switches_on_at_the_threshold():
    faces = npe_points(YE_SCAN)
    state = LocalState(1e9, faces, superfluid=False)
    allowed = state.direct_urca_allowed
    q = direct_urca_emissivity(state)
    assert np.all(q[~allowed] == 0) and np.all(q[allowed] > 0)
    q2 = direct_urca_emissivity(LocalState(2e9, faces, superfluid=False))
    np.testing.assert_allclose(q2[allowed] / q[allowed], 64.0, rtol=1e-12)  # T^6
    expected = (
        4e27
        * state.mstar_n
        * state.mstar_p
        * (state.T * 1e-9) ** 6
        * (state.ye * state.rho / RHO0) ** (1 / 3)
    )
    np.testing.assert_allclose(q[allowed], expected[allowed], rtol=1e-6)
    np.testing.assert_array_equal(
        neutrino_emissivity(state)["direct_urca"], q
    )  # it is the registered process


def test_the_threshold_with_muons_uses_the_proton_momentum(bsk22_model):
    cells = bsk22_model.background.cells
    state = LocalState(1e9, cells)
    nb = cells.rho / 1.66e-24

    def kf(y):
        return (3 * np.pi**2 * nb * y * 1e-39) ** (1 / 3)

    np.testing.assert_allclose(
        state.kf_p, kf(cells.yp), rtol=1e-6
    )  # kf_e (yp/ye)^(1/3) with muons present
    assert (
        state.kf_p[0] > 1.1 * state.kf_e[0]
    )  # the protons exceed the electrons: muons carry charge too
    margin = kf(cells.yp) + kf(cells.ye) - kf(cells.yn)
    expected = (cells.xh == 0) & (margin > 0)
    clear = np.abs(margin) > 1e-5
    np.testing.assert_array_equal(state.direct_urca_allowed[clear], expected[clear])


def pressure_where(eos_table, margin):
    """Pressure (log-linear between core rows) where the margin of a core row array turns positive."""
    core = eos_table.phase == 0
    p, m = np.log(eos_table.pressure[core]), margin[core]
    j = int(np.argmax(m > 0))
    return float(np.exp(p[j - 1] + (-m[j - 1] / (m[j] - m[j - 1])) * (p[j] - p[j - 1])))


def du_threshold_pressure(eos_table):
    """kf_n <= kf_p + kf_e in the table (kf ~ (n Y)^(1/3) at equal n), muons included in Yp."""
    t = eos_table
    return pressure_where(t, t.yp ** (1 / 3) + t.ye ** (1 / 3) - t.yn ** (1 / 3))


def test_the_threshold_mass_matches_the_published_value(eos_table):
    """The CompOSE description of BSk22 quotes M_DU,e = 1.151 Msun for the electron Direct Urca."""
    mass = tov.integrate_star(eos_table, du_threshold_pressure(eos_table)).mass_solar
    assert mass == pytest.approx(
        1.151, rel=1e-2
    )  # 1.148 Msun: the muons of the table are accounted for
    naive = pressure_where(
        eos_table, eos_table.yp - 1 / 9
    )  # the npe rule Yp = 1/9 ignores the muons
    assert (
        tov.integrate_star(eos_table, naive).mass_solar < 0.95
    )  # 0.88 Msun: wrong for this matter


def test_direct_urca_region_of_the_eos_star(bsk22_model, eos_table):
    cells = bsk22_model.background.cells
    allowed = LocalState(1e9, cells).direct_urca_allowed
    rho_threshold = float(eos_table.at_pressure(du_threshold_pressure(eos_table)).rho_baryon)
    far = (
        np.abs(cells.rho / rho_threshold - 1) > 0.02
    )  # cells whose centre is not within 2% of the threshold
    np.testing.assert_array_equal(
        allowed[far], ((cells.xh == 0) & (cells.rho > rho_threshold))[far]
    )
    k = int(allowed.sum())
    assert k >= 3 and np.all(allowed[:k]) and not allowed[k:].any()  # a central region
    assert cells.rho.max() < 10.18e14  # a fixed-density threshold would switch it off everywhere
    assert LocalState(1e9, bsk22_model.background.faces).direct_urca_allowed[0]


def test_the_default_star_never_reaches_the_threshold(background):
    """BSk26, 1.4 Msun: neither the electron nor the muon channel is open anywhere."""
    for points in (background.faces, background.cells):
        state = LocalState(1e9, points)
        assert not state.direct_urca_allowed.any() and not state.direct_urca_muon_allowed.any()


# ------------------------------------------------------------------ surface
def test_the_luminosity_factor_scales_the_flux_and_not_the_temperature():
    a = surface_boundary(2.0, 1.5, 0.8, "accreted")
    b = surface_boundary(2.0, 1.5, 0.8, "accreted", luminosity_factor=1.07)
    assert b.flux_offset == pytest.approx(1.07 * a.flux_offset, rel=1e-14)
    assert b.flux_slope == pytest.approx(1.07 * a.flux_slope, rel=1e-14)
    assert b.surface_temperature_K == a.surface_temperature_K
    assert surface_boundary(2.0, 1.5, 0.8, "accreted", luminosity_factor=1.0) == a


@pytest.mark.parametrize("which", ["table", "eos"])
def test_the_flux_leaving_the_base_is_the_luminosity_that_is_reported(
    which, background, bsk22_model
):
    """Energy conservation across the envelope: the heat removed by the thermal matrix at the
    base equals L_inf, with and without a surface above the last face."""
    b = background if which == "table" else bsk22_model.background
    u = np.full(b.n_cells, 5.0)  # T_inf / 1e8 K
    geometry = thermal_geometry(b)
    boundary = evaluate_coefficients(u, b, geometry=geometry).boundary
    removed = (
        geometry.areas[-1] * (boundary.flux_offset + boundary.flux_slope * u[-1]) * 1e40
    )  # erg/s
    assert removed == pytest.approx(surface_luminosity(u[-1], b, geometry), rel=1e-9)


def test_the_surface_factor_enters_luminosity_and_boundary(bsk22_model):
    b = bsk22_model.background
    plain = replace(b, surface_radius_km=None, surface_ephi=None)
    u = np.full(b.n_cells, 5.0)
    geometry = thermal_geometry(b)
    assert surface_luminosity(u[-1], b, geometry) / surface_luminosity(
        u[-1], plain, geometry
    ) == pytest.approx(b.surface_factor, rel=1e-13)
    c1, c0 = (
        evaluate_coefficients(u, b, geometry=geometry).boundary,
        evaluate_coefficients(u, plain, geometry=geometry).boundary,
    )
    assert c1.flux_offset == pytest.approx(b.surface_factor * c0.flux_offset, rel=1e-13)
    assert c1.flux_slope == pytest.approx(b.surface_factor * c0.flux_slope, rel=1e-13)
    assert c1.surface_temperature_K == c0.surface_temperature_K

