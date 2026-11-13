"""Structure integration: analytic limit, the structure equations, matching and convergence."""

import numpy as np
import pytest

import tov
from constants import C, G, M_SUN, PRESSURE_TO_GEOM, RHO_TO_GEOM

from conftest import CENTRAL_PRESSURE_14


class Polytrope:
    """P = K rho^2 (n = 1) with e/c^2 = rho: the Newtonian limit has an analytic solution."""

    def __init__(self, K, rho_min, rho_max):
        self.K = K
        self.pressure_min, self.pressure_max = K * rho_min**2, K * rho_max**2

    def energy_density(self, pressure):
        return np.sqrt(np.asarray(pressure) / self.K)


R_ANALYTIC = 1.0e6  # cm; for n = 1 the radius is pi sqrt(K / (2 pi G)), independent of the mass
K_N1 = 2 * np.pi * G * (R_ANALYTIC / np.pi) ** 2


def polytrope_star(rho_c):
    eos = Polytrope(K_N1, 1e-6 * rho_c, rho_c * 1.0001)
    star = tov.integrate_star(eos, K_N1 * rho_c**2, pressure_outer=eos.pressure_min * 1.0001)
    mass_analytic = 4 / np.pi * rho_c * R_ANALYTIC**3  # Lane-Emden n = 1: M = 4 pi^2 rho_c alpha^3
    compactness = G * mass_analytic / (R_ANALYTIC * C**2)
    return star, mass_analytic, compactness


def test_newtonian_polytrope_has_the_analytic_radius_and_mass():
    deviations = []
    for rho_c in (1e11, 1e12, 1e13):
        star, mass, compactness = polytrope_star(rho_c)
        dr, dm = star.radius_cm / R_ANALYTIC - 1, star.mass_solar * M_SUN / mass - 1
        assert abs(dr) < 1.5 * compactness + 3e-5 and abs(dm) < 4 * compactness + 3e-5
        deviations.append(abs(dr))
    # the departure is the post-Newtonian correction: it grows in proportion to the compactness
    assert 7 < deviations[2] / deviations[1] < 13


def test_the_radius_of_the_polytrope_does_not_depend_on_the_mass():
    radii = [polytrope_star(rho_c)[0].radius_cm for rho_c in (3e10, 1e11)]
    assert radii[0] == pytest.approx(radii[1], rel=5e-5)


def clustered_profile(star, n=40001):
    x = (star.x_outer * np.linspace(0, 1, n) ** 4)[1:-1]
    return x, star.profile(x)


def test_mass_equation(eos_table, bsk22_star):
    """M(r) = integral of 4 pi r^2 e/c^2 dr, evaluated independently of the integration."""
    x, p = clustered_profile(bsk22_star)
    r = p["r_cm"]
    rho = eos_table.energy_density(bsk22_star.pressure_at_x(x))
    f = 4 * np.pi * r**2 * rho
    mass = np.concatenate([[0.0], np.cumsum(0.5 * (f[1:] + f[:-1]) * np.diff(r))])
    assert mass[-1] / p["mass_g"][-1] == pytest.approx(1.0, abs=1e-6)
    assert np.max(np.abs(mass[r > 1e5] / p["mass_g"][r > 1e5] - 1)) < 1e-6


def test_hydrostatic_equilibrium(eos_table, bsk22_star):
    """dP/dr from finite differences of the profile against the TOV expression.

    The agreement is checked in the bulk: isolated points next to a row of the tabulated
    EOS (a kink of the interpolated e(P)) are excluded from the strict bound, because the
    derivative of the dense-output polynomial is less accurate there than its values.
    """
    x, p = clustered_profile(bsk22_star)
    r, m = p["r_cm"], p["mass_g"] * G / C**2
    pressure = bsk22_star.pressure_at_x(x)
    pg, rg = pressure * PRESSURE_TO_GEOM, eos_table.energy_density(pressure) * RHO_TO_GEOM
    stride = 10
    i = np.arange(stride, len(r) - stride)
    dp = (pg[i + stride] - pg[i - stride]) / (r[i + stride] - r[i - stride])
    expected = (
        -(rg[i] + pg[i]) * (m[i] + 4 * np.pi * r[i] ** 3 * pg[i]) / (r[i] * (r[i] - 2 * m[i]))
    )
    residual = np.abs(dp / expected - 1)[r[i] > 1e5]
    assert np.median(residual) < 1e-4
    assert np.mean(residual > 1e-3) < 0.01


def test_metric_function_equation(bsk22_star):
    """d Phi / d r = (m + 4 pi r^3 P) / (r (r - 2 m)) in geometric units."""
    x, p = clustered_profile(bsk22_star)
    r, m = p["r_cm"], p["mass_g"] * G / C**2
    pg = bsk22_star.pressure_at_x(x) * PRESSURE_TO_GEOM
    sel = (r > 1e5) & (np.arange(len(r)) > 5) & (np.arange(len(r)) < len(r) - 5)
    numeric = np.gradient(np.log(p["ephi"]), r)
    expected = (m + 4 * np.pi * r**3 * pg) / (r * (r - 2 * m))
    assert np.max(np.abs(numeric[sel] / expected[sel] - 1)) < 1e-4


def test_center_values(bsk22_star):
    c = bsk22_star.profile(0.0)
    assert float(c["r_cm"][0]) == 0.0 and float(c["mass_g"][0]) == 0.0
    assert float(c["elam"][0]) == 1.0 and float(c["g_cm_s2"][0]) == 0.0
    assert 0.5 < float(c["ephi"][0]) < bsk22_star.surface_ephi


def test_surface_gravity_has_the_expected_size(bsk22_star):
    g_newton = G * bsk22_star.mass_solar * M_SUN / bsk22_star.radius_cm**2
    assert (
        1.0 < bsk22_star.surface_gravity_cm_s2 / g_newton < 1.5
    )  # relativistic enhancement e^Lambda (1 + ...)


def test_inverse_radius_mapping(bsk22_star):
    r = np.linspace(1e4, bsk22_star.radius_cm * 0.999, 9)
    back = bsk22_star.profile(bsk22_star.x_at_radius(r))["r_cm"]
    assert np.max(np.abs(back / r - 1)) < 1e-7


def test_profile_is_limited_to_the_integrated_range(bsk22_star):
    for bad in (bsk22_star.x_outer * 1.01, 0.1):
        with pytest.raises(ValueError, match="outside the integrated range"):
            bsk22_star.profile(bad)


def test_convergence_with_the_integration_tolerance(eos_table, bsk22_star):
    fine = tov.integrate_star(eos_table, CENTRAL_PRESSURE_14, rtol=1e-11, max_step=0.02)
    assert bsk22_star.mass_solar == pytest.approx(fine.mass_solar, rel=1e-6)
    assert bsk22_star.radius_cm == pytest.approx(fine.radius_cm, rel=1e-6)
    assert bsk22_star.surface_ephi == pytest.approx(fine.surface_ephi, rel=1e-6)
    assert bsk22_star.surface_gravity_cm_s2 == pytest.approx(fine.surface_gravity_cm_s2, rel=2e-6)


def test_impossible_target_masses_are_refused(eos_table):
    with pytest.raises(ValueError, match="exceeds the maximum mass"):
        tov.central_pressure_for_mass(eos_table, 2.5)
    with pytest.raises(ValueError, match="below the minimum mass 0.08"):
        tov.central_pressure_for_mass(eos_table, 0.05)
    for bad in (-1.0, 0.0, np.nan):
        with pytest.raises(ValueError):
            tov.central_pressure_for_mass(eos_table, bad)


def test_mass_radius_curve_matches_the_published_one(eos_table):
    """The curve of the authors (eos.mr) is the independent reference of the whole pipeline."""
    sweep = []
    for pc in np.geomspace(5e33, 1.1e36, 11):
        s = tov.integrate_star(eos_table, pc, rtol=1e-8, max_step=0.1)
        sweep.append((s.mass_solar, s.radius_cm / 1e5))
    sweep = np.array(sweep)
    stable = sweep[: np.argmax(sweep[:, 0]) + 1]
    mr = eos_table.mass_radius
    k = int(np.argmax(mr[:, 1]))
    assert stable[:, 0].max() == pytest.approx(mr[:, 1].max(), rel=2e-3)
    for mass in (0.8, 1.4, 2.0):
        ours = np.interp(mass, stable[:, 0], stable[:, 1])
        reference = np.interp(mass, mr[: k + 1, 1], mr[: k + 1, 0])
        assert ours == pytest.approx(
            reference, rel=3e-3
        )  # <= 0.13% observed; the outer layers are neglected


def test_invalid_settings_are_refused(eos_table):
    pc = CENTRAL_PRESSURE_14
    with pytest.raises(ValueError, match="Pmin <= Pouter"):
        tov.integrate_star(eos_table, pc, pressure_outer=pc * 2)
    with pytest.raises(ValueError, match="Pmin <= Pouter"):
        tov.integrate_star(eos_table, eos_table.pressure_max * 2)
    with pytest.raises(ValueError, match="Positive finite"):
        tov.integrate_star(eos_table, pc, rtol=0.0)
    with pytest.raises(ValueError, match="Center radius too large"):
        tov.integrate_star(eos_table, pc, center_radius_cm=1e5)


# ------------------------------------------------------------------ the stable branch
def test_the_stable_branch_of_bsk22(eos_table, branch):
    assert branch.mass_max == pytest.approx(
        eos_table.mass_radius[:, 1].max(), rel=1e-3
    )  # eos.mr: 2.265
    assert 0.08 < branch.mass_min < 0.09 and 9e32 < branch.pressure_min_mass < 1.0e33
    assert 1.0e36 < branch.pressure_max_mass < 1.03e36

    def mass(pc):
        return tov.integrate_star(eos_table, pc, rtol=1e-8, max_step=0.1).mass_solar

    for factor in (0.9, 1.1):  # both are extrema: the mass is larger/smaller on either side
        assert mass(branch.pressure_min_mass * factor) > branch.mass_min
        assert mass(branch.pressure_max_mass * factor) < branch.mass_max


def test_the_extended_stars_below_the_stable_branch(eos_table, branch):
    """BSk22 has a near-plateau of P between its last crust row and its first core row, so a star
    centred there has an enormous tenuous envelope: a valid but unstable solution (M falls with Pc).
    """
    assert eos_table.pressure_core < branch.pressure_min_mass
    star = tov.integrate_star(eos_table, eos_table.pressure_core)
    assert star.radius_cm / 1e5 > 1000.0 and star.mass_solar == pytest.approx(0.672, abs=0.005)


def test_light_targets_are_on_the_stable_branch(eos_table):
    for target, radius in ((0.2, 16.0), (0.5, 13.1)):
        pc = tov.central_pressure_for_mass(eos_table, target)
        star = tov.integrate_star(eos_table, pc)
        assert star.mass_solar == pytest.approx(
            target, abs=2e-6
        ) and star.radius_cm / 1e5 == pytest.approx(radius, abs=0.1)
        heavier = tov.integrate_star(eos_table, pc * 1.02, rtol=1e-8, max_step=0.1)
        assert heavier.mass_solar > star.mass_solar  # the mass grows with Pc: the stable side


