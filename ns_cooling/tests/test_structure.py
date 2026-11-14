"""Tests of the structure area (decision D-028). Run them with `make test-structure` when this part of the code changes, or on request.

Code covered: eos_table.py, tov.py, stellar_structure.py, mesh.py, background.py, input_model.py.
Every test states a hypothesis, a reference and a criterion, in the style of the validation steps of the GitHub code.
"""

import warnings

import numpy as np
import pytest
import tov
from background import composition_residuals
from config import CoolingConfig
from conftest import CENTRAL_PRESSURE_14
from constants import C, G, M_SUN, PRESSURE_TO_GEOM, RHO_TO_GEOM
from cooling import run_cooling
from dataclasses import replace
from eos_table import PHASE_CORE
from mesh import mesh_spec, place_faces
from stellar_structure import StructureOptions, build_stellar_model


def test_thermodynamic_identity_from_the_raw_columns(eos_dir, eos_table):
    """Hypothesis: the tabulated EOS obeys e + p = n_b mu_b at T = 0.
    Reference: the thermodynamic identity, evaluated from the raw columns independently of the
        reader.
    Criterion: relative deviation below 1e-6 in every row.
    """
    rows = np.array(
        [
            [float(v) for v in line.split()[:10]]
            for line in (eos_dir / "eos.thermo").read_text().splitlines()[1:]
        ]
    )
    m_n = float((eos_dir / "eos.thermo").read_text().split()[0])
    nb = eos_table.nb_fm3
    p, mu_b, e = rows[:, 3] * nb, rows[:, 5], nb * m_n * (1 + rows[:, 9])
    assert np.max(np.abs((e + p) / (nb * m_n * (1 + mu_b)) - 1)) < 1e-6


def test_rows_are_reproduced_exactly(eos_table):
    """Hypothesis: the reader returns the rows of the package unchanged when it is evaluated at their
        own pressures (it is also the smoke test that the package loads).
    Reference: the table itself.
    Criterion: densities to 1e-12, Ye to 1e-9, and the phase flag exactly.
    """
    state = eos_table.at_pressure(eos_table.pressure)
    np.testing.assert_allclose(state.rho_energy, eos_table.rho_energy, rtol=1e-12)
    np.testing.assert_allclose(state.nb_fm3, eos_table.nb_fm3, rtol=1e-12)
    np.testing.assert_allclose(state.ye, eos_table.ye, rtol=1e-9, atol=1e-15)
    np.testing.assert_array_equal(state.core, eos_table.phase == PHASE_CORE)


@pytest.fixture(scope="module")
def limits(eos_table, bsk22_star):
    pc = bsk22_star.pressure_central
    x_base = float(np.log(float(eos_table.pressure_at_rho_baryon(1e10)) / pc))
    x_interface = float(np.log(eos_table.pressure_core / pc))
    return x_base, x_interface


@pytest.fixture(scope="module")
def faces(bsk22_star, limits):
    return {n: place_faces(bsk22_star, limits[0], limits[1], mesh_spec(n)) for n in "ABCD"}


@pytest.mark.parametrize("name", ["A"])
def test_the_interface_is_an_exact_face(faces, limits, name):
    """Hypothesis: the crust-core interface is a face of the mesh, so that no cell mixes crust and
        core.
    Reference: the definition of the mesh.
    Criterion: one face equals the interface value exactly (mesh A).
    """
    assert np.any(faces[name] == limits[1])


PC_BSK26_14 = 1.2541308424984512e35  # dyn/cm^2: 1.4 Msun in BSk26 (input_bsk26)


@pytest.fixture(scope="module")
def bsk26_model(root, input_config):
    options = StructureOptions("BSk26", central_pressure=PC_BSK26_14)
    with warnings.catch_warnings():  # the table has a pressure reversal, reported by the reader
        warnings.simplefilter("ignore")
        return build_stellar_model(options, root / "eos", input_config)


# ------------------------------------------------------------------ the EOS stars
def test_core_charge_neutrality_and_baryon_number(bsk26_model):
    """Hypothesis: the star built from the EOS is charge neutral and conserves baryon number in the
        core and in the crust.
    Reference: the identities Yp = Ye + Ymu and Yn + Yp = 1 in the core, and Yn + Xh = 1 in the
        crust.
    Criterion: residuals below 1e-7 in the core and 1e-12 in the crust, on faces and cells (BSk26).
    """
    for points in (bsk26_model.background.faces, bsk26_model.background.cells):
        r = composition_residuals(points)
        core, crust = points.xh == 0, points.xh > 0
        assert np.abs(r["charge"][core]).max() < 1e-7  # yp = ye + ymu
        assert np.abs(r["baryon"][core]).max() < 1e-7  # yn + yp = 1
        assert np.abs(r["baryon"][crust]).max() < 1e-12  # yn + xh = 1 by construction


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


def clustered_profile(star, n=40001):
    x = (star.x_outer * np.linspace(0, 1, n) ** 4)[1:-1]
    return x, star.profile(x)


def test_newtonian_polytrope_has_the_analytic_radius_and_mass():
    """Hypothesis: for an n = 1 polytrope the integration reproduces the analytic radius and mass up to
        the post-Newtonian correction.
    Reference: R = pi sqrt(K / 2 pi G), independent of the mass, and the analytic mass.
    Criterion: deviations below 1.5 C + 3e-5 (radius) and 4 C + 3e-5 (mass), with C the compactness,
        and a radius deviation that grows in proportion to C.
    """
    deviations = []
    for rho_c in (1e11, 1e12, 1e13):
        star, mass, compactness = polytrope_star(rho_c)
        dr, dm = star.radius_cm / R_ANALYTIC - 1, star.mass_solar * M_SUN / mass - 1
        assert abs(dr) < 1.5 * compactness + 3e-5 and abs(dm) < 4 * compactness + 3e-5
        deviations.append(abs(dr))
    # the departure is the post-Newtonian correction: it grows in proportion to the compactness
    assert 7 < deviations[2] / deviations[1] < 13


def test_mass_equation(eos_table, bsk22_star):
    """Hypothesis: the mass of the profile is the integral of 4 pi r^2 e / c^2.
    Reference: an independent trapezoid integral of the EOS energy density over the profile.
    Criterion: agreement to 1e-6 at the surface and at every point beyond 1 km.
    """
    x, p = clustered_profile(bsk22_star)
    r = p["r_cm"]
    rho = eos_table.energy_density(bsk22_star.pressure_at_x(x))
    f = 4 * np.pi * r**2 * rho
    mass = np.concatenate([[0.0], np.cumsum(0.5 * (f[1:] + f[:-1]) * np.diff(r))])
    assert mass[-1] / p["mass_g"][-1] == pytest.approx(1.0, abs=1e-6)
    assert np.max(np.abs(mass[r > 1e5] / p["mass_g"][r > 1e5] - 1)) < 1e-6


def test_hydrostatic_equilibrium(eos_table, bsk22_star):
    """Hypothesis: the profile satisfies the Tolman-Oppenheimer-Volkoff equation.
    Reference: dP/dr from finite differences of the profile, against the TOV expression.
    Criterion: median deviation below 1e-4 and fewer than 1 % of the points above 1e-3.


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


def test_convergence_with_the_integration_tolerance(eos_table, bsk22_star):
    """Hypothesis: the star does not change when the integration is made tighter.
    Reference: a run with rtol = 1e-11 and a smaller step.
    Criterion: mass, radius and surface redshift within 1e-6, surface gravity within 2e-6.
    """
    fine = tov.integrate_star(eos_table, CENTRAL_PRESSURE_14, rtol=1e-11, max_step=0.02)
    assert bsk22_star.mass_solar == pytest.approx(fine.mass_solar, rel=1e-6)
    assert bsk22_star.radius_cm == pytest.approx(fine.radius_cm, rel=1e-6)
    assert bsk22_star.surface_ephi == pytest.approx(fine.surface_ephi, rel=1e-6)
    assert bsk22_star.surface_gravity_cm_s2 == pytest.approx(fine.surface_gravity_cm_s2, rel=2e-6)


def test_mass_radius_curve_matches_the_published_one(eos_table):
    """Hypothesis: the mass-radius curve of BSk22 reproduces the curve of the authors (the one
        structure test against a reference, decision D-022).
    Reference: eos.mr of the CompOSE package.
    Criterion: M_max within 0.2 % and the radii at 0.8, 1.4 and 2.0 Msun within 0.3 % (observed:
        0.03 % and 0.13 %).
    """
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


# ------------------------------------------------------------------ independence of the mesh
@pytest.fixture(scope="module")
def model_d(root, input_config):
    options = StructureOptions("BSk22", central_pressure=CENTRAL_PRESSURE_14, mesh=mesh_spec("D"))
    return build_stellar_model(options, root / "eos", input_config)


# ------------------------------------------------------------------ resolution of the Direct Urca edge
@pytest.fixture(scope="module")
def model_c(root, input_config):
    options = StructureOptions("BSk22", central_pressure=CENTRAL_PRESSURE_14, mesh=mesh_spec("C"))
    return build_stellar_model(options, root / "eos", input_config)


def test_effective_masses_come_from_the_eos_in_the_core(bsk22_model):
    """Hypothesis: in the core the effective masses are those of the package, with no fixed values
        (decisions D-006 and D-009).
    Reference: the columns of the package.
    Criterion: equality to 1e-12 on the core faces, and transport masses equal to them.
    """
    m = bsk22_model
    f = m.background.faces
    core = f.xh == 0
    i = int(np.argmax(~core)) - 1  # the interface face
    inside = np.arange(i)  # core faces away from the interface
    state = m.eos.at_pressure(m.star.pressure_central * np.exp(m.x_faces[inside]))
    np.testing.assert_allclose(f.mstar_n[inside], state.mstar_n, rtol=1e-12)
    np.testing.assert_allclose(f.mstar_p[inside], state.mstar_p, rtol=1e-12)
    n_crust = int((m.eos.phase != 0).sum())
    assert f.mstar_n[i] == pytest.approx(m.eos.mstar_n[n_crust], rel=1e-12)  # the first core row
    for name in ("transport_mstar_n", "mstar_n"):
        np.testing.assert_array_equal(getattr(f, name)[core], f.mstar_n[core])
    np.testing.assert_array_equal(f.transport_mstar_p[core], f.mstar_p[core])


def test_the_star_does_not_depend_on_the_mesh(bsk22_model, model_d):
    """Hypothesis: the surface quantities and the redshift profile do not depend on the mesh.
    Reference: the same star on mesh A and on the finest mesh D.
    Criterion: identical surface values; the redshift profile agrees to 1e-5.
    """
    a, d = bsk22_model.background, model_d.background
    assert d.n_faces > 3 * a.n_faces
    assert (a.g14, a.surface_radius_km, a.surface_ephi) == (
        d.g14,
        d.surface_radius_km,
        d.surface_ephi,
    )
    assert (
        a.radius_km[-1] == d.radius_km[-1]
        and a.faces.ephi[-1] == d.faces.ephi[-1]
        and a.surface_factor == d.surface_factor
    )
    interpolated = np.interp(a.radius_km, d.radius_km, d.faces.ephi)
    assert np.max(np.abs(interpolated / a.faces.ephi - 1)) < 1e-5


@pytest.mark.slow
def test_the_cooling_of_a_star_with_direct_urca_needs_the_edge_to_be_resolved(
    bsk22_model, model_c, model_d
):
    """Hypothesis: on the unresolved mesh A the late luminosity of a star with direct Urca is off by
        about 16 %, while mesh C (core cells below 300 m) agrees with the finest mesh D. Slow (17
        s): run it with `make test-slow`.
    Reference: the finest mesh D.
    Criterion: |L_C / L_D - 1| < 1 % and |L_A / L_D - 1| > 5 % at 100 and 300 yr.
    """
    config = replace(CoolingConfig(), final_years=300.0)
    runs = {
        n: run_cooling(m.background, config)
        for n, m in (("A", bsk22_model), ("C", model_c), ("D", model_d))
    }

    def luminosity(run, t):
        return float(
            np.exp(np.interp(np.log(t), np.log(run.times_years), np.log(run.luminosity_erg_s)))
        )

    for t in (100.0, 300.0):
        assert abs(luminosity(runs["C"], t) / luminosity(runs["D"], t) - 1) < 0.01
        assert abs(luminosity(runs["A"], t) / luminosity(runs["D"], t) - 1) > 0.05
