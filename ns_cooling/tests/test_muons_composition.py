"""Muons, step 1: composition, number densities, Fermi momenta and the consistency checks."""

import warnings

import numpy as np
import pytest

from background import (
    Points,
    check_composition,
    composition_residuals,
)
from state import LocalState
from synthetic import point_fields, points_from_table
from stellar_structure import StructureOptions, build_stellar_model

PC_BSK26_14 = 1.2541308424984512e35  # dyn/cm^2: 1.4 Msun in BSk26 (input_bsk26)
M_U = 1.66e-24


@pytest.fixture(scope="module")
def bsk26_model(root, input_config):
    options = StructureOptions("BSk26", central_pressure=PC_BSK26_14)
    with warnings.catch_warnings():  # the table has a pressure reversal, reported by the reader
        warnings.simplefilter("ignore")
        return build_stellar_model(options, root / "eos", input_config)


@pytest.fixture(params=["bsk22_model", "bsk26_model"])
def model(request):
    return request.getfixturevalue(request.param)


def kf_from_density(rho, fraction):
    return (3 * np.pi**2 * rho * fraction / M_U * 1e-39) ** (1 / 3)


# ------------------------------------------------------------------ the EOS stars
def test_core_charge_neutrality_and_baryon_number(model):
    for points in (model.background.faces, model.background.cells):
        r = composition_residuals(points)
        core, crust = points.xh == 0, points.xh > 0
        assert np.abs(r["charge"][core]).max() < 1e-7  # yp = ye + ymu
        assert np.abs(r["baryon"][core]).max() < 1e-7  # yn + yp = 1
        assert np.abs(r["baryon"][crust]).max() < 1e-12  # yn + xh = 1 by construction


def test_the_crust_charge_residual_by_region(model):
    """ye = xh zz/aa: exact in the outer crust, ~1e-4 in the inner crust (the interpolation of the
    product between table rows) and up to 8% (BSk22) or 15% (BSk26) in a boundary layer next to the
    interface, down to about half of its density (a property of the tables, where the nuclei give
    way to the core)."""
    for points in (model.background.faces, model.background.cells):
        r = composition_residuals(points)["charge_relative"]
        crust = np.flatnonzero(points.xh > 0)  # ordered from the interface outwards
        outer = crust[points.rho[crust] < 4e11]
        assert np.median(r[outer]) < 1e-7 and r[outer].max() < 1e-3
        bad = crust[r[crust] > 1e-3]
        assert bad.size <= 15 and np.all(points.rho[bad] > 0.4 * points.rho[crust[0]])
        assert r[crust].max() < 0.2


def test_charge_neutrality_in_momentum_space(model):
    for points in (model.background.faces, model.background.cells):
        core = points.xh == 0
        np.testing.assert_allclose(
            points.kf_p[core] ** 3, points.kf_e[core] ** 3 + points.kf_mu[core] ** 3, rtol=1e-6
        )
        s = LocalState(1e9, points)
        np.testing.assert_allclose(s.n_p[core], (s.n_e + s.n_mu)[core], rtol=1e-6)
        np.testing.assert_allclose((s.n_n + s.n_p)[core], s.nb[core], rtol=1e-6)  # baryon number


def test_every_fermi_momentum_comes_from_its_own_density(model):
    for points in (model.background.faces, model.background.cells):
        for species, kf in (("yn", "kf_n"), ("yp", "kf_p"), ("ye", "kf_e"), ("ymu", "kf_mu")):
            np.testing.assert_allclose(
                getattr(points, kf),
                kf_from_density(points.rho, getattr(points, species)),
                rtol=1e-6,
            )
        assert np.all(points.kf_p[points.xh > 0] == 0)  # no free protons in the crust


def test_muons_exist_only_in_the_dense_core(model):
    eos, f = model.eos, model.background.faces
    assert np.all(model.background.cells.ymu[model.background.cells.xh > 0] == 0)
    assert np.all(np.diff(f.ymu) <= 1e-12)  # the muon fraction falls outwards
    onset = (
        eos.nb_fm3[np.argmax(eos.ymu > 0)] * 1e39 * M_U
    )  # baryon mass density of the first row with muons
    c = model.background.cells
    assert np.all(c.ymu[c.rho < 0.97 * onset] == 0) and np.all(c.ymu[c.rho > 1.3 * onset] > 0)
    assert (
        f.ymu[0] > 0.02 and f.kf_mu[0] > 0.6
    )  # at the centre (BSk22: 0.062, 0.89 fm^-1; BSk26: 0.021, 0.69 fm^-1)


def test_the_check_reports_the_residuals(model):
    worst = check_composition(model.background.faces)
    assert (
        worst["core_charge"] < 1e-7
        and worst["baryon"] < 1e-7
        and worst["crust_charge_relative"] < 0.2
    )


# ------------------------------------------------------------------ synthetic points
def table(n=6, ymu=None, ye=0.06):
    z = np.zeros(n)
    ymu_values = np.zeros(n) if ymu is None else np.full(n, ymu)
    t = dict(
        radius_km=np.arange(n, dtype=float),
        rho=np.linspace(9e14, 5e14, n),
        ye=np.full(n, ye),
        yn=np.full(n, 1 - ye - ymu_values[0]),
        yp=np.full(n, ye) + ymu_values,
        xh=z,
        aa=np.ones(n),
        zz=np.full(n, ye),
        ephi=np.full(n, 0.8),
        elam=np.full(n, 1.2),
        tc_crust=z,
        tc_n=z,
        tc_p=z,
    )
    return t, ymu_values


def test_invalid_muon_fractions_are_refused():
    fields = point_fields(table()[0])
    for bad, message in ((-0.01, "ymu must lie"), (1.5, "ymu must lie")):
        with pytest.raises(ValueError, match=message):
            Points(**dict(fields, ymu=np.full(6, bad)))
    # more muons than protons: no charge to neutralize
    with pytest.raises(ValueError, match="cannot exceed yp"):
        Points(**dict(fields, ymu=np.full(6, 0.2)))


# ------------------------------------------------------------------ the checks themselves
def test_check_composition_refuses_broken_cores():
    t, mu = table(ymu=0.03)
    t["ymu"] = mu
    good = points_from_table(t)
    assert check_composition(good)["core_charge"] < 1e-12
    broken = dict(t, yp=t["yp"] + 1e-3)  # a proton excess no electron or muon neutralizes
    with pytest.raises(ValueError, match="charge neutrality in the core: residual 1.000e-03"):
        check_composition(points_from_table(broken))
    lost = dict(t, yn=t["yn"] - 1e-3)  # baryons that disappeared
    with pytest.raises(ValueError, match="baryon number"):
        check_composition(points_from_table(lost))


def test_a_grossly_inconsistent_crust_is_reported():
    t, _ = table()
    crust = dict(
        t,
        xh=np.ones(6),
        aa=np.full(6, 56.0),
        zz=np.full(6, 26.0),
        yn=np.zeros(6),
        yp=np.zeros(6),
        ye=np.full(6, 0.2),
    )  # a crust with ye = 0.2 where 26/56 = 0.464 is expected
    points = points_from_table(crust)
    with pytest.warns(
        UserWarning, match="Crust charge neutrality: ye differs from xh\\*zz/aa by 132%"
    ):
        worst = check_composition(points)
    assert worst["crust_charge_relative"] == pytest.approx(
        abs(0.2 - 26 / 56) / 0.2, rel=1e-12
    )  # relative to ye
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        check_composition(points, crust_charge_warning=2.0)  # a looser limit accepts the 132%
