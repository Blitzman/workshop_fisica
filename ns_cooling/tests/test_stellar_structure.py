"""Stellar model from an EOS: options, input file, conventions of the background, surface values."""

import shutil
from dataclasses import replace

import numpy as np
import pytest

from background import Background, derived_fermi_quantities
from constants import C, G, M_SUN
from config import CoolingConfig
from cooling import run_cooling
from energy_balance import energy_budget, energy_closure
from input_model import load_background, load_backgrounds, read_parameters
import tov
from mesh import mesh_spec
from state import LocalState
from stellar_structure import StructureOptions, build_stellar_model, direct_urca_boundary_widths

from conftest import CENTRAL_PRESSURE_14

PC_LINE = f"central_pressure = {CENTRAL_PRESSURE_14}"


INPUT_HEADER = """Structure file (not used)
none
Initial temperature T_8 [10^8 K]
1.d2
Unused line
800
Superfluidity
a) crust neutrons
1 7
b) core neutrons
1 12
c) core protons
1 8
envelope_model = accreted"""


def write_input(tmp_path, root, extra, copy_eos=False):
    path = tmp_path / "input_test"
    path.write_text(INPUT_HEADER + "\n" + extra + "\n")
    if copy_eos:
        shutil.copytree(root / "eos" / "BSk22", tmp_path / "eos" / "BSk22")
    return path


# ------------------------------------------------------------------ options and input file
def test_options_need_exactly_one_of_mass_and_central_pressure():
    with pytest.raises(ValueError, match="exactly one"):
        StructureOptions("BSk22")
    with pytest.raises(ValueError, match="exactly one"):
        StructureOptions("BSk22", mass_solar=1.4, central_pressure=1e35)
    with pytest.raises(ValueError, match="positive"):
        StructureOptions("BSk22", mass_solar=-1.0)
    with pytest.raises(ValueError, match="base_density"):
        StructureOptions("BSk22", mass_solar=1.4, base_density=0.0)


def test_options_from_a_configuration():
    o = StructureOptions.from_config(
        dict(eos="BSk22", mass=1.4, mesh="C", refine=2.0, dlnp_crust=0.04)
    )
    assert (
        o.mesh.name == "C"
        and o.mesh.refine == 2.0
        and o.mesh.dlnp_crust == 0.04
        and o.mesh.dr_core_max_m == 300.0
    )
    assert StructureOptions.from_config(dict(eos="BSk22", mass=1.4)).mesh == mesh_spec("A")
    with pytest.raises(ValueError, match="no 'eos'"):
        StructureOptions.from_config(dict(mass=1.4))


def test_input_file_options_are_parsed(root, tmp_path):
    cfg = read_parameters(
        write_input(
            tmp_path, root, "eos = BSk22\nmass = 1.4\nmesh = b\nrefine = 2\nbase_density = 2e10"
        )
    )
    assert (
        cfg["eos"] == "BSk22" and cfg["mass"] == (1.4,) and cfg["mesh"] == "B"
    )  # the folder keeps its case
    assert (
        cfg["refine"] == 2.0 and cfg["base_density"] == 2e10 and cfg["envelope_model"] == "accreted"
    )


@pytest.mark.parametrize(
    "extra, message",
    [
        ("foo = 1", "Unknown input option"),
        ("eos = BSk22\nmass = 1.4\nmass = 1.5", "Duplicate mass"),
        ("eos = BSk22\nmass = abc", "needs a number"),
        ("eos = BSk22", "exactly one"),
        ("eos = BSk22\nmass = 1.4\ncentral_pressure = 1e35", "exactly one"),
        ("mass = 1.4", "no 'eos'"),
        ("eos = BSk22\nmass = 1.4\nmesh = Q", "Unknown mesh"),
        ("eos =\nmass = 1.4", "needs a value"),
    ],
)
def test_input_file_errors(root, tmp_path, extra, message):
    with pytest.raises(ValueError, match=message):
        read_parameters(write_input(tmp_path, root, extra))


@pytest.mark.filterwarnings("ignore:The Direct Urca")
def test_load_background_builds_the_star_from_the_input_file(root, tmp_path, bsk22_model):
    path = write_input(tmp_path, root, f"eos = BSk22\n{PC_LINE}\nmesh = A", copy_eos=True)
    background = load_background(path)
    assert background.n_faces == bsk22_model.background.n_faces
    np.testing.assert_array_equal(background.faces.rho, bsk22_model.background.faces.rho)
    np.testing.assert_array_equal(background.radius_km, bsk22_model.background.radius_km)
    assert "BSk22" in background.description and background.envelope_model == "accreted"


# ------------------------------------------------------------------ the background
def test_geometry_and_regions(bsk22_model):
    b = bsk22_model.background
    f, c = b.faces, b.cells
    assert (
        b.radius_km[0] == 0 and b.n_faces == len(bsk22_model.x_faces) and b.n_cells == b.n_faces - 1
    )
    crust = f.xh > 0
    i = int(np.argmax(crust)) - 1  # the last core face is the interface face
    assert f.xh[i] == 0 and np.all(f.xh[: i + 1] == 0) and np.all(crust[i + 1 :])
    assert np.all(c.xh[:i] == 0) and np.all(c.xh[i:] > 0)  # no cell mixes crust and core values
    assert f.rho[-1] == pytest.approx(1e10, rel=1e-9)


def test_composition_conventions(bsk22_model):
    f = bsk22_model.background.faces
    core = f.xh == 0
    assert np.all(f.aa[core] == 1) and np.allclose(f.zz[core], f.ye[core])
    assert np.allclose(f.xh[~core], 1 - f.yn[~core], rtol=1e-12) and np.all(f.yp[~core] == 0)
    assert np.allclose(f.yn[core] + f.yp[core], 1, atol=1e-6)
    assert (
        np.all(f.yp[core] >= f.ye[core] - 1e-12) and f.yp[0] > 1.3 * f.ye[0]
    )  # muons at the centre
    assert np.all(f.aa[~core] > 0) and np.all(f.zz[~core] > 0)


def test_rho_is_the_baryon_density_of_the_microphysics(bsk22_model):
    b = bsk22_model.background
    for points in (b.faces, b.cells):
        nn = (
            points.rho * points.yn / 1.66e-24
        )  # n = rho / 1.66e-24 is how the microphysics converts
        kf_n = (3 * np.pi**2 * nn * 1e-39) ** (1 / 3)
        np.testing.assert_allclose(points.kf_n, kf_n, rtol=1e-6)
        ne = points.rho * points.ye / 1.66e-24
        np.testing.assert_allclose(points.kf_e, (3 * np.pi**2 * ne * 1e-39) ** (1 / 3), rtol=1e-6)
    eos = bsk22_model.eos
    centre = eos.at_pressure(bsk22_model.star.pressure_central)
    assert b.faces.rho[0] == pytest.approx(float(centre.nb_fm3) * 1e39 * 1.66e-24, rel=1e-12)


def test_cell_values_lie_between_the_neighbouring_faces(bsk22_model):
    b = bsk22_model.background
    f, c = b.faces, b.cells
    assert np.all(f.rho[:-1] > c.rho) and np.all(c.rho > f.rho[1:])  # the density falls outwards
    assert np.all(f.ephi[:-1] < c.ephi) and np.all(c.ephi < f.ephi[1:])  # exp(Phi) grows outwards
    # exp(Lambda) is not monotonic (it peaks inside the star); the mass function m = r (1 - exp(-2 Lambda)) / 2 is
    assert b.elam[0] == 1.0 and np.all(b.elam >= 1.0)
    mass = 0.5 * (b.radius_km * 1e5) * (1 - 1 / b.elam**2) * C**2 / G / M_SUN
    assert np.all(np.diff(mass) >= 0) and mass[-1] < bsk22_model.star.mass_solar
    assert (
        1 - mass[-1] / bsk22_model.star.mass_solar < 1e-5
    )  # the layers above the base carry almost no mass


def test_effective_masses_come_from_the_eos_in_the_core(bsk22_model):
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


def test_the_crust_keeps_the_tabulated_effective_masses(bsk22_model):
    f = bsk22_model.background.faces
    crust = f.xh > 0
    derived = derived_fermi_quantities(f.rho, f.ye, f.yn, f.yp)
    for name in ("mstar_n", "mstar_p", "transport_mstar_n", "transport_mstar_p"):
        np.testing.assert_array_equal(getattr(f, name)[crust], derived[name][crust])


def test_critical_temperatures_follow_the_regions(bsk22_model, input_config):
    f = bsk22_model.background.faces
    core, crust = f.xh == 0, f.xh > 0
    assert (
        np.all(f.tc_n[crust] == 0) and np.all(f.tc_p[crust] == 0) and np.all(f.tc_crust[core] == 0)
    )
    assert f.tc_n[core].max() > 1e8 and f.tc_p[core].max() > 1e8 and f.tc_crust[crust].max() > 1e8
    tc = pairing_tc(f, input_config)
    for name in ("tc_crust", "tc_n", "tc_p"):
        np.testing.assert_array_equal(getattr(f, name), tc[name])


def pairing_tc(points, cfg):
    from superfluid import pairing_critical_temperatures

    return pairing_critical_temperatures(
        cfg, dict(rho=points.rho, yn=points.yn, yp=points.yp, xh=points.xh)
    )


# ------------------------------------------------------------------ surface values
def test_surface_values_come_from_the_outer_edge(bsk22_model):
    b, s = bsk22_model.background, bsk22_model.star
    assert b.g14 == pytest.approx(s.surface_gravity_cm_s2 / 1e14, rel=1e-15)
    assert (
        b.surface_radius_km == pytest.approx(s.radius_cm / 1e5, rel=1e-15)
        and b.surface_ephi == s.surface_ephi
    )
    assert b.surface_radius_km > b.radius_km[-1] and b.surface_ephi > b.faces.ephi[-1]
    expected = ((b.surface_radius_km * b.surface_ephi) / (b.radius_km[-1] * b.faces.ephi[-1])) ** 2
    assert b.surface_factor == pytest.approx(expected, rel=1e-14) and 1.0 < b.surface_factor < 1.06


# ------------------------------------------------------------------ independence of the mesh
@pytest.fixture(scope="module")
def model_d(root, input_config):
    options = StructureOptions("BSk22", central_pressure=CENTRAL_PRESSURE_14, mesh=mesh_spec("D"))
    return build_stellar_model(options, root / "eos", input_config)


def test_the_star_does_not_depend_on_the_mesh(bsk22_model, model_d):
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


# ------------------------------------------------------------------ errors
def test_the_base_density_must_lie_inside_the_star(root, input_config):
    for rho_base in (1.5e15, 1e5):  # above the central density / below the table
        options = StructureOptions(
            "BSk22", central_pressure=CENTRAL_PRESSURE_14, base_density=rho_base
        )
        with pytest.raises(ValueError, match="density"):
            build_stellar_model(options, root / "eos", input_config)


def test_an_unknown_eos_folder_is_reported(root, input_config):
    with pytest.raises(FileNotFoundError):
        build_stellar_model(
            StructureOptions("Nope", central_pressure=1e35), root / "eos", input_config
        )


def test_the_csv_of_an_eos_background_is_for_inspection_only(bsk22_model, tmp_path):
    path = tmp_path / "background.csv"
    bsk22_model.background.to_csv(path)
    header = "\n".join(path.read_text().splitlines()[:8])
    assert (
        "exact = False" in header
        and "description = EOS BSk22" in header
        and "surface_radius_km" in header
    )
    with pytest.raises(ValueError, match="inspection only"):
        Background.from_csv(path)


# ------------------------------------------------------------------ the whole pipeline
def test_cooling_runs_on_the_eos_star_and_conserves_energy(bsk22_model):
    background = bsk22_model.background
    run = run_cooling(background, replace(CoolingConfig(), final_years=300.0))
    assert (
        run.stop_reason == "final_time"
        and np.all(np.isfinite(run.luminosity_erg_s))
        and np.all(run.luminosity_erg_s > 0)
    )
    assert run.luminosity_erg_s[0] > run.luminosity_erg_s[-1]
    closure = energy_closure(energy_budget(run, background), start_years=1.0, tolerance=0.03)
    assert closure.ok, closure.format()


# ------------------------------------------------------------------ resolution of the Direct Urca edge
@pytest.fixture(scope="module")
def model_c(root, input_config):
    options = StructureOptions("BSk22", central_pressure=CENTRAL_PRESSURE_14, mesh=mesh_spec("C"))
    return build_stellar_model(options, root / "eos", input_config)


def test_an_unresolved_direct_urca_edge_is_reported(root, input_config):
    options = StructureOptions("BSk22", central_pressure=CENTRAL_PRESSURE_14, mesh=mesh_spec("A"))
    with pytest.warns(UserWarning, match="Direct Urca threshold lies in cells up to 6[0-9][0-9] m"):
        build_stellar_model(options, root / "eos", input_config)


def test_a_resolved_direct_urca_edge_is_silent(root, input_config, bsk22_model, model_c, model_d):
    assert direct_urca_boundary_widths(bsk22_model.background) == pytest.approx(
        [0.66], abs=0.03
    )  # mesh A
    for model, width in ((model_c, 0.30), (model_d, 0.15)):
        assert direct_urca_boundary_widths(model.background) == pytest.approx([width], abs=0.02)
    import warnings

    with warnings.catch_warnings():
        warnings.simplefilter("error")
        build_stellar_model(model_c.options, root / "eos", input_config)


def test_a_star_without_direct_urca_has_no_edge_to_resolve(root, input_config):
    import warnings

    options = StructureOptions(
        "BSk22", central_pressure=2e34, mesh=mesh_spec("A")
    )  # below the threshold central pressure
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        model = build_stellar_model(options, root / "eos", input_config)
    assert direct_urca_boundary_widths(model.background).size == 0


def test_boundary_widths_of_a_synthetic_region():
    from test_urca_masses_surface import YE_SCAN, npe_background

    widths = direct_urca_boundary_widths(
        npe_background(YE_SCAN)
    )  # cells are 1 km wide; Direct Urca on for the last four faces
    np.testing.assert_allclose(widths, [1.0])


def test_the_cooling_of_a_star_with_direct_urca_needs_the_edge_to_be_resolved(
    bsk22_model, model_c, model_d
):
    """Why the guard exists: on the unresolved mesh A the late luminosity is off by ~16%;
    mesh C (core cells <= 300 m) agrees with the finest mesh D to better than 1%."""
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


# ------------------------------------------------------------------ several stars in one input file
def test_lists_of_stars_are_parsed(root, tmp_path):
    cfg = read_parameters(
        write_input(tmp_path, root, "eos = BSk22\nmass = 1.0, 1.4 1.8")
    )  # commas or blanks
    assert cfg["mass"] == (1.0, 1.4, 1.8)
    cfg = read_parameters(write_input(tmp_path, root, "eos = BSk22\ncentral_pressure = 4e34, 8e34"))
    assert cfg["central_pressure"] == (4e34, 8e34)
    options = StructureOptions.list_from_config(cfg)
    assert [o.central_pressure for o in options] == [4e34, 8e34] and all(
        o.mass_solar is None for o in options
    )
    assert len({o.mesh for o in options}) == 1  # the same mesh for all the stars
    assert [
        o.mass_solar for o in StructureOptions.list_from_config(dict(eos="BSk22", mass=(1.0, 1.8)))
    ] == [1.0, 1.8]


def test_several_stars_need_list_from_config(root, tmp_path):
    cfg = read_parameters(write_input(tmp_path, root, "eos = BSk22\nmass = 1.0, 1.4"))
    with pytest.raises(ValueError, match="defines 2 stars"):
        StructureOptions.from_config(cfg)
    with pytest.raises(ValueError, match="exactly one"):
        StructureOptions.list_from_config(dict(eos="BSk22", mass=(1.0,), central_pressure=(1e35,)))
    with pytest.raises(ValueError, match="needs a number"):
        read_parameters(write_input(tmp_path, root, "eos = BSk22\nmass = 1.0, x"))


def test_load_background_refuses_several_stars_before_building_them(root):
    with pytest.raises(ValueError, match="defines 3 stars; use load_backgrounds"):
        load_background(root / "input_bsk22")  # raises at once: nothing is integrated


@pytest.mark.filterwarnings("ignore:The Direct Urca")
def test_the_three_default_stars_and_their_direct_urca(root):
    stars = load_backgrounds(root / "input_bsk22")
    assert [label for label, _ in stars] == ["1.00 Msun", "1.40 Msun", "1.80 Msun"]
    backgrounds = [b for _, b in stars]
    allowed = [LocalState(1e9, b.cells).direct_urca_allowed for b in backgrounds]
    assert [bool(a.any()) for a in allowed] == [
        False,
        True,
        True,
    ]  # off below ~1.15 Msun, on above it
    radius = [b.radius_km[int(a.sum())] if a.any() else 0.0 for b, a in zip(backgrounds, allowed)]
    assert radius[1] == pytest.approx(5.3, abs=0.4) and radius[2] == pytest.approx(
        8.1, abs=0.5
    )  # km (exactly 5.5 and 8.2)
    assert backgrounds[0].g14 < backgrounds[1].g14 < backgrounds[2].g14
    assert all(b.surface_factor > 1.0 for b in backgrounds)
    assert [b.n_faces for b in backgrounds] == pytest.approx([154, 161, 168], abs=3)
