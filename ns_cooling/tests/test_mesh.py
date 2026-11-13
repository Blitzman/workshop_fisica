"""Placement of the faces: exact interface and base, no slivers, caps, refinement."""

import numpy as np
import pytest

import mesh
import tov
from mesh import MeshSpec, describe_mesh, mesh_spec, place_faces


@pytest.fixture(scope="module")
def limits(eos_table, bsk22_star):
    pc = bsk22_star.pressure_central
    x_base = float(np.log(float(eos_table.pressure_at_rho_baryon(1e10)) / pc))
    x_interface = float(np.log(eos_table.pressure_core / pc))
    return x_base, x_interface


@pytest.fixture(scope="module")
def faces(bsk22_star, limits):
    return {n: place_faces(bsk22_star, limits[0], limits[1], mesh_spec(n)) for n in "ABCD"}


def test_presets_and_overrides():
    assert mesh_spec("A").dlnp_crust == 0.101 and mesh_spec("C").dr_core_max_m == 300.0
    assert (
        mesh_spec("B", dlnp_crust=0.04, blend=None).dlnp_crust == 0.04
        and mesh_spec("B", blend=None).blend == 0.5
    )
    with pytest.raises(ValueError, match="Unknown mesh"):
        mesh_spec("Z")
    with pytest.raises(ValueError, match="Unknown mesh parameters"):
        mesh_spec("A", spacing=1.0)


@pytest.mark.parametrize(
    "kwargs",
    [
        dict(dlnp_core=0.0),
        dict(dlnp_crust=-1.0),
        dict(refine=0.0),
        dict(blend=-0.1),
        dict(dr_core_max_m=0.0),
        dict(dlnp_core=np.nan),
    ],
)
def test_invalid_specs_are_refused(kwargs):
    with pytest.raises(ValueError):
        MeshSpec(**kwargs)


@pytest.mark.parametrize("name", list("ABCD"))
def test_faces_run_from_the_center_to_the_base(faces, limits, name):
    x = faces[name]
    assert x[0] == 0.0 and x[-1] == limits[0]
    assert np.all(np.diff(x) < 0) and np.all(np.isfinite(x))


@pytest.mark.parametrize("name", list("ABCD"))
def test_the_interface_is_an_exact_face(faces, limits, name):
    assert np.any(faces[name] == limits[1])


def test_no_sliver_cells_next_to_the_interface_or_the_base(faces):
    steps = -np.diff(faces["A"])
    assert steps.min() > 0.45 * 0.101 and steps.max() < 1.55 * 0.101  # the regular step is 0.101
    assert (-np.diff(faces["B"])).min() > 0.45 * 0.05


def test_face_counts(faces):
    assert [len(faces[n]) for n in "ABCD"] == pytest.approx([161, 274, 290, 581], abs=4)


def test_crust_is_twice_as_fine_in_B(bsk22_star, faces, limits):
    a, b = (describe_mesh(bsk22_star, faces[n], limits[1]) for n in "AB")
    assert b["n_crust_cells"] / a["n_crust_cells"] == pytest.approx(2.0, rel=0.1)
    assert b["crust_dr_m"][1] / a["crust_dr_m"][1] == pytest.approx(0.5, rel=0.1)
    assert b["n_core_cells"] == pytest.approx(
        a["n_core_cells"], abs=4
    )  # the core keeps its resolution


def test_radial_cap_in_the_core(bsk22_star, faces, limits):
    for name, cap in (("C", 300.0), ("D", 150.0)):
        dr = np.diff(bsk22_star.profile(faces[name])["r_cm"]) / 100.0
        core = 0.5 * (faces[name][1:] + faces[name][:-1]) > limits[1]
        core_dr = dr[core]
        assert np.sort(core_dr)[-2] <= cap * (
            1 + 1e-3
        )  # every core cell, but possibly the one at the interface
        assert core_dr.max() <= 1.5 * cap  # which may be moved by up to half a step
    assert (
        describe_mesh(bsk22_star, faces["A"], limits[1])["core_dr_m"][2] > 2000.0
    )  # without a cap: km-wide


def test_neighbouring_cells_vary_smoothly(bsk22_star, faces, limits):
    ratio = {
        n: describe_mesh(bsk22_star, faces[n], limits[1])["worst_neighbour_ratio"] for n in "ABCD"
    }
    assert ratio["A"] < 3.0 and ratio["C"] < 1.6 and ratio["D"] < 2.0


def test_refine_doubles_the_resolution(bsk22_star, limits):
    n1 = len(place_faces(bsk22_star, limits[0], limits[1], mesh_spec("A")))
    n2 = len(place_faces(bsk22_star, limits[0], limits[1], mesh_spec("A", refine=2.0)))
    assert 1.9 < n2 / n1 < 2.1


def test_describe_mesh_is_consistent(bsk22_star, faces, limits):
    d = describe_mesh(bsk22_star, faces["A"], limits[1])
    assert d["n_core_cells"] + d["n_crust_cells"] == d["n_faces"] - 1 == len(faces["A"]) - 1


def test_a_star_without_a_core_is_all_crust(eos_table):
    star = tov.integrate_star(eos_table, 1e31)  # below the pressure of the first core row
    x_base = float(np.log(float(eos_table.pressure_at_rho_baryon(1e10)) / 1e31))
    x = place_faces(star, x_base, np.inf, mesh_spec("B"))
    d = describe_mesh(star, x, np.inf)
    assert d["n_core_cells"] == 0 and d["n_crust_cells"] == len(x) - 1
    steps = -np.diff(x)
    assert np.allclose(steps[:-1], 0.05, rtol=1e-12)  # the crust step everywhere...
    assert 0.5 * 0.05 <= steps[-1] <= 1.5 * 0.05  # ...except the last one, which lands on the base


def test_the_base_must_lie_below_the_center(bsk22_star, limits):
    with pytest.raises(ValueError, match="below the center"):
        place_faces(bsk22_star, 0.5, limits[1], mesh_spec("A"))


def test_runaway_meshes_are_stopped(bsk22_star, limits, monkeypatch):
    monkeypatch.setattr(mesh, "_MAX_FACES", 50)
    with pytest.raises(RuntimeError, match="Too many faces"):
        place_faces(bsk22_star, limits[0], limits[1], mesh_spec("B"))
