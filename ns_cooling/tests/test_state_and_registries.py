"""LocalState and the three registries (neutrino processes, conductivity channels, heat capacity)."""

import numpy as np
import pytest

from conductivity import CONDUCTIVITY_CHANNELS, thermal_transport
from emissivity import NEUTRINO_PROCESSES, neutrino_emissivity, process_factors
from heat_capacity import HEAT_CAPACITY_TERMS, heat_capacity
from registry import Registry
from state import LocalState


@pytest.fixture
def cell_state(background):
    return LocalState(np.full(background.n_cells, 3e8), background.cells)


# ------------------------------------------------------------------ LocalState
def test_state_broadcasts_a_temperature_grid_against_the_cells(background):
    state = LocalState(np.array([1e8, 2e8, 4e8])[:, None], background.cells)
    assert state.shape == (3, background.n_cells)
    assert state.rho.shape == state.shape and state.T.shape == state.shape


def test_state_rejects_invalid_temperatures(background):
    for bad in (0.0, -1.0, np.nan):
        with pytest.raises(ValueError):
            LocalState(np.full(background.n_cells, bad), background.cells)


def test_tau_follows_t_over_tc_and_is_ten_without_superfluidity(background):
    T = np.full(background.n_cells, 3e8)
    s = LocalState(T, background.cells)
    assert np.array_equal(s.tau_n, T / np.maximum(background.cells.tc_n, 0.1 * T))
    off = LocalState(T, background.cells, superfluid=False)
    assert np.all(off.tau_crust == 10.0) and np.all(off.tau_n == 10.0) and np.all(off.tau_p == 10.0)


def test_masks_split_crust_and_core(cell_state):
    assert np.array_equal(cell_state.crust, cell_state.xh > 0)
    assert np.array_equal(cell_state.core, ~cell_state.crust)
    assert cell_state.crust.any() and cell_state.core.any()


def test_with_temperature_keeps_tc_ratio_and_shares_fixed_data(cell_state):
    factors = process_factors(cell_state)
    hot = cell_state.with_temperature(1.01 * cell_state.T)
    assert np.allclose(hot.T, 1.01 * cell_state.T)
    assert np.array_equal(hot.tau_n, cell_state.tau_n)  # T/Tc is pinned
    assert process_factors(hot) is factors  # shared, not recomputed
    cell_state.memo("x", lambda: 1)
    assert "x" not in hot._memo  # temperature-dependent memo is not shared


# ------------------------------------------------------------------ Registry
def test_registry_order_duplicates_and_replacement():
    reg = Registry("thing")
    reg.register("a")(lambda: 1)
    reg.register("b")(lambda: 2)
    with pytest.raises(ValueError, match="already registered"):
        reg.register("a")(lambda: 3)
    reg.register("a", replace=True)(lambda: 4)
    assert reg.names() == ("a", "b") and reg.items()[0][1]() == 4  # replacement keeps the position
    reg.unregister("a")
    reg.register("a")(lambda: 5)
    assert reg.names() == ("b", "a")  # removing and registering again moves it to the end
    with pytest.raises(KeyError):
        reg.unregister("zzz")


def test_default_registry_contents_and_order():
    assert NEUTRINO_PROCESSES.names() == (
        "modified_urca", "brems_nn", "brems_np", "brems_pp", "brems_ep", "pbf_core",
        "direct_urca", "brems_e_ions", "plasma", "pbf_crust", "pair",
    )  # fmt: skip
    assert CONDUCTIVITY_CHANNELS.names() == (
        "core_electrons",
        "core_muons",
        "core_neutrons",
        "crust_electrons",
    )
    assert HEAT_CAPACITY_TERMS.names() == ("electrons", "muons", "neutrons", "protons", "ions")


# ------------------------------------------------------------------ new physics plugs in
def test_new_neutrino_process_adds_exactly_its_emissivity(cell_state):
    before = neutrino_emissivity(cell_state)["total"].copy()

    @NEUTRINO_PROCESSES.register("toy")
    def toy(state):
        return 1e20 * (state.T * 1e-9) ** 8 * (state.rho > 2.8e14)

    try:
        after = neutrino_emissivity(cell_state)
        assert np.allclose(after["total"] - before, toy(cell_state), rtol=1e-12)
        assert "toy" in after
    finally:
        NEUTRINO_PROCESSES.unregister("toy")
    assert np.array_equal(neutrino_emissivity(cell_state)["total"], before)


def test_new_conductivity_channel_and_heat_capacity_term_add_exactly(cell_state):
    k0 = thermal_transport(cell_state)["kappa_total"].copy()
    c0 = heat_capacity(cell_state)["total"].copy()

    @CONDUCTIVITY_CHANNELS.register("toy")
    def toy_k(state):
        return np.where(state.crust, 1e18, 0.0)

    @HEAT_CAPACITY_TERMS.register("toy")
    def toy_c(state):
        return 1e18 * state.T * (state.xh == 0)

    try:
        assert np.allclose(
            thermal_transport(cell_state)["kappa_total"] - k0, toy_k(cell_state), rtol=1e-9
        )
        assert np.allclose(heat_capacity(cell_state)["total"] - c0, toy_c(cell_state), rtol=1e-9)
    finally:
        CONDUCTIVITY_CHANNELS.unregister("toy")
        HEAT_CAPACITY_TERMS.unregister("toy")
    assert np.array_equal(thermal_transport(cell_state)["kappa_total"], k0)
    assert np.array_equal(heat_capacity(cell_state)["total"], c0)


# ------------------------------------------------------------------ microphysics sanity
def test_components_sum_to_total_and_are_physical(cell_state):
    for result in (neutrino_emissivity(cell_state), heat_capacity(cell_state)):
        parts = [v for k, v in result.items() if k != "total"]
        assert np.allclose(sum(parts), result["total"], rtol=1e-13)
        assert all(np.isfinite(v).all() and (v >= 0).all() for v in result.values())


def test_superfluid_heat_capacity_regimes(background):
    """Normal above Tc, strongly reduced deep inside the gap, and the (deliberate) jump just below Tc."""

    def ratio(T):
        T = np.full(background.n_cells, T)
        on = heat_capacity(LocalState(T, background.cells, superfluid=True))["total"]
        off = heat_capacity(LocalState(T, background.cells, superfluid=False))["total"]
        return on / off

    c = background.cells
    assert np.all(ratio(1e10) == 1.0)  # hotter than every Tc: identical to the normal state
    hot = ratio(2e8)
    deep = (c.xh == 0) & (c.tc_n > 3 * 2e8)
    assert deep.any() and np.all(hot[deep] < 0.5)
    near = (c.xh == 0) & (c.tc_n > 2e8) & (c.tc_n < 1.5 * 2e8)
    assert near.any() and np.all(hot[near] > 1.0)  # the heat-capacity jump at Tc is kept


def test_transport_is_evaluated_on_face_states(background):
    faces = background.faces.interior()
    k = thermal_transport(LocalState(np.full(len(faces), 3e8), faces))
    assert k["transport_internal"].shape == (background.n_faces - 2,)
    assert np.all(k["transport_internal"] > 0)
