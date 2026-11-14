"""Energy balance dE/dt = -(L_nu + L_gamma): the diagnostic and what it detects."""

from dataclasses import replace

import numpy as np
import pytest

import thermal_matrix
from config import CoolingConfig
from cooling import run_cooling
from emissivity import NEUTRINO_PROCESSES
from energy_balance import energy_budget, energy_closure, thermal_energy_tables
from heat_capacity import heat_capacity
from state import LocalState


def test_default_scheme_closes_the_energy_balance(short_run, background):
    # The largest error sits in the window with the arrival of the cooling wave at the surface (about
    # 400 yr), where the luminosity drops steeply; it halves with every halving of max_change.
    closure = energy_closure(energy_budget(short_run, background), start_years=1.0, tolerance=0.03)
    assert closure.ok, closure.format()
    assert closure.max_error < 0.025  # measured 0.9%; bounds a regression of the scheme


def test_plain_backward_euler_is_flagged(background):
    run = run_cooling(
        background,
        replace(CoolingConfig(), final_years=1000.0, theta=1.0, predictor_corrector=False),
    )
    closure = energy_closure(energy_budget(run, background), start_years=1.0, tolerance=0.02)
    assert not closure.ok  # first-order lag of about 4-10% per decade
    assert all(w.closure_error > 0 for w in closure.windows)  # the star loses less than L says


def test_budget_is_consistent(short_run, background):
    b = energy_budget(short_run, background)
    assert len(b.times_years) == short_run.completed_steps + 1
    assert b.residual[0] == 0.0 and b.neutrino_energy[0] == 0.0
    assert np.all(np.diff(b.thermal_energy) < 0)  # the star only cools
    assert np.all(b.neutrino_luminosity > 0) and np.all(b.photon_luminosity > 0)
    # photon luminosity agrees with the luminosity curve recorded by the run itself
    late = b.times_years > 0.01
    assert np.array_equal(b.photon_luminosity[late], short_run.luminosity_erg_s)
    names, columns = b.as_columns()
    assert columns.shape == (len(b.times_years), len(names))


def test_window_edges_are_accepted_steps_and_cover_decades(short_run, background):
    c = energy_closure(energy_budget(short_run, background), start_years=1.0)
    times = set(short_run.step_times_years)
    assert all(w.start_years in times and w.end_years in times for w in c.windows)
    assert [round(np.log10(w.start_years)) for w in c.windows] == [0, 1, 2]


def test_thermal_energy_table_integrates_the_heat_capacity(background):
    """U(T) must satisfy dU/dT = C_v: compare a centered difference of the table with C_v."""
    ln_t, energy = thermal_energy_tables(background)
    t = np.exp(ln_t)
    for i in (5, 80, 140, 160):  # crust and core cells
        dudt = np.gradient(energy[:, i], t)
        for j in (1200, 2000, 3000):
            cv = heat_capacity(LocalState(np.array([t[j]]), background.cells))["total"][i]
            assert dudt[j] == pytest.approx(cv, rel=2e-3), (i, j)


def test_budget_requires_every_profile(background):
    run = run_cooling(background, replace(CoolingConfig(), final_years=10.0), snapshots={0})
    with pytest.raises(ValueError, match="every accepted step"):
        energy_budget(run, background)


# ------------------------------------------------------------------ the check has teeth, and new physics passes
def test_new_physics_passes_the_check_automatically(background):
    """A new neutrino process is part of both the solver and the diagnostic: the balance still closes."""

    @NEUTRINO_PROCESSES.register("toy")
    def toy(state):
        return 2e21 * (state.T * 1e-9) ** 8 * (state.xh == 0)  # about twice modified Urca

    try:
        run = run_cooling(background, replace(CoolingConfig(), final_years=300.0))
        closure = energy_closure(energy_budget(run, background), start_years=1.0, tolerance=0.02)
    finally:
        NEUTRINO_PROCESSES.unregister("toy")
    assert closure.ok, closure.format()


def test_an_inconsistent_solver_is_detected(background, monkeypatch):
    """If the solver loses 10% less neutrino energy than the diagnostic's L_nu, the check must fail."""
    original = thermal_matrix.linearized_neutrino_source

    def weakened(state):
        q0, q1 = original(state)
        return 0.9 * q0, 0.9 * q1

    monkeypatch.setattr(thermal_matrix, "linearized_neutrino_source", weakened)
    run = run_cooling(background, replace(CoolingConfig(), final_years=300.0))
    closure = energy_closure(energy_budget(run, background), start_years=1.0, tolerance=0.02)
    assert not closure.ok and closure.max_error > 0.05


@pytest.mark.parametrize("value", [0.0, -1.0, float("nan"), float("inf"), -float("inf")])
@pytest.mark.parametrize("option", ["start_years", "tolerance"])
def test_invalid_energy_options_fail_before_accessing_budget(option, value):
    with pytest.raises(ValueError, match="positive and finite"):
        energy_closure(None, **{option: value})


def test_closure_avoids_cancellation_from_initial_radiation():
    from energy_balance import EnergyBudget
    from thermal_matrix import SECONDS_PER_YEAR

    # Each late window loses exactly one energy unit per year. A huge
    # unresolved initial ledger must not erase that well-resolved balance.
    times = np.array([0.0, 1.0, 10.0, 100.0])
    budget = EnergyBudget(
        times,
        200.0 - times,
        np.ones(4) / SECONDS_PER_YEAR,
        np.zeros(4),
        np.array([0.0, 1e30, 1e30, 1e30]),
        np.zeros(4),
    )
    closure = energy_closure(budget)
    assert len(closure.windows) == 2
    assert closure.max_error < 1e-14


def test_cli_rejects_invalid_energy_options_before_loading_model(root, tmp_path):
    import subprocess
    import sys

    for option in ("--energy-start-years", "--energy-tolerance"):
        result = subprocess.run(
            [
                sys.executable,
                str(root / "cooling_driver.py"),
                "--check-energy",
                "--input",
                str(tmp_path / "missing_input"),
                option,
                "0",
                "--no-animation",
            ],
            capture_output=True,
            text=True,
        )
        assert result.returncode == 2
        assert "positive and finite" in result.stderr
