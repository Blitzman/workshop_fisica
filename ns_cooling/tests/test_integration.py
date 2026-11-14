"""Tests of the integration area (decision D-028). Run them with `make test-integration` when this part of the code changes, or on request.

Code covered: cooling.py, thermal_solver.py, thermal_matrix.py, thermal_geometry.py, config.py, energy_balance.py.
Every test states a hypothesis, a reference and a criterion, in the style of the validation steps of the GitHub code.
"""

import numpy as np
import pytest
import thermal_matrix
import types
from config import CoolingConfig
from cooling import run_cooling
from dataclasses import replace
from energy_balance import energy_budget, energy_closure, thermal_energy_tables
from heat_capacity import heat_capacity
from scipy.linalg import expm, solve
from state import LocalState
from thermal_geometry import thermal_geometry
from thermal_matrix import assemble_matrix
from thermal_solver import solve_thermal_system


def test_default_scheme_closes_the_energy_balance(short_run, background):
    """Hypothesis: the default scheme conserves energy: dE/dt = -(L_nu + L_gamma).
    Reference: the thermal energy of the star, integrated independently of the solver
        (energy_budget).
    Criterion: largest closure error below 2.5 % in every window (about 2 % at the default step in
        the window with the arrival of the cooling wave, decision D-021).
    """
    # The largest error sits in the window with the arrival of the cooling wave at the surface (about
    # 400 yr), where the luminosity drops steeply; it halves with every halving of max_change.
    closure = energy_closure(energy_budget(short_run, background), start_years=1.0, tolerance=0.03)
    assert closure.ok, closure.format()
    assert closure.max_error < 0.025  # measured 0.9%; bounds a regression of the scheme


def test_thermal_energy_table_integrates_the_heat_capacity(background):
    """Hypothesis: the thermal-energy table satisfies dU/dT = C_v.
    Reference: centred differences of the table, against the heat capacity of the code.
    Criterion: relative 2e-3 in crust and core cells.
    """
    ln_t, energy = thermal_energy_tables(background)
    t = np.exp(ln_t)
    for i in (5, 80, 140, 160):  # crust and core cells
        dudt = np.gradient(energy[:, i], t)
        for j in (1200, 2000, 3000):
            cv = heat_capacity(LocalState(np.array([t[j]]), background.cells))["total"][i]
            assert dudt[j] == pytest.approx(cv, rel=2e-3), (i, j)


def test_an_inconsistent_solver_is_detected(background, monkeypatch):
    """Hypothesis: a solver that loses 10 % less neutrino energy than the diagnostic's L_nu breaks the
        energy check (negative control).
    Reference: the same run with a weakened source term.
    Criterion: the check fails, with an error above 5 %.
    """
    original = thermal_matrix.linearized_neutrino_source

    def weakened(state):
        q0, q1 = original(state)
        return 0.9 * q0, 0.9 * q1

    monkeypatch.setattr(thermal_matrix, "linearized_neutrino_source", weakened)
    run = run_cooling(background, replace(CoolingConfig(), final_years=300.0))
    closure = energy_closure(energy_budget(run, background), start_years=1.0, tolerance=0.02)
    assert not closure.ok and closure.max_error > 0.05


# ------------------------------------------------------------------ order of the theta-method
def linear_problem(background):
    """Constant-coefficient problem on the real geometry with an exact (matrix exponential) solution."""
    geo = thermal_geometry(background)
    n = geo.n_cells
    x = np.linspace(0, 1, n)
    cv, q0, q1 = 5.0 + 3 * np.sin(3 * x), -0.02 * (1 + x), 0.05 * (1 + 0.5 * np.cos(2 * x))
    k = 1e-2 * (1 + x[:-1])
    boundary = types.SimpleNamespace(flux_offset=-0.3, flux_slope=0.8)
    u0 = 1.0 + 0.5 * x
    probe = assemble_matrix(u0, 1.0, geo, cv, q0, q1, k, boundary, 1.0)
    c, volume = probe.face_conductance, geo.volumes
    M = np.zeros((n, n))
    b = -volume * q0
    for i in range(n):
        M[i, i] -= volume[i] * q1[i]
        if i < n - 1:
            M[i, i] -= c[i]
            M[i, i + 1] += c[i]
        if i > 0:
            M[i, i] -= c[i - 1]
            M[i, i - 1] += c[i - 1]
    M[-1, -1] -= geo.areas[-1] * boundary.flux_slope
    b[-1] -= geo.areas[-1] * boundary.flux_offset
    M /= (cv * volume)[:, None]
    b /= cv * volume
    total = 2.0 / np.sort(np.abs(np.linalg.eigvals(M)))[0]  # two slowest e-foldings
    exact = expm(M * total) @ u0 + solve(M, (expm(M * total) - np.eye(n)) @ b)
    return geo, (cv, q0, q1, k, boundary), u0, probe.dt_seconds, total, exact


@pytest.fixture(scope="module")
def long_run(background):
    return run_cooling(background, replace(CoolingConfig(), final_years=1e5))


def test_backward_euler_is_first_order_and_crank_nicolson_second_order(background):
    """Hypothesis: backward Euler is first order and Crank-Nicolson second order in time.
    Reference: the exact solution of the linearised problem.
    Criterion: error ratios of 2 (within 3 %) and 4 (within 5 %) when the step is halved, and a
        Crank-Nicolson error 1000 times smaller.
    """
    geo, coeffs, u0, seconds_per_year, total, exact = linear_problem(background)

    def error(theta, steps):
        u = u0.copy()
        dt_years = total / steps / seconds_per_year
        for _ in range(steps):
            u = solve_thermal_system(assemble_matrix(u, dt_years, geo, *coeffs, theta))
        return np.max(np.abs(u - exact) / np.abs(exact))

    # With enough steps to resolve the fastest modes the orders are exactly 1 and 2.
    be = [error(1.0, n) for n in (1024, 2048, 4096)]
    cn = [error(0.5, n) for n in (1024, 2048, 4096)]
    assert be[0] / be[1] == pytest.approx(2.0, rel=0.03) and be[1] / be[2] == pytest.approx(
        2.0, rel=0.03
    )
    assert cn[0] / cn[1] == pytest.approx(4.0, rel=0.05) and cn[1] / cn[2] == pytest.approx(
        4.0, rel=0.05
    )
    assert cn[2] < 1e-3 * be[2]


def test_cooling_is_monotonic_with_the_default_start(long_run):
    """Hypothesis: the backward-Euler start-up prevents the ringing of the outer crust.
    Reference: physical expectation: a cooling star never re-heats.
    Criterion: no profile increases in time and the luminosity decreases at every step.
    """
    u = np.array([long_run.profiles[n] for n in sorted(long_run.profiles)])
    assert not (np.diff(u, axis=0) > 0).any()
    assert np.all(np.diff(long_run.luminosity_erg_s) < 0)
