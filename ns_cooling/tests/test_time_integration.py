"""Time integration: order of the theta-method, default scheme, regression and stability."""

import types
from dataclasses import replace

import numpy as np
import pytest
from scipy.linalg import expm, solve

from config import CoolingConfig
from cooling import run_cooling
from thermal_geometry import thermal_geometry
from thermal_matrix import assemble_matrix
from thermal_solver import solve_thermal_system


# ------------------------------------------------------------------ configuration
def test_default_scheme_is_crank_nicolson_predictor_corrector_with_backward_euler_start():
    cfg = CoolingConfig()
    assert cfg.theta == 0.5 and cfg.predictor_corrector and cfg.fixed_steps_backward_euler
    assert cfg.n_fixed_steps == 20 and cfg.dt_initial_years == 0.01


@pytest.mark.parametrize(
    "bad",
    [dict(theta=0.4), dict(theta=1.1), dict(initial_temperature_8=0.0), dict(final_years=0.001),
     dict(minimum_temperature_K=-1.0), dict(n_fixed_steps=-1)],
)  # fmt: skip
def test_invalid_configurations_are_rejected(bad):
    with pytest.raises(ValueError):
        CoolingConfig(**bad)


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


def test_backward_euler_is_first_order_and_crank_nicolson_second_order(background):
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


def test_theta_must_lie_between_crank_nicolson_and_backward_euler(background):
    geo, coeffs, u0, *_ = linear_problem(background)
    with pytest.raises(ValueError):
        assemble_matrix(u0, 1.0, geo, *coeffs, 0.3)


# ------------------------------------------------------------------ regression and stability
def test_run_bookkeeping(short_run):
    assert short_run.stop_reason == "final_time"
    t = short_run.step_times_years
    assert (
        t[0] == 0.0 and t[-2] < 1000.0 <= t[-1]
    )  # the first step that reaches the final time ends the run
    assert len(short_run.step_times_years) == short_run.completed_steps + 1
    assert np.all(np.diff(short_run.step_times_years) > 0)
    assert set(short_run.profiles) == set(range(short_run.completed_steps + 1))
    # first 20 steps are fixed at dt = 0.01 yr
    assert np.allclose(np.diff(short_run.step_times_years)[:20], 0.01)


def test_run_stops_at_the_minimum_temperature(background):
    run = run_cooling(
        background, replace(CoolingConfig(), minimum_temperature_K=3e8, final_years=1e5)
    )
    assert run.stop_reason == "minimum_local_temperature"
    assert run.step_times_years[-1] < 1e5


@pytest.fixture(scope="module")
def long_run(background):
    return run_cooling(background, replace(CoolingConfig(), final_years=1e5))


def test_cooling_is_monotonic_with_the_default_start(long_run):
    """Guards the backward-Euler start-up: without it the outer crust rings and re-heats."""
    u = np.array([long_run.profiles[n] for n in sorted(long_run.profiles)])
    assert not (np.diff(u, axis=0) > 0).any()
    assert np.all(np.diff(long_run.luminosity_erg_s) < 0)


def test_the_last_step_is_not_shortened(long_run):
    """The run ends at the first step that reaches the final time, without trimming it."""
    t = long_run.step_times_years
    assert t[-2] < 1e5 <= t[-1]
    dt = np.diff(t)
    assert dt[-1] >= dt[-2] - 1e-9 * t[-1]  # the same step as before, or the grown one


def test_plain_backward_euler_is_available(background):
    cfg = replace(CoolingConfig(), final_years=100.0, theta=1.0, predictor_corrector=False)
    default = run_cooling(background, replace(CoolingConfig(), final_years=100.0))
    euler = run_cooling(background, cfg)
    assert euler.completed_steps != default.completed_steps or not np.allclose(
        euler.luminosity_erg_s[-1], default.luminosity_erg_s[-1]
    )
    assert np.isfinite(euler.luminosity_erg_s).all()


def test_rejected_steps_leave_the_state_untouched(background):
    """A tight limit on the relative change forces a rejection; accepted steps must respect it."""
    cfg = replace(
        CoolingConfig(),
        final_years=2.0,
        max_change=0.02,
        grow_below_change=0.005,
        target_change=0.01,
    )
    run = run_cooling(background, cfg)
    assert run.rejected_steps > 0
    assert np.all(run.relative_changes[20:] <= cfg.max_change * (1 + 1e-12))
