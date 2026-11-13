"""Numerical settings of the cooling run (time stepping and stop conditions)."""

from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True)
class CoolingConfig:
    """Settings of the time integration.

    Time stepping
        The first ``n_fixed_steps`` accepted steps use ``dt_initial_years``
        and are never rejected. From the next step on, a step is rejected and
        repeated with a smaller dt if the largest relative temperature change
        between cells exceeds ``max_change``. If it is below ``grow_below_change``,
        the next dt is increased (by at most ``max_growth_factor``) aiming at
        ``target_change``. Otherwise dt is kept. Before late_growth_start_years, dt never exceeds
        ``dt_max_years``. Afterwards each accepted step proposes 1.05 times its dt,
        without that cap; the same temperature-change rejection remains active; going below ``dt_min_years`` is an error.

    Time integration scheme
        ``theta`` selects the theta-method: 0.5 is Crank-Nicolson (default),
        1.0 is backward Euler. Coefficients are frozen at the old state unless
        the predictor-corrector below is used.

    Predictor-corrector
        With ``predictor_corrector=True`` (default) each step first predicts the
        new state (backward Euler, coefficients at the old state), then
        re-evaluates the coefficients at the midpoint and solves the theta-method
        corrector (``theta`` = 0.5 is the intended choice). Two microphysics
        evaluations per step instead of one. Plain backward
        Euler is ``theta=1.0, predictor_corrector=False``.

    Start-up damping
        With ``fixed_steps_backward_euler=True`` (default) the first ``n_fixed_steps`` steps
        use plain backward Euler whatever ``theta`` and ``predictor_corrector``
        say. Backward Euler damps completely the fast modes that the
        unresolved first steps (dt = 0.01 yr from 1e10 K) would otherwise seed
        and that Crank-Nicolson does not remove.

    Initial state
        ``initial_temperature_8`` is the uniform initial T_infinity in units of
        1e8 K (100 = 1e10 K). The input file provides it to the driver.

    Stop conditions
        Default: stop at the first accepted state with L_gamma_infinity < 1e31 erg/s.
        The optional time and temperature limits below are disabled by default.
        The run ends at the first accepted step that reaches or passes ``final_years``
        (the last step is not shortened, so the final time can exceed ``final_years`` by
        up to one step) or when the minimum local temperature drops to
        ``minimum_temperature_K``.

    Output
        Luminosity samples are recorded only for times later than
        ``luminosity_start_years`` (the initial state has no luminosity sample).
    """

    initial_temperature_8: float = 100.0
    final_years: float = float("inf")
    minimum_temperature_K: float = 0.0
    minimum_luminosity_erg_s: float = 1e31
    late_growth_start_years: float = 1e5
    late_growth_factor: float = 1.05
    dt_initial_years: float = 0.01
    dt_min_years: float = 1e-12
    dt_max_years: float = 1000.0
    n_fixed_steps: int = 20
    max_change: float = 0.05
    grow_below_change: float = 0.01
    target_change: float = 0.025
    min_step_factor: float = 0.2
    max_shrink_factor: float = 0.5
    max_growth_factor: float = 2.0
    luminosity_start_years: float = 0.01
    theta: float = 0.5
    predictor_corrector: bool = True
    fixed_steps_backward_euler: bool = True

    def __post_init__(self):
        if not np.isfinite(self.initial_temperature_8) or self.initial_temperature_8 <= 0:
            raise ValueError("initial_temperature_8 must be positive and finite")
        if not np.isfinite(self.minimum_temperature_K) or self.minimum_temperature_K < 0:
            raise ValueError("minimum_temperature_K must be nonnegative and finite (0 disables it)")
        if np.isnan(self.final_years) or self.final_years <= self.luminosity_start_years:
            raise ValueError(
                "final_years must be greater than luminosity_start_years (infinity disables it) "
                f"({self.luminosity_start_years:g})"
            )
        if not np.isfinite(self.minimum_luminosity_erg_s) or self.minimum_luminosity_erg_s <= 0:
            raise ValueError("minimum_luminosity_erg_s must be positive and finite")
        if not np.isfinite(self.late_growth_start_years) or self.late_growth_start_years < 0:
            raise ValueError("late_growth_start_years must be nonnegative and finite")
        if not np.isfinite(self.late_growth_factor) or self.late_growth_factor < 1:
            raise ValueError("late_growth_factor must be finite and >= 1")
        if not 0.5 <= self.theta <= 1.0:
            raise ValueError("theta must be between 0.5 (Crank-Nicolson) and 1.0 (backward Euler)")
        if self.n_fixed_steps < 0:
            raise ValueError("n_fixed_steps must be non-negative")
