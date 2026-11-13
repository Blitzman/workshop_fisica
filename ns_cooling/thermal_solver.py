"""One time step of the thermal evolution (frozen-coefficient theta-method, optionally with a
predictor-corrector); no time loop."""

import numpy as np
from scipy.linalg import solve_banded
from thermal_geometry import thermal_geometry
from thermal_matrix import build_thermal_system, evaluate_coefficients

TEMPERATURE_FLOOR_8 = 1e-5


def solve_thermal_system(system):
    """Return the raw solution in redshifted 1e8 K, without clipping.

    Diagonals follow ThermalSystem's row-index convention. Inputs are not
    overwritten. SciPy checks finite coefficients and reports singularity.
    """
    lower, diagonal, upper, rhs = [
        np.asarray(a, dtype=float)
        for a in (system.lower, system.diagonal, system.upper, system.rhs)
    ]
    if diagonal.ndim != 1 or diagonal.size == 0:
        raise ValueError("Expected a nonempty one-dimensional diagonal")
    if any(a.shape != diagonal.shape for a in (lower, upper, rhs)):
        raise ValueError("All four arrays must have the same shape")
    if lower[0] != 0 or upper[-1] != 0:
        raise ValueError("Nonzero coefficients outside the physical domain")
    bands = np.zeros((3, diagonal.size))
    bands[0, 1:] = upper[:-1]
    bands[1] = diagonal
    bands[2, :-1] = lower[1:]
    result = solve_banded((1, 1), bands, rhs, check_finite=True)
    if not np.isfinite(result).all():
        raise FloatingPointError("Non-finite thermal solution")
    return result


def thermal_step(
    Tinf_8,
    dt_years,
    background,
    *,
    geometry=None,
    theta=0.5,
    predictor_corrector=True,
    coefficients=None,
):
    """Advance one step from the state Tinf_8; apply the temperature floor.

    Single stage (predictor_corrector=False): coefficients frozen at the OLD
    state, theta-method with theta = 1 (backward Euler) or 1/2 (Crank-Nicolson).

    Predictor-corrector (default): (1) a backward-Euler predictor with coefficients at
    the old state gives a first estimate of the new state; (2) the coefficients
    are re-evaluated at the MIDPOINT of the old and predicted states, and the
    theta-method corrector (normally theta = 1/2) is solved from the old state.
    The loss, heat capacity and conductivity are then centered in time, which
    is what makes Crank-Nicolson second order for the nonlinear problem. Cost:
    two microphysics evaluations per step.

    ``coefficients`` optionally supplies the microphysics already evaluated at
    the old state (reused when a rejected step is repeated with a smaller dt).
    Returns (new_temperature, system); the system is the one of the last
    solve. No clock/output update.
    """
    if geometry is None:
        geometry = thermal_geometry(background)
    if coefficients is None:
        coefficients = evaluate_coefficients(Tinf_8, background, geometry=geometry)
    if not predictor_corrector:
        system = build_thermal_system(
            Tinf_8, dt_years, background, geometry=geometry, theta=theta, coefficients=coefficients
        )
        return np.maximum(solve_thermal_system(system), TEMPERATURE_FLOOR_8), system
    predictor = build_thermal_system(
        Tinf_8, dt_years, background, geometry=geometry, theta=1.0, coefficients=coefficients
    )
    predicted = np.maximum(solve_thermal_system(predictor), TEMPERATURE_FLOOR_8)
    midpoint = evaluate_coefficients(0.5 * (Tinf_8 + predicted), background, geometry=geometry)
    system = build_thermal_system(
        Tinf_8, dt_years, background, geometry=geometry, theta=theta, coefficients=midpoint
    )
    return np.maximum(solve_thermal_system(system), TEMPERATURE_FLOOR_8), system
