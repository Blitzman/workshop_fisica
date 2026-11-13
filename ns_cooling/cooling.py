"""Time integration of the cooling: adaptive time steps and luminosity at infinity."""

from dataclasses import dataclass

import numpy as np

from config import CoolingConfig
from state import LocalState
from emissivity import NEUTRINO_PROCESSES
from heat_capacity import HEAT_CAPACITY_TERMS
from thermal_boundary import surface_boundary
from thermal_fluxes import SIGMA
from thermal_geometry import interior_face_temperatures, thermal_geometry
from thermal_matrix import evaluate_coefficients
from thermal_solver import thermal_step


def surface_luminosity(u_outer, background, geometry):
    """Luminosity at infinity [erg/s] for the outermost cell temperature ``u_outer`` (T_inf/1e8 K).

    L = 4 pi R^2 sigma T_s^4 exp(2 Phi), with T_s from the envelope relation (R and Phi of the
    last face, times ``background.surface_factor`` to refer them to the outer edge).
    """
    surface = surface_boundary(
        u_outer, background.g14, geometry.ephi_faces[-1], background.envelope_model
    )
    radius_cm = geometry.radius_faces[-1] * 1e5
    ts = surface.surface_temperature_K
    base = 4.0 * np.pi * radius_cm**2 * SIGMA * ts**4 * geometry.ephi_faces[-1] ** 2
    return base * background.surface_factor  # 1 unless the outer edge lies above the last face


@dataclass
class CoolingRun:
    """Result of a cooling run."""

    times_years: np.ndarray  # times of the luminosity samples [yr]
    luminosity_erg_s: np.ndarray  # luminosity at infinity [erg/s]
    profiles: dict  # accepted step -> T_inf/1e8 K in every cell
    completed_steps: int
    stop_reason: str  # "minimum_luminosity", "final_time" or "minimum_local_temperature"
    step_times_years: np.ndarray  # time after every accepted step, starting at t=0
    relative_changes: np.ndarray  # largest relative T change of every accepted step
    rejected_steps: int
    neutrino_processes: tuple | None = None  # ordered callable snapshot for diagnostics
    heat_capacity_terms: tuple | None = None
    neutrino_luminosities: np.ndarray | None = None  # all step_times, including t=0
    dominant_neutrino_process: np.ndarray | None = None


def run_cooling(background, config=None, *, snapshots=None):
    """Evolve the star from a uniform initial temperature.

    ``snapshots`` is an optional collection of step numbers whose temperature
    profiles are stored (default: all steps). The last accepted state is
    always stored.
    """
    neutrino_processes = NEUTRINO_PROCESSES.items()
    heat_capacity_terms = HEAT_CAPACITY_TERMS.items()
    cfg = config if config is not None else CoolingConfig()
    dt = cfg.dt_initial_years
    time = 0.0
    step_times = [0.0]
    changes = []
    rejected = 0
    geo = thermal_geometry(background)
    u = np.full(geo.n_cells, cfg.initial_temperature_8)  # T_inf / 1e8 K
    weights = geo.volumes * 1e15 * background.cells.ephi**2
    def neutrino_monitor(temperature):
        state = LocalState(temperature * 1e8 / background.cells.ephi, background.cells)
        values = np.array([np.sum(fn(state) * weights) for _, fn in neutrino_processes])
        if not np.all(np.isfinite(values)) or np.any(values < 0):
            raise RuntimeError("Invalid integrated neutrino luminosity")
        name = neutrino_processes[int(np.argmax(values))][0] if values.size and values.max() > 0 else "none"
        return values, name
    first_lnu, first_dominant = neutrino_monitor(u)
    neutrino_history = [first_lnu]
    dominant_history = [first_dominant]
    profiles = {0: u.copy()} if snapshots is None or 0 in snapshots else {}
    times = []
    luminosities = []
    completed = 0
    reason = "final_time"
    coefficients = None  # microphysics at the current accepted state; reused if a step is retried
    while time < cfg.final_years:
        if surface_luminosity(u[-1], background, geo) < cfg.minimum_luminosity_erg_s:
            reason = "minimum_luminosity"
            break
        n = completed + 1
        trial_dt = dt
        if time + trial_dt == time:
            raise RuntimeError("Time step cannot advance floating-point time")
        local_min = (
            min(
                np.min(u / background.cells.ephi),
                np.min(interior_face_temperatures(u, geo)),
            )
            * 1e8
        )
        if cfg.minimum_temperature_K > 0 and local_min <= cfg.minimum_temperature_K:
            reason = "minimum_local_temperature"
            print(
                f"Stop: min local T={local_min:.6g} K <= {cfg.minimum_temperature_K:g} K",
                flush=True,
            )
            break
        if coefficients is None:
            coefficients = evaluate_coefficients(u, background, geometry=geo)
        startup = cfg.fixed_steps_backward_euler and n <= cfg.n_fixed_steps
        trial, _ = thermal_step(
            u,
            float(trial_dt),
            background,
            geometry=geo,
            theta=1.0 if startup else cfg.theta,
            predictor_corrector=cfg.predictor_corrector and not startup,
            coefficients=coefficients,
        )
        change = float(np.max(np.abs((trial - u) / u)))
        if not np.isfinite(change):
            raise RuntimeError("Non-finite temperature change")
        # Steps 1..n_fixed_steps are never rejected; step n_fixed_steps+1 is
        # the first adaptive one.
        if n > cfg.n_fixed_steps and change > cfg.max_change:
            rejected += 1
            factor = float(
                np.clip(cfg.target_change / change, cfg.min_step_factor, cfg.max_shrink_factor)
            )
            dt = trial_dt * factor
            if dt < cfg.dt_min_years:
                raise RuntimeError(f"Adaptive time step below {cfg.dt_min_years:g} yr")
            continue  # Neither time nor accepted state/output has changed.
        u = trial
        coefficients = None
        time += trial_dt
        step_times.append(time)
        changes.append(change)
        # The decision taken after step n_fixed_steps sets dt of the first adaptive step.
        if time >= cfg.late_growth_start_years and n >= cfg.n_fixed_steps:
            dt = trial_dt * cfg.late_growth_factor
        elif n >= cfg.n_fixed_steps and change < cfg.grow_below_change:
            if change == 0:
                factor = cfg.max_growth_factor
            else:
                factor = float(
                    np.clip(
                        cfg.target_change / change,
                        cfg.min_step_factor,
                        cfg.max_growth_factor,
                    )
                )
            dt = min(cfg.dt_max_years, trial_dt * factor)
        else:
            dt = trial_dt
        # Output belongs to the updated state; the matrix was assembled with the old one.
        if time > cfg.luminosity_start_years:
            lum = surface_luminosity(u[-1], background, geo)
            times.append(time)
            luminosities.append(lum)
        lnu, dominant = neutrino_monitor(u)
        neutrino_history.append(lnu)
        dominant_history.append(dominant)
        completed = n
        if snapshots is None or n in snapshots:
            profiles[n] = u.copy()
        if n % 100 == 0 or time >= cfg.final_years:
            print(f"Step {n}, t={time:.6g} yr, dt={trial_dt:.6g} yr", flush=True)
    profiles[completed] = u.copy()
    return CoolingRun(
        np.asarray(times),
        np.asarray(luminosities),
        profiles,
        completed,
        reason,
        np.asarray(step_times),
        np.asarray(changes),
        rejected,
        neutrino_processes,
        heat_capacity_terms,
        np.asarray(neutrino_history),
        np.asarray(dominant_history),
    )
