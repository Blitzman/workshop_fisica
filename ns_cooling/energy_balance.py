"""Energy-balance check of a cooling run:  dE/dt = -(L_nu + L_gamma).

The thermal energy at infinity, the neutrino luminosity and the photon
luminosity are computed from the temperature profiles of every accepted step,
with exactly the registered microphysics of the run (so any new heat-capacity
term, neutrino process or envelope is checked automatically):

    E       = sum_i V_i e^{Phi_i} U_i(T_i),   U(T) = integral of C_v dT from 0 to T
    L_nu    = sum_i V_i e^{2 Phi_i} Q_i(T_i)
    L_gamma = 4 pi R^2 sigma T_s^4 e^{2 Phi}

with V the proper volume of the cell (CGS) and T the local temperature. The
ledger integrals use the trapezoid rule over the accepted steps. ``C_v T`` is not
the energy: U is tabulated per cell from the heat capacity itself.

The closure error of a time window is (dE + int L dt) / (-dE): the fraction of
the energy lost in the window that the ledger does not account for. A
positive value means that the star loses less energy than the luminosity
integral says (the usual lag of first-order time stepping).

The unresolved start (first ~0.2 yr from 1e10 K, where the neutrino cooling time
is far shorter than the first step) never closes; the report therefore starts
at ``start_years``.
"""

from dataclasses import dataclass

import numpy as np

from cooling import surface_luminosity
from emissivity import neutrino_emissivity_total
from heat_capacity import heat_capacity
from state import LocalState
from thermal_geometry import thermal_geometry
from thermal_matrix import SECONDS_PER_YEAR

KM3_TO_CM3 = 1e15


@dataclass
class EnergyBudget:
    """Energy, luminosities and cumulative radiated energy at every accepted step."""

    times_years: np.ndarray
    thermal_energy: np.ndarray  # E [erg]
    neutrino_luminosity: np.ndarray  # L_nu [erg/s]
    photon_luminosity: np.ndarray  # L_gamma [erg/s]
    neutrino_energy: np.ndarray  # integral of L_nu dt from t = 0 [erg]
    photon_energy: np.ndarray  # integral of L_gamma dt from t = 0 [erg]

    @property
    def residual(self):
        """E + integral of (L_nu + L_gamma) dt - E(0): zero if the energy is conserved."""
        return (
            self.thermal_energy + self.neutrino_energy + self.photon_energy - self.thermal_energy[0]
        )

    def as_columns(self):
        names = (
            "time_years",
            "thermal_energy_erg",
            "L_neutrino_erg_s",
            "L_photon_erg_s",
            "radiated_neutrino_erg",
            "radiated_photon_erg",
            "residual_erg",
        )
        values = (
            self.times_years,
            self.thermal_energy,
            self.neutrino_luminosity,
            self.photon_luminosity,
            self.neutrino_energy,
            self.photon_energy,
            self.residual,
        )
        return names, np.column_stack(values)


@dataclass
class ClosureWindow:
    start_years: float
    end_years: float
    energy_lost: float  # -dE over the window [erg]
    closure_error: float  # (dE + int L dt) / (-dE)


@dataclass
class EnergyClosure:
    windows: list
    start_years: float
    tolerance: float

    @property
    def max_error(self):
        return max((abs(w.closure_error) for w in self.windows), default=float("nan"))

    @property
    def ok(self):
        return bool(self.windows) and self.max_error <= self.tolerance

    def format(self):
        lines = [
            f"Energy balance dE/dt = -(L_nu + L_gamma), windows from t >= {self.start_years:g} yr",
            "  window [yr]            energy lost [erg]   closure error",
        ]
        for w in self.windows:
            lines.append(
                f"  {w.start_years:9.3g} - {w.end_years:9.3g}   {w.energy_lost:16.3e}   "
                f"{100 * w.closure_error:+10.3f}%"
            )
        verdict = "OK" if self.ok else "WARNING: above the tolerance"
        lines.append(
            f"  largest closure error {100 * self.max_error:.3f}%  (tolerance {100 * self.tolerance:g}%): {verdict}"
        )
        return "\n".join(lines)


# t_min_K, t_max_K and points: free numerical choices (range and size of the thermal-energy tables).
def thermal_energy_tables(background, *, t_min_K=1e4, t_max_K=4e10, points=4000, terms=None):
    """Per-cell table of U(T) = integral of C_v dT [erg cm^-3] on a logarithmic temperature grid.

    The (tiny) energy below ``t_min_K`` is approximated by C_v(t_min) t_min / 2.
    Returns (ln T grid, U of shape (points, n_cells)).
    """
    ln_t = np.linspace(np.log(t_min_K), np.log(t_max_K), points)
    grid = np.exp(ln_t)[:, None]
    cv = heat_capacity(LocalState(grid, background.cells), terms=terms)["total"]
    integrand = cv * grid  # integral in d ln T
    energy = np.zeros_like(cv)
    energy[1:] = np.cumsum(0.5 * (integrand[1:] + integrand[:-1]) * np.diff(ln_t)[:, None], axis=0)
    energy[0] = integrand[0] / 2.0
    energy[1:] += energy[0]
    return ln_t, energy


def energy_budget(run, background, *, chunk=500):
    """Energy, L_nu and L_gamma at every accepted step of a ``cooling.CoolingRun``.

    The run must have stored the profile of every step (the default).
    Saved Q and Cv callable snapshots are used without changing the active registries.
    """
    if run.neutrino_processes is None or run.heat_capacity_terms is None:
        raise ValueError("The energy budget needs the microphysics snapshot saved by run_cooling")
    steps = list(range(run.completed_steps + 1))
    if any(n not in run.profiles for n in steps):
        raise ValueError("The energy budget needs the profile of every accepted step")
    geometry = thermal_geometry(background)
    profiles = np.array([run.profiles[n] for n in steps])  # T_inf / 1e8 K, (steps, cells)
    phi = background.cells.ephi
    t_local = profiles * 1e8 / phi
    ln_t, energy = thermal_energy_tables(
        background, t_min_K=min(1e4, 0.5 * t_local.min()), t_max_K=max(4e10, 2.0 * t_local.max()),
        terms=run.heat_capacity_terms,
    )
    volume = geometry.volumes * KM3_TO_CM3  # proper cell volumes [cm^3]
    ln_u = np.log(energy)
    u_cells = np.empty_like(t_local)
    for i in range(t_local.shape[1]):
        u_cells[:, i] = np.exp(np.interp(np.log(t_local[:, i]), ln_t, ln_u[:, i]))
    thermal = (u_cells * (volume * phi)).sum(axis=1)
    neutrino = np.empty(len(steps))
    for a in range(0, len(steps), chunk):
        q = neutrino_emissivity_total(
            LocalState(t_local[a : a + chunk], background.cells), processes=run.neutrino_processes
        )
        neutrino[a : a + chunk] = (q * (volume * phi**2)).sum(axis=1)
    photon = np.array([surface_luminosity(p[-1], background, geometry) for p in profiles])
    times = np.asarray(run.step_times_years)[: len(steps)]
    dt = np.diff(times) * SECONDS_PER_YEAR

    def cumulative(values):
        return np.concatenate([[0.0], np.cumsum(0.5 * (values[1:] + values[:-1]) * dt)])

    return EnergyBudget(times, thermal, neutrino, photon, cumulative(neutrino), cumulative(photon))


def validate_energy_options(start_years, tolerance):
    """Reject invalid diagnostic options before integration or window iteration."""
    if not np.isfinite(start_years) or start_years <= 0:
        raise ValueError("energy start_years must be positive and finite")
    if not np.isfinite(tolerance) or tolerance <= 0:
        raise ValueError("energy tolerance must be positive and finite")


# start_years = 1 yr and tolerance = 2 % are free choices.
def energy_closure(budget, *, start_years=1.0, tolerance=0.02):
    """Closure error per decade of time, from ``start_years`` to the end of the run.

    Window edges are accepted steps, so no interpolation is involved.
    """
    validate_energy_options(start_years, tolerance)
    t, energy = budget.times_years, budget.thermal_energy
    luminosity = budget.neutrino_luminosity + budget.photon_luminosity
    windows = []
    edge = start_years
    while edge < t[-1]:
        a = int(np.searchsorted(t, edge))
        b = int(np.searchsorted(t, min(10.0 * edge, t[-1])))
        b = min(b, len(t) - 1)
        lost = energy[a] - energy[b]
        if b > a and lost > 0:
            # Integrate locally: subtracting cumulative residuals can lose
            # precision because the unresolved first step dominates the ledger.
            radiated = np.sum(
                0.5
                * (luminosity[a:b] + luminosity[a + 1 : b + 1])
                * np.diff(t[a : b + 1])
                * SECONDS_PER_YEAR
            )
            windows.append(ClosureWindow(t[a], t[b], lost, (radiated - lost) / lost))
        edge *= 10.0
    return EnergyClosure(windows, start_years, tolerance)
