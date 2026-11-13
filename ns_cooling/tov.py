"""Relativistic hydrostatic structure of a cold star (Tolman-Oppenheimer-Volkoff).

The integration runs once, from the center to the outer edge of the equation of
state, in the variable x = ln(P/Pc), with the radius r, the gravitational mass
m (geometric units, G m / c^2 in cm) and the metric function Phi as unknowns.
The solution keeps its dense output, so any quantity can be evaluated afterwards
at any pressure, so a mesh of any kind can be built afterwards, independently of the integration.

The start uses the regular expansion at a small radius (``center_radius_cm``)
to avoid the apparent singularity at r = 0. The metric function is normalized
globally by matching exp(Phi) to the Schwarzschild exterior value at the outer
edge of the table. That edge is the first row of the EOS, not a computed
photosphere (the layers above it are neglected).

An equation of state is any object with ``pressure_min``, ``pressure_max`` and
``energy_density(P)`` [g cm^-3] (see ``eos_table.CompOSETable``).
"""

from dataclasses import dataclass

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq, minimize_scalar

from constants import C, G, M_SUN, PRESSURE_TO_GEOM, RHO_TO_GEOM

_GRID_POINTS = 200001  # sampling used to invert r(x)


def _rhs(x, y, eos, pc):
    r, m, _ = y
    pressure = pc * np.exp(x)
    lo, hi = eos.pressure_min, eos.pressure_max
    # Guard end-point rounding only (a few ulps); never extrapolate the EOS.
    if lo * (1 - 1e-13) <= pressure < lo:
        pressure = lo
    if hi < pressure <= hi * (1 + 1e-13):
        pressure = hi
    rho = float(eos.energy_density(pressure)) * RHO_TO_GEOM
    p = pressure * PRESSURE_TO_GEOM
    if r <= 0 or m < 0 or r <= 2 * m:
        raise ValueError("Invalid TOV state or trapped surface")
    dp_dr = -(rho + p) * (m + 4 * np.pi * r**3 * p) / (r * (r - 2 * m))
    dr_dx = p / dp_dr
    return np.array([dr_dx, 4 * np.pi * r * r * rho * dr_dx, -p / (rho + p)])


class StarSolution:
    """The structure of one star: global quantities and the dense solution.

    Attributes: pressure_central, pressure_outer [dyn cm^-2]; mass_solar [Msun] and
    radius_cm of the outer edge; surface_ephi = exp(Phi) there and
    surface_gravity_cm_s2 = local proper acceleration there (with the pressure term).
    """

    def __init__(self, eos, pc, po, y0, x0, ode, shift, rtol, r0):
        self._ode, self._x0, self._shift, self._r0 = ode, x0, shift, r0
        self._rho_c = float(eos.energy_density(pc)) * RHO_TO_GEOM
        self._p_c = pc * PRESSURE_TO_GEOM
        self.pressure_central, self.pressure_outer = pc, po
        self.x_outer = float(np.log(po / pc))
        self.rtol, self.center_radius_cm = rtol, r0
        edge = self.profile(self.x_outer)
        self.radius_cm = float(edge["r_cm"][0])
        self.mass_solar = float(edge["mass_g"][0]) / M_SUN
        self.surface_ephi = float(edge["ephi"][0])
        self.surface_gravity_cm_s2 = float(edge["g_cm_s2"][0])
        self._grid = None

    def profile(self, x):
        """Structure at x = ln(P/Pc) (scalar or array, 0 at the center, negative outwards).

        Returns a dict of arrays: r_cm, mass_g, ephi, elam, g_cm_s2 (proper acceleration
        of a static observer, including the pressure term).
        """
        x = np.atleast_1d(np.asarray(x, dtype=float))
        if np.any(x > 1e-12) or np.any(x < self.x_outer * (1 + 1e-12) - 1e-12):
            raise ValueError("x = ln(P/Pc) outside the integrated range")
        x = np.clip(x, self.x_outer, 0.0)
        r, m, phi = np.empty_like(x), np.empty_like(x), np.empty_like(x)
        seed = x > self._x0  # between the center and the first integrated point
        ratio = np.sqrt(np.clip(x[seed] / self._x0, 0.0, 1.0))
        r[seed] = self._r0 * ratio
        m[seed] = 4 * np.pi * self._rho_c * r[seed] ** 3 / 3
        phi[seed] = 2 * np.pi * (self._rho_c / 3 + self._p_c) * r[seed] ** 2
        if (~seed).any():
            r[~seed], m[~seed], phi[~seed] = self._ode.sol(x[~seed])
        pressure_geom = self.pressure_central * np.exp(x) * PRESSURE_TO_GEOM
        safe = np.where(r > 0, r, 1.0)
        elam = np.where(r > 0, 1.0 / np.sqrt(1.0 - 2.0 * m / safe), 1.0)
        g = np.where(r > 0, C**2 * (m + 4 * np.pi * safe**3 * pressure_geom) / safe**2 * elam, 0.0)
        return dict(
            r_cm=r, mass_g=m * C**2 / G, ephi=np.exp(phi + self._shift), elam=elam, g_cm_s2=g
        )

    def _radius_grid(self):
        if self._grid is None:
            t = np.linspace(0.0, 1.0, _GRID_POINTS)
            x = self.x_outer * t**4  # clustered at the center, where r(x) is steepest
            self._grid = (x, self.profile(x)["r_cm"])
        return self._grid

    def x_at_radius(self, r_cm):
        """x = ln(P/Pc) at a radius [cm] (inverse of the monotonic r(x))."""
        x, r = self._radius_grid()
        return np.interp(r_cm, r, x)

    def pressure_at_x(self, x):
        return self.pressure_central * np.exp(x)


def integrate_star(
    eos, pressure_central, *, pressure_outer=None, center_radius_cm=10.0, rtol=1e-9, max_step=0.05
):
    """Integrate the structure from the center (Pc) to the outer pressure (default: the EOS minimum)."""
    pc = float(pressure_central)
    if eos.pressure_max < pc <= eos.pressure_max * (1 + 1e-12):  # rounding of a log10 grid at the table edge
        pc = eos.pressure_max
    po = eos.pressure_min if pressure_outer is None else float(pressure_outer)
    if not eos.pressure_min <= po < pc <= eos.pressure_max:
        raise ValueError("Require EOS Pmin <= Pouter < Pcentral <= EOS Pmax")
    values = np.array([center_radius_cm, rtol, max_step], dtype=float)
    if not np.all(np.isfinite(values)) or np.any(values <= 0):
        raise ValueError("Positive finite numerical settings required")
    rho_c = float(eos.energy_density(pc)) * RHO_TO_GEOM
    p_c = pc * PRESSURE_TO_GEOM
    r0 = float(center_radius_cm)
    # Regular expansion at the center to the leading nonzero order; Phi(0) = 0 gauge.
    delta = 2 * np.pi * (rho_c + p_c) * (rho_c / 3 + p_c) * r0 * r0
    if delta >= p_c * 1e-4:
        raise ValueError("Center radius too large for the leading-order expansion")
    p_start = (p_c - delta) / PRESSURE_TO_GEOM
    if p_start <= po:
        raise ValueError("The center seed lies beyond the outer pressure")
    y0 = np.array([r0, 4 * np.pi * rho_c * r0**3 / 3, 2 * np.pi * (rho_c / 3 + p_c) * r0 * r0])
    x0, x1 = np.log(p_start / pc), np.log(po / pc)
    ode = solve_ivp(
        _rhs,
        (x0, x1),
        y0,
        args=(eos, pc),
        method="DOP853",
        rtol=rtol,
        atol=[1e-5, 1e-6, 1e-12],
        max_step=max_step,
        first_step=min(max_step, abs(x1 - x0), max(1e-16, 0.1 * abs(x0))),
        dense_output=True,
    )
    if not ode.success:
        raise RuntimeError(ode.message)
    r_o, m_o, phi_o = ode.y[:, -1]
    shift = 0.5 * np.log1p(-2 * m_o / r_o) - phi_o  # exp(Phi) = Schwarzschild value at the edge
    return StarSolution(eos, pc, po, y0, x0, ode, shift, rtol, r0)


@dataclass(frozen=True)
class StableBranch:
    """The stable branch of the mass-central pressure relation of an EOS.

    M(Pc) has a minimum and a maximum; the stars between them (M increasing with Pc) are
    stable. Below the minimum lie configurations with a huge, tenuous envelope (for a table
    with a near-plateau of P at the crust-core transition, radii of thousands of km), above the
    maximum the stars collapse.
    """

    pressure_min_mass: float
    mass_min: float
    pressure_max_mass: float
    mass_max: float


def _mass_scan(eos, rtol, max_step, points):
    lower = (
        eos.pressure_core
        if hasattr(eos, "pressure_core") and np.isfinite(eos.pressure_core)
        else None
    )
    if lower is None:
        lower = max(eos.pressure_min * 100.0, 1e30)
    grid = np.linspace(np.log10(lower), np.log10(eos.pressure_max), points)

    def mass(log_pc):
        return integrate_star(eos, 10.0**log_pc, rtol=rtol, max_step=max_step).mass_solar

    return grid, np.array([mass(g) for g in grid]), mass


def _refine_extremum(mass, grid, index, sign, hi_limit):
    """log10 Pc and mass of the extremum of M near grid[index] (sign +1: minimum, -1: maximum)."""
    bounds = (grid[max(index - 1, 0)], grid[min(index + 1, hi_limit)])
    best = minimize_scalar(
        lambda g: sign * mass(g), bounds=bounds, method="bounded", options=dict(xatol=1e-4)
    )
    return float(best.x), float(sign * best.fun)


def stable_branch(eos, *, rtol=1e-8, max_step=0.1, points=16):
    """Pressures and masses of the minimum and maximum mass of the stable branch (a ``StableBranch``).

    A scan of ``points`` central pressures from the first core row of the table (or a low
    pressure if it has none) to its maximum locates the extrema, which are then refined.
    """
    grid, masses, mass = _mass_scan(eos, rtol, max_step, points)
    peak = int(np.argmax(masses))
    trough = int(np.argmin(masses[: peak + 1]))
    g_max, m_max = _refine_extremum(mass, grid, peak, -1.0, len(grid) - 1)
    g_min, m_min = _refine_extremum(mass, grid, trough, +1.0, peak)
    return StableBranch(10.0**g_min, m_min, 10.0**g_max, m_max)


def central_pressure_for_mass(eos, mass_solar, *, rtol=1e-8, max_step=0.1, xtol=1e-9):
    """Central pressure [dyn cm^-2] of the star of the given mass on the stable branch.

    On the stable branch the mass grows with Pc between the minimum and the maximum mass
    (see ``StableBranch``). A coarse scan brackets the target inside that branch, the
    extremum is refined only when the target lies beyond the scanned masses, and Brent's method
    refines log10 Pc to ``xtol``. The searches use looser integration settings than the final
    structure. Targets outside [M_min, M_max] raise ValueError.
    """
    if not (np.isfinite(mass_solar) and mass_solar > 0):
        raise ValueError("The target mass must be positive and finite")
    grid, masses, mass = _mass_scan(eos, rtol, max_step, 14)
    peak = int(np.argmax(masses))
    trough = int(np.argmin(masses[: peak + 1]))
    above = np.nonzero(masses[trough : peak + 1] >= mass_solar)[0]
    if above.size and above[0] > 0:  # the target is crossed inside the scanned stable branch
        i = trough + int(above[0])
        lo, hi = grid[i - 1], grid[i]
    elif above.size:  # at or below the lightest scanned stable star: find the true minimum mass
        g_min, m_min = _refine_extremum(mass, grid, trough, +1.0, peak)
        if mass_solar < m_min:
            raise ValueError(
                f"Target mass {mass_solar:g} Msun is below the minimum mass {m_min:.4f} Msun of the stable branch"
            )
        lo, hi = g_min, grid[trough]
    else:  # above every scanned stable star: find the true maximum mass
        g_max, m_max = _refine_extremum(mass, grid, peak, -1.0, len(grid) - 1)
        if mass_solar > m_max:
            raise ValueError(
                f"Target mass {mass_solar:g} Msun exceeds the maximum mass {m_max:.4f} Msun"
            )
        lo, hi = grid[max(peak - 1, trough)], g_max
    return 10.0 ** brentq(lambda g: mass(g) - mass_solar, lo, hi, xtol=xtol)
