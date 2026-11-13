"""Frozen-coefficient theta-method system for the thermal evolution.

theta = 1 is backward Euler (fully implicit), theta = 1/2 is Crank-Nicolson.

This module assembles the system; it does not solve, advance time, iterate,
apply a temperature floor.
"""

from dataclasses import dataclass
import numpy as np
from conductivity import thermal_transport
from heat_capacity import heat_capacity_internal
from emissivity import linearized_neutrino_source
from state import LocalState
from thermal_geometry import thermal_geometry, cell_temperatures, interior_face_temperatures
from thermal_boundary import surface_boundary

SECONDS_PER_YEAR = 3.15576e7  # Julian year, 365.25 d (decision D-007)


@dataclass(frozen=True)
class ThermalSystem:
    lower: np.ndarray
    diagonal: np.ndarray
    upper: np.ndarray
    rhs: np.ndarray
    heat_capacity: np.ndarray
    source_offset: np.ndarray
    source_slope: np.ndarray
    face_transport_internal: np.ndarray
    face_conductance: np.ndarray
    boundary: object
    geometry: object
    dt_seconds: float
    theta: float = 1.0

    def matvec(self, values):
        v = np.asarray(values, dtype=float)
        if v.shape != self.diagonal.shape:
            raise ValueError("Vector length differs from number of cells")
        result = self.diagonal * v
        result[1:] += self.lower[1:] * v[:-1]
        result[:-1] += self.upper[:-1] * v[1:]
        return result


def interior_face_transport(Tinf_8, background, geometry):
    """Conductivity (internal units) at the interior faces, evaluated once each."""
    local_8 = interior_face_temperatures(Tinf_8, geometry)
    # Local temperature and critical temperatures are both in Kelvin.
    state = LocalState(local_8 * 1e8, background.faces.interior())
    return thermal_transport(state)["transport_internal"]


def assemble_matrix(Tinf_8, dt_years, geometry, cv, q0, q1, face_k, boundary, theta=1.0):
    """Assemble with supplied frozen microphysics; all arrays use internal units.

    With N the (backward-Euler) operator and s the constant source terms, the
    theta-method solves  (I + theta*N) u_new = (I - (1-theta)*N) u_old - s.
    theta = 1 gives backward Euler; theta = 1/2 gives Crank-Nicolson, whose
    explicit half uses the loss evaluated exactly at the old state.
    """
    if not 0.5 <= theta <= 1.0:
        raise ValueError("theta must be between 1/2 (Crank-Nicolson) and 1 (backward Euler)")
    u = cell_temperatures(Tinf_8, geometry)
    if not np.isfinite(dt_years) or dt_years <= 0:
        raise ValueError("dt_years must be positive and finite")
    cv, q0, q1 = [np.asarray(a, dtype=float) for a in (cv, q0, q1)]
    if any(a.shape != u.shape or not np.isfinite(a).all() for a in (cv, q0, q1)):
        raise ValueError("Cell coefficients must be finite and match the cell vector")
    if np.any(cv <= 0):
        raise ValueError("Heat capacity must be positive")
    k = np.asarray(face_k, dtype=float)
    if k.shape != (len(u) - 1,) or not np.isfinite(k).all() or np.any(k < 0):
        raise ValueError("One nonnegative finite conductivity per interior face required")
    if not np.isfinite([boundary.flux_offset, boundary.flux_slope]).all():
        raise ValueError("Non-finite surface coefficients")
    dt = dt_years * SECONDS_PER_YEAR
    metric_k = k * geometry.ephi_faces[1:-1] / geometry.elam_faces[1:-1]
    conductance = metric_k * geometry.areas[1:-1] / geometry.center_distances
    factor = dt / cv / geometry.volumes
    tf = theta * factor  # implicit weight; for theta = 1 this is exactly factor
    lower = -tf * np.r_[0.0, conductance]
    upper = -tf * np.r_[conductance, 0.0]
    diagonal = 1 - (lower + upper) + tf * geometry.volumes * q1
    rhs = u - factor * geometry.volumes * q0
    # Center: zero area, hence zero flux. Surface: linearized radiative flux.
    diagonal[-1] += tf[-1] * geometry.areas[-1] * boundary.flux_slope
    rhs[-1] -= factor[-1] * geometry.areas[-1] * boundary.flux_offset
    if theta < 1.0:
        # Explicit part: -(1-theta) * N u_old, with N = (implicit operator) / theta.
        n_u = (diagonal - 1) * u
        n_u[1:] += lower[1:] * u[:-1]
        n_u[:-1] += upper[:-1] * u[1:]
        rhs -= ((1.0 - theta) / theta) * n_u
    return ThermalSystem(
        lower, diagonal, upper, rhs, cv, q0, q1, k, conductance, boundary, geometry, dt, theta
    )


@dataclass(frozen=True)
class FrozenCoefficients:
    """Microphysics evaluated at ONE state, in internal units.

    The neutrino source and the surface flux are linearized about that same
    state (offset + slope * u), with T/Tc held fixed in the source derivative.
    """

    heat_capacity: np.ndarray
    source_offset: np.ndarray
    source_slope: np.ndarray
    face_transport: np.ndarray
    boundary: object


def evaluate_coefficients(Tinf_8, background, *, geometry=None):
    """Evaluate the microphysics at the state Tinf_8 (one value per physical cell).

    Conductivity: degenerate electrons in the crust, npe matter in the core.
    Heat capacity and neutrino emissivity include the fixed-T/Tc source
    linearization.
    """
    if geometry is None:
        geometry = thermal_geometry(background)
    u = cell_temperatures(Tinf_8, geometry)
    cells = LocalState(u * 1e8 / background.cells.ephi, background.cells)
    cv = heat_capacity_internal(cells)
    q0, q1 = linearized_neutrino_source(cells)
    k = interior_face_transport(u, background, geometry)
    boundary = surface_boundary(
        u[-1],
        background.g14,
        geometry.ephi_faces[-1],
        background.envelope_model,
        luminosity_factor=background.surface_factor,
    )
    return FrozenCoefficients(cv, q0, q1, k, boundary)


def build_thermal_system(
    Tinf_8, dt_years, background, *, geometry=None, theta=1.0, coefficients=None
):
    """Assemble the theta-method system for a step starting from the state Tinf_8.

    ``coefficients`` default to the microphysics evaluated at that same (old)
    state. A predictor-corrector passes coefficients evaluated at another
    state (the midpoint of the step) while the right-hand side still starts
    from Tinf_8.
    """
    if geometry is None:
        geometry = thermal_geometry(background)
    if coefficients is None:
        coefficients = evaluate_coefficients(Tinf_8, background, geometry=geometry)
    c = coefficients
    u = cell_temperatures(Tinf_8, geometry)
    return assemble_matrix(
        u,
        dt_years,
        geometry,
        c.heat_capacity,
        c.source_offset,
        c.source_slope,
        c.face_transport,
        c.boundary,
        theta,
    )
