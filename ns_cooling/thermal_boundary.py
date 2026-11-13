"""Nonmagnetic envelope with consistent exponent and flux derivative.

Gudmundsson (corrected), Fe or accreted envelopes; local Kelvin.
Flux coefficients retain the evolution equation's redshift convention.
"""

from dataclasses import dataclass
import numpy as np
from thermal_fluxes import thermal_fluxes


@dataclass(frozen=True)
class SurfaceBoundary:
    base_temperature_8: float
    surface_temperature_K: float
    surface_temperature_infinity_K: float
    flux_local_internal: float
    flux_offset: float
    flux_slope: float


def surface_boundary(
    Tinf_outer_8, g14, ephi_outer, envelope_model="gudmundsson", luminosity_factor=1.0
):
    if not all(np.isfinite(x) and x > 0 for x in (Tinf_outer_8, g14, ephi_outer)):
        raise ValueError("Boundary temperature, gravity and redshift must be positive and finite")
    x = float(Tinf_outer_8 / ephi_outer)
    result = thermal_fluxes(x * 1e8, g14, envelope_model)
    ts = float(result.surface_temperature_K)
    return SurfaceBoundary(
        x,
        ts,
        ts * ephi_outer,
        float(result.flux_internal),
        float(result.flux_offset) * ephi_outer**2 * luminosity_factor,
        float(result.flux_slope) * ephi_outer * luminosity_factor,
    )
