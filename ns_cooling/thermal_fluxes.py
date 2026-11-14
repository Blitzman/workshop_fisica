"""Standalone B=0 envelope (Tb-Ts) relations.

Inputs Tb [local K], g14=g/(1e14 cm s^-2). No stellar model is needed.
Fe and accreted fits: Potekhin et al. (2003), ApJ 594, 404. 'accreted' is that
paper's accreted-envelope fit, not an arbitrary pure-H envelope model.
"""

from dataclasses import dataclass
import numpy as np

from constants import SIGMA_SB

MODELS = ("gudmundsson", "fe", "accreted")
SIGMA = SIGMA_SB  # Stefan-Boltzmann constant [erg cm^-2 s^-1 K^-4]


@dataclass(frozen=True)
class ThermalFluxes:
    surface_temperature_K: np.ndarray
    flux_cgs: np.ndarray
    flux_internal: np.ndarray
    flux_offset: np.ndarray
    flux_slope: np.ndarray


def thermal_fluxes(Tb_K, g14=1.36, model="gudmundsson"):
    """Return Ts and local flux F ≈ offset + slope*Tb8 about the input Tb.

    Internal F = CGS F * 1e-30; slope is dF_internal/dTb8.
    Slopes are analytic for each model, not the Gudmundsson slope reused
    for Fe/accreted. Inputs may be scalars or broadcastable arrays.
    """
    T, g = np.broadcast_arrays(np.asarray(Tb_K, float), np.asarray(g14, float))
    if np.any(~np.isfinite(T)) or np.any(T <= 0) or np.any(~np.isfinite(g)) or np.any(g <= 0):
        raise ValueError("Tb and g14 must be positive and finite")
    if model not in MODELS:
        raise ValueError(f"model must be one of {MODELS}")
    x = T / 1e8
    if model == "gudmundsson":
        # Gudmundsson-type relation with exponent 1/1.82.
        ts = 1e6 * g**0.25 * (x / 1.288) ** (1 / 1.82)
        flux = SIGMA * ts**4 * 1e-30
        slope = 4 / 1.82 * flux / x
    else:
        t = T / 1e9
        z = t - 1e-3 * g**0.25 * np.sqrt(7 * t)
        if np.any(z <= 0):
            raise ValueError("Fe/accreted fit requires zeta>0; increase Tb")
        dz = 1 - 1e-3 * g**0.25 * np.sqrt(7) / (2 * np.sqrt(t))
        fe = g * ((7 * z) ** 2.25 + (z / 3) ** 1.25)
        dfe = g * (2.25 * 7 * (7 * z) ** 1.25 + 1.25 / 3 * (z / 3) ** 0.25) * dz
        if model == "fe":
            q, dq = fe, dfe
        else:
            h = (6.2 * t) ** 4
            corr = 0.447 + 0.075 * np.log10(1e9 * t) / (1 + h)
            dcorr = (
                0.075 * ((1 + h) / (t * np.log(10)) - np.log10(1e9 * t) * 4 * h / t) / (1 + h) ** 2
            )
            a = g * (18.1 * t) ** 2.42 * corr
            da = g * (18.1 * t) ** 2.42 * (2.42 * corr / t + dcorr)
            b = 3.2 * t**1.67
            db = 1.67 * b / t
            q = (a + b * fe) / (1 + b)
            dq = (da + db * fe + b * dfe - q * db) / (1 + b)
        ts = 1e6 * q**0.25
        flux = SIGMA * 1e24 * q * 1e-30
        slope = SIGMA * 1e24 * dq * 0.1 * 1e-30
    if np.any(~np.isfinite(ts)) or np.any(~np.isfinite(slope)):
        raise ValueError("Nonfinite envelope result")
    return ThermalFluxes(ts, flux * 1e30, flux, flux - slope * x, slope)
