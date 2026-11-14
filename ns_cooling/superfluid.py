"""Pairing gap models and superfluid functions."""

import numpy as np
from constants import BOLTZMANN_MEV_PER_K
from numerical_validation import checked

# Gap models 1-16 (1-7 crust neutron 1S0, 8-10 proton 1S0, 11-16 core neutron 3P2) in the four-parameter form of Ho et al.
# (2015), Eq. 2. The numbering follows NSCool (README, to be checked); the source of each row is not recorded. Checked against the Ho et
# al. models (item SF-2): model 11 equals EEHO (triplet) after the conversion of the normalisation; models 11 and 15, and 2 and 7,
# are identical; the other rows do not match any Ho et al. model within 3 %.
# imod: (Delta_scale [MeV], k1 [fm^-1], k2 [fm^-2], k3 [fm^-1], k4 [fm^-2])
GAP_MODELS = {
    1: (9.3, 0.02, 0.6, 1.55, 0.32),
    2: (68.0, 0.1, 4.0, 1.7, 4.0),
    3: (4.0, 0.4, 1.5, 1.65, 0.05),
    4: (22.0, 0.3, 0.09, 1.05, 4.0),
    5: (2.9, 0.3, 0.017, 1.3, 0.07),
    6: (82.11, 0.05, 2.31, 1.63, 5.81),
    7: (68.0, 0.1, 4.0, 1.7, 4.0),
    8: (120.0, 0.0, 9.0, 1.3, 1.8),
    9: (55.0, 0.15, 4.0, 1.27, 4.0),
    10: (2.27, 0.1, 0.07, 1.05, 0.25),
    11: (0.068, 1.28, 0.1, 2.37, 0.02),
    12: (0.15, 2.0, 0.1, 3.1, 0.02),
    13: (2.2, 1.05, 1.0, 2.82, 0.6),
    14: (0.425, 1.1, 0.5, 2.7, 0.5),
    15: (0.068, 1.28, 0.1, 2.37, 0.02),
    16: (2.9, 1.21, 0.5, 1.62, 0.5),
}
# k_B T_c / Delta: 0.5669 (singlet) and 0.8416 (triplet) for the Delta of this table; the Ho et al. table uses 0.1187 for its
# triplet Delta0 (item SF-3). Type 3 is not used.
TC_GAP_RATIO = {1: 0.5669, 2: 0.8416, 3: 0.4926}


def fermi_quantities(rho, particle_fraction):
    """Return the particle number density [cm^-3] and the Fermi momentum kF [fm^-1]."""
    rho, particle_fraction = np.broadcast_arrays(
        np.asarray(rho, float), np.asarray(particle_fraction, float)
    )
    if not (np.isfinite(rho).all() and np.isfinite(particle_fraction).all()):
        raise ValueError("rho and particle_fraction must be finite")
    if np.any(rho < 0) or np.any((particle_fraction < 0) | (particle_fraction > 1)):
        raise ValueError("rho must be non-negative and particle_fraction in [0, 1]")
    with np.errstate(over="raise", invalid="raise"):
        number_density = rho / 1.66e-24 * particle_fraction
        kF = (3.0 * np.pi**2 * number_density) ** (1.0 / 3.0) * 1e-13
    return number_density, kF


def pairing_gap_from_kF(gap_model, kF):
    """Zero-temperature gap Delta(T=0) in MeV, with each model's own normalization."""
    if gap_model not in GAP_MODELS:
        raise ValueError("Unknown gap model")
    scale, k1, k2, k3, k4 = GAP_MODELS[gap_model]
    kF = np.asarray(kF, float)
    if not np.isfinite(kF).all() or np.any(kF < 0):
        raise ValueError("kF must be finite and non-negative")
    if scale < 0 or k2 <= 0 or k4 <= 0 or k1 >= k3:
        raise ValueError("Invalid gap-model parameters")
    gap_MeV = np.zeros_like(kF)
    mask = (kF > k1) & (kF <= k3)
    left, right = (kF[mask] - k1) ** 2, (kF[mask] - k3) ** 2
    gap_MeV[mask] = scale * (left / (left + k2)) * (right / (right + k4))
    return gap_MeV


def pairing_gap(gap_model, pairing_type, rho, particle_fraction):
    """Return (Delta0 [MeV], Tc [K]) for a gap model and pairing channel."""
    if pairing_type not in TC_GAP_RATIO:
        raise ValueError("Unknown pairing type")
    _, kF = fermi_quantities(rho, particle_fraction)
    gap_MeV = pairing_gap_from_kF(gap_model, kF)
    return gap_MeV, TC_GAP_RATIO[pairing_type] * gap_MeV / BOLTZMANN_MEV_PER_K


def pairing_critical_temperatures(cfg, star):
    """Critical temperatures [K] at the faces of crust neutrons, core neutrons and core protons.

    ``cfg`` holds the on/off switches and gap models of the input file; ``star``
    is a dictionary with the arrays ``rho``, ``yn``, ``yp`` and ``xh`` at the faces. Returns ``tc_crust``, ``tc_n`` and
    ``tc_p`` (zero outside the region of each species or when it is switched off).
    """
    result = {}
    for key, name, pairing_type, fraction, region in (
        ("crust_n", "tc_crust", 1, star["yn"], star["xh"] > 0.0),
        ("core_n", "tc_n", 2, star["yn"], star["xh"] == 0.0),
        ("core_p", "tc_p", 1, star["yp"], star["xh"] == 0.0),
    ):
        _, Tc_K = pairing_gap(cfg[key]["model"], pairing_type, star["rho"], fraction)
        result[name] = np.where(region & cfg[key]["enabled"], Tc_K, 0.0)
    return result


def dimensionless_pairing_gap(tau, pairing_type):
    """Dimensionless gap v for pairing types 1 and 2, valid for 0 <= tau <= 1.

    For tau <= 0.1 a constant 1e-10 is added to tau, so that tau = 0 is finite. Outside
    this domain the normal phase must be used instead of this function.
    """
    tau = checked(tau, "tau")
    if np.any(tau > 1) or pairing_type not in (1, 2):
        raise ValueError("tau <= 1 and pairing_type 1 or 2 required")
    alpha = 0.5669 if pairing_type == 1 else 0.8416
    v = np.empty_like(tau)
    low = tau <= 0.1
    v[low] = 1.0 / alpha / (tau[low] + 1e-10)
    x = tau[~low]
    if pairing_type == 1:
        v[~low] = np.sqrt(1 - x) * (1.456 - 0.157 / np.sqrt(x) + 1.764 / x)
    else:
        v[~low] = np.sqrt(1 - x) * (0.7893 + 1.188 / x)
    return v
