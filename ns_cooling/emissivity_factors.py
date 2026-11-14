"""Scalar pairing-reduction and PBF factors for the neutrino processes.

The function names identify the process and the paired species. These are NOT the heat-capacity reduction factors, which
belong to a different module (heat_capacity.py).
"""

import numpy as np
from numpy import sqrt, exp


# The fits of this module are those of Yakovlev et al. (2001), Eqs. (199)-(229) (README: Yakovlev and Levenfish 1995, Levenfish and
# Yakovlev 1994, believed). They equal the transcription of the GitHub code to 1e-15, except the asymptotic form for v > 15 of the
# neutron branch of the modified Urca (item NU-11).
def exp_with_floor(z):
    """exp(z) with the argument floored at -150.

    Used only where the formula applies this cutoff explicitly.
    """
    return exp(np.maximum(z, -150.0))


def direct_urca_singlet_factor(v):
    """Dimensionless pairing-reduction factor; v = Delta(T)/(k_B T)."""
    z = 3.427 - sqrt(3.427**2 + v**2)
    return (0.2312 + sqrt(0.7688**2 + (0.1438 * v) ** 2)) ** 5.5 * exp_with_floor(z)


def direct_urca_triplet_m0_factor(v):
    """Dimensionless pairing-reduction factor; v = Delta(T)/(k_B T)."""
    z = 2.701 - sqrt(2.701**2 + v**2)
    return (0.2546 + sqrt(0.7454**2 + (0.1284 * v) ** 2)) ** 5.0 * exp_with_floor(z)


def direct_urca_triplet_m2_factor(v):
    """Dimensionless pairing-reduction factor; v = Delta(T)/(k_B T)."""
    z = 1.0 - sqrt(1.0 + (0.4129 * v) ** 2)
    return (0.5 + (0.09226 * v) ** 2) / (
        1.0 + (0.1821 * v) ** 2 + (0.16736 * v) ** 4
    ) + 0.5 * exp_with_floor(z)


def direct_urca_combined_factor(vn, vp):
    """Dimensionless pairing-reduction factor for neutron and proton gaps (vn, vp)."""
    vn, vp = np.broadcast_arrays(np.asarray(vn, float), np.asarray(vp, float))
    out = np.empty_like(vn)
    small = vn**2 + vp**2 <= 25.0
    n2, p2 = vn[small] ** 2, vp[small] ** 2
    n4, p4 = vn[small] ** 4, vp[small] ** 4
    out[small] = (1e4 - 2.839 * p4 - 5.022 * n4) / (
        1e4 + 757.0 * p2 + 1494.0 * n2 + 211.1 * p2 * n2 + 0.4832 * p4 * n4
    )
    out[~small] = np.minimum(
        direct_urca_triplet_m0_factor(vn[~small]), direct_urca_singlet_factor(vp[~small])
    )
    return out


def modified_urca_n_proton_factor(v):
    """Dimensionless pairing-reduction factor; v = Delta(T)/(k_B T)."""
    a = 0.1477 + sqrt(0.8523**2 + (0.1175 * v) ** 2)
    b = 0.1477 + sqrt(0.8523**2 + (0.1297 * v) ** 2)
    return 0.5 * (a**7.5 + b**5.5) * exp_with_floor(3.4370 - sqrt(3.4370**2 + v**2))


def modified_urca_p_proton_factor(v):
    """Dimensionless pairing-reduction factor; v = Delta(T)/(k_B T)."""
    c = 0.2414 + sqrt(0.7586**2 + (0.1318 * v) ** 2)
    # Unlike the other factors, this one deliberately has no -150 exponent floor.
    return c**7 * exp(5.339 - sqrt(5.339**2 + (2.0 * v) ** 2))


def modified_urca_n_neutron_singlet_factor(v):
    """Dimensionless pairing-reduction factor; v = Delta(T)/(k_B T)."""
    return modified_urca_p_proton_factor(v)


def modified_urca_p_neutron_singlet_factor(v):
    """Dimensionless pairing-reduction factor; v = Delta(T)/(k_B T)."""
    return modified_urca_n_proton_factor(v)


def modified_urca_n_neutron_triplet_factor(v):
    """Dimensionless pairing-reduction factor; v = Delta(T)/(k_B T)."""
    v = np.asarray(v, float)
    out = np.empty_like(v)
    # For v > 15 an asymptotic form is used, whose source is not recorded; it falls a factor ~52 below the Yakovlev et al. (2001) fit
    # at v = 15 (item NU-11, open).
    large = v > 15.0
    prefactor = 1.20960 * 2.0 * 0.209 / 0.11513 / np.pi**8 / 3.0 / np.sqrt(3.0)
    out[large] = prefactor * v[large] ** 6 * exp_with_floor(-2.0 * v[large])
    out[~large] = modified_urca_p_proton_factor(v[~large])
    return out


def modified_urca_p_neutron_triplet_factor(v):
    """Dimensionless pairing-reduction factor; v = Delta(T)/(k_B T)."""
    a = 0.1612 + sqrt(0.8388**2 + (0.1117 * v) ** 2)
    b = 0.1612 + sqrt(0.8388**2 + (0.1274 * v) ** 2)
    return 0.5 * (a**7 + b**5) * exp_with_floor(2.398 - sqrt(2.398**2 + v**2))


def modified_urca_n_combined_factor(vn, vp):
    """Dimensionless pairing-reduction factor for neutron and proton gaps (vn, vp)."""
    return (
        direct_urca_combined_factor(2.0 * vn, vp)
        * modified_urca_n_proton_factor(vp)
        / direct_urca_singlet_factor(vp)
    )


def modified_urca_p_combined_factor(vn, vp):
    """Dimensionless pairing-reduction factor for neutron and proton gaps (vn, vp)."""
    return (
        direct_urca_combined_factor(vn, 2.0 * vp)
        * modified_urca_p_neutron_triplet_factor(vn)
        / direct_urca_triplet_m0_factor(vn)
    )


def brems_np_proton_factor(v):
    """Dimensionless pairing-reduction factor; v = Delta(T)/(k_B T)."""
    a = 0.9982 + sqrt(0.0018**2 + (0.3815 * v) ** 2)
    b = 0.3949 + sqrt(0.6051**2 + (0.2666 * v) ** 2)
    z1 = 1.306 - sqrt(1.306**2 + v**2)
    z2 = 3.303 - sqrt(3.303**2 + 4.0 * v**2)
    return (a * exp_with_floor(z1) + 1.732 * b**7 * exp_with_floor(z2)) / 2.732


def brems_pp_proton_factor(v):
    """Dimensionless pairing-reduction factor; v = Delta(T)/(k_B T)."""
    c = 0.1747 + sqrt(0.8253**2 + (0.07933 * v) ** 2)
    d = 0.7333 + sqrt(0.2667**2 + (0.1678 * v) ** 2)
    z1 = 4.228 - sqrt(4.228**2 + 4.0 * v**2)
    z2 = 7.762 - sqrt(7.762**2 + 9.0 * v**2)
    return (c**2 * exp_with_floor(z1) + d**7.5 * exp_with_floor(z2)) / 2.0


def brems_np_neutron_singlet_factor(v):
    """Dimensionless pairing-reduction factor; v = Delta(T)/(k_B T)."""
    return brems_np_proton_factor(v)


def brems_np_neutron_triplet_factor(v):
    """Dimensionless pairing-reduction factor; v = Delta(T)/(k_B T)."""
    return brems_np_neutron_singlet_factor(
        v
    )  # Same expression as the corresponding singlet/proton channel, as in the model.


def brems_nn_neutron_singlet_factor(v):
    """Dimensionless pairing-reduction factor; v = Delta(T)/(k_B T)."""
    return brems_pp_proton_factor(v)


def brems_nn_neutron_triplet_factor(v):
    """Dimensionless pairing-reduction factor; v = Delta(T)/(k_B T)."""
    return brems_pp_proton_factor(
        v
    )  # Same expression as the corresponding singlet/proton channel, as in the model.


def brems_np_combined_factor(vn, vp):
    """Dimensionless pairing-reduction factor for neutron and proton gaps (vn, vp)."""
    return (
        direct_urca_combined_factor(vn, vp)
        * brems_np_proton_factor(vp)
        / direct_urca_singlet_factor(vp)
    )


def pbf_singlet_function(v):
    """Dimensionless PBF (pair breaking and formation) function; v = Delta(T)/(k_B T)."""
    a = 0.602 * v**2 + 0.5942 * v**4 + 0.288 * v**6
    b = 0.5547 + sqrt(0.4453**2 + 0.0113 * v**2)
    return a * sqrt(b) * exp_with_floor(2.245 - sqrt(2.245**2 + 4.0 * v**2))


def pbf_triplet_m0_function(v):
    """Dimensionless PBF (pair breaking and formation) function; v = Delta(T)/(k_B T)."""
    a = 1.204 * v**2 + 3.733 * v**4 + 0.3191 * v**6
    b = 0.7591 + sqrt(0.2409**2 + 0.3145 * v**2)
    c = 1.0 + 0.3511 * v**2
    return a * b**2 / c * exp_with_floor(0.4616 - sqrt(0.4616**2 + 4.0 * v**2))
