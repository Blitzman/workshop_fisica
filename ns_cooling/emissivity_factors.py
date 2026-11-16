"""Scalar pairing-reduction and PBF factors for the neutrino processes.

The function names identify the process and the paired species. These are NOT the heat-capacity reduction factors, which
belong to a different module (heat_capacity.py). v = Delta(T) / (k_B T) throughout.

Sources (decision D-034; the equation numbers were checked in the papers by the literature audit of 15 November 2026):
- Direct Urca, modified Urca, bremsstrahlung and pair-breaking functions: Yakovlev, Kaminker, Gnedin and Haensel (2001), Phys. Rep.
  354, 1, Eqs. (199)-(229) and (236)-(241) (Y01). The two-species factors use their approximations (204)-(205), (225)-(227), (229).
- Neutron branch of the modified Urca with neutron triplet pairing and normal protons: the near-origin and intermediate expressions of
  Gusakov (2002), A&A 389, 702, Table 2 and Eqs. (A1), (A2), (A7), (A10), on the neutron-only axis, joined to the asymptote of
  Y01 Eq. (216) with xi = 0.130 (Yakovlev and Levenfish 1995, Eqs. 26 and 29). The joins (v = 3-4 and 8-20) are cubic interpolations
  in log R; no independent accuracy is claimed inside them. The Y01 fit that the code used below v = 15 is a proton-singlet
  surrogate meant for moderate pairing (v up to about 10); it is up to 53 times larger than this factor at v = 15.
- Two-species direct-Urca factor: Eq. (205) adjusted by the ratios of the exact single-species fits to its axis values; beyond radius 5
  the strong-pairing approximation (204) is matched at r = 5 along each ray, with the matching offset in log R decaying as
  exp[-(r-5)/5]. The prescription is local, continuous in value and recovers (204); its two-species accuracy outside the original fit
  domain is not quantified.
All the factors are evaluated in log space and have no exponent floor, so that strongly suppressed rates tend to zero.
Known limit (open item NU-10): the neutron-only factor of the modified Urca is not continuous with the combined factor at the
temperature where the protons become superfluid (the combined factor, at v_p -> 0, is 1.4 times the neutron-only one at v_n = 2, 4.7 times
at 5 and 17 times at 15).
"""

import numpy as np
from numpy import sqrt, exp


# Unless otherwise indicated, the fits are from Yakovlev et al. (2001).


def exp_with_floor(z):
    """exp(max(z, -150)). Kept only for the plasmon-decay rate of emissivity.py; no pairing factor uses it."""
    return exp(np.maximum(z, -150.0))


def _log_direct_urca_singlet(v):
    z = 3.427 - np.hypot(3.427, v)
    return 5.5 * np.log(0.2312 + np.hypot(0.7688, 0.1438 * v)) + z


def _log_direct_urca_triplet_m0(v):
    z = 2.701 - np.hypot(2.701, v)
    return 5.0 * np.log(0.2546 + np.hypot(0.7454, 0.1284 * v)) + z


def direct_urca_singlet_factor(v):
    """Dimensionless pairing-reduction factor; v = Delta(T)/(k_B T)."""
    return exp(_log_direct_urca_singlet(v))


def direct_urca_triplet_m0_factor(v):
    """Dimensionless pairing-reduction factor; v = Delta(T)/(k_B T)."""
    return exp(_log_direct_urca_triplet_m0(v))


def direct_urca_triplet_m2_factor(v):
    """Dimensionless pairing-reduction factor; v = Delta(T)/(k_B T)."""
    z = 1.0 - sqrt(1.0 + (0.4129 * v) ** 2)
    return (0.5 + (0.09226 * v) ** 2) / (
        1.0 + (0.1821 * v) ** 2 + (0.16736 * v) ** 4
    ) + 0.5 * exp(z)


def _log_direct_urca_ba_weak(vn, vp):
    """Log of Eq. (205), normalized to the exact single-species fits."""
    n2, p2 = vn**2, vp**2
    n4, p4 = n2**2, p2**2
    fit = (1e4 - 2.839 * p4 - 5.022 * n4) / (
        1e4 + 757.0 * p2 + 1494.0 * n2 + 211.1 * p2 * n2 + 0.4832 * p4 * n4
    )
    n_fit = (1e4 - 5.022 * n4) / (1e4 + 1494.0 * n2)
    p_fit = (1e4 - 2.839 * p4) / (1e4 + 757.0 * p2)
    log_axis_correction = (
        _log_direct_urca_triplet_m0(vn) + _log_direct_urca_singlet(vp) - np.log(n_fit * p_fit)
    )
    return np.log(fit) + log_axis_correction


def _log_direct_urca_combined(vn, vp):
    vn, vp = np.broadcast_arrays(np.asarray(vn, float), np.asarray(vp, float))
    radius = np.hypot(vn, vp)
    out = np.empty_like(radius)
    weak = radius <= 5.0
    out[weak] = _log_direct_urca_ba_weak(vn[weak], vp[weak])

    # Match the strong-pairing approximation (204) to the fit (205) at
    # radius 5 along each ray. Relax the log offset over five gap units.
    strong = ~weak
    boundary_scale = 5.0 / radius[strong]
    n_edge = vn[strong] * boundary_scale
    p_edge = vp[strong] * boundary_scale
    edge_min = np.minimum(_log_direct_urca_triplet_m0(n_edge), _log_direct_urca_singlet(p_edge))
    edge_offset = _log_direct_urca_ba_weak(n_edge, p_edge) - edge_min
    strong_min = np.minimum(
        _log_direct_urca_triplet_m0(vn[strong]), _log_direct_urca_singlet(vp[strong])
    )
    out[strong] = strong_min + edge_offset * exp(-(radius[strong] - 5.0) / 5.0)

    # The exact single-species expressions apply on the axes.
    return np.where(
        vp == 0, _log_direct_urca_triplet_m0(vn),
        np.where(vn == 0, _log_direct_urca_singlet(vp), out)
    )


def direct_urca_combined_factor(vn, vp):
    """Dimensionless pairing-reduction factor for neutron and proton gaps."""
    return exp(_log_direct_urca_combined(vn, vp))


def _log_modified_urca_n_proton(v):
    a = 0.1477 + sqrt(0.8523**2 + (0.1175 * v) ** 2)
    b = 0.1477 + sqrt(0.8523**2 + (0.1297 * v) ** 2)
    return np.log(0.5 * (a**7.5 + b**5.5)) + 3.4370 - np.hypot(3.4370, v)


def modified_urca_n_proton_factor(v):
    """Dimensionless pairing-reduction factor; v = Delta(T)/(k_B T)."""
    return exp(_log_modified_urca_n_proton(v))


def modified_urca_p_proton_factor(v):
    """Dimensionless pairing-reduction factor; v = Delta(T)/(k_B T)."""
    c = 0.2414 + sqrt(0.7586**2 + (0.1318 * v) ** 2)
    return exp(7.0 * np.log(c) + 5.339 - np.hypot(5.339, 2.0 * v))


def modified_urca_n_neutron_singlet_factor(v):
    """Dimensionless pairing-reduction factor; v = Delta(T)/(k_B T)."""
    return modified_urca_p_proton_factor(v)


def modified_urca_p_neutron_singlet_factor(v):
    """Dimensionless pairing-reduction factor; v = Delta(T)/(k_B T)."""
    return modified_urca_n_proton_factor(v)


def modified_urca_n_neutron_triplet_factor(v):
    """Neutron modified Urca branch with type-B neutrons and normal protons."""
    v = np.asarray(v, float)
    v2 = v**2

    # Gusakov (2002), Table 2, Eqs. (A2), (A10), evaluated at vp=0.
    log_weak = -v2 * (0.565001 + 1.667194e-4 * v2**2) / sqrt(
        1.0 + v2 * (0.173165 + v2**2 * (7.710124e-8 + 3.782805e-6 * v2))
    )

    # Eqs. (A1), (A7) on the same axis (phi=pi/2, t=0).
    z = np.sin(0.275434)**2
    y = np.cos(0.330574)**2
    a = -0.719681 + 1.260056
    b = 0.148464
    c = 0.375796 + 0.096843 / (y * (1.0 + 3.100942 * z**2)**3)
    log_intermediate = -a * v2 / (1.0 + b * v2)**c

    # Blend the two fits over 3 <= v <= 4 near Gusakov's region boundary.
    weight = np.clip(v - 3.0, 0.0, 1.0)
    weight = weight**2 * (3.0 - 2.0 * weight)
    log_factor = (1.0 - weight) * log_weak + weight * log_intermediate

    # Beyond the range of stated fit accuracy use the Y01 Eq. (216)
    # asymptote, with xi=0.130 from Yakovlev & Levenfish (1995).
    # The transition from v=8 to 20 is an interpolation in log R.
    tail = v > 8.0
    vt = v[tail]
    prefactor = (120960.0 / (11513.0 * np.pi**8)) * (2.0 * 0.130 / (3.0 * sqrt(3.0)))
    log_asymptote = np.log(prefactor) + 6.0 * np.log(vt) - 2.0 * vt
    weight = np.clip((vt - 8.0) / 12.0, 0.0, 1.0)
    weight = weight**2 * (3.0 - 2.0 * weight)
    out = np.empty_like(v)
    out[~tail] = log_factor[~tail]
    out[tail] = (1.0 - weight) * log_intermediate[tail] + weight * log_asymptote
    return exp(out)


def _log_modified_urca_p_neutron_triplet(v):
    a = 0.1612 + sqrt(0.8388**2 + (0.1117 * v) ** 2)
    b = 0.1612 + sqrt(0.8388**2 + (0.1274 * v) ** 2)
    return np.log(0.5 * (a**7 + b**5)) + 2.398 - np.hypot(2.398, v)


def modified_urca_p_neutron_triplet_factor(v):
    """Dimensionless pairing-reduction factor; v = Delta(T)/(k_B T)."""
    return exp(_log_modified_urca_p_neutron_triplet(v))


def modified_urca_n_combined_factor(vn, vp):
    """Dimensionless pairing-reduction factor for neutron and proton gaps (vn, vp)."""
    log_ratio = _log_direct_urca_combined(2.0 * vn, vp) - _log_direct_urca_singlet(vp)
    return exp(log_ratio + _log_modified_urca_n_proton(vp))


def modified_urca_p_combined_factor(vn, vp):
    """Dimensionless pairing-reduction factor for neutron and proton gaps (vn, vp)."""
    log_ratio = _log_direct_urca_combined(vn, 2.0 * vp) - _log_direct_urca_triplet_m0(vn)
    return exp(log_ratio + _log_modified_urca_p_neutron_triplet(vn))


def _log_brems_np_proton(v):
    a = 0.9982 + sqrt(0.0018**2 + (0.3815 * v) ** 2)
    b = 0.3949 + sqrt(0.6051**2 + (0.2666 * v) ** 2)
    z1 = 1.306 - sqrt(1.306**2 + v**2)
    z2 = 3.303 - sqrt(3.303**2 + 4.0 * v**2)
    return np.logaddexp(np.log(a) + z1, np.log(1.732) + 7.0 * np.log(b) + z2) - np.log(2.732)


def brems_np_proton_factor(v):
    """Dimensionless pairing-reduction factor; v = Delta(T)/(k_B T)."""
    return exp(_log_brems_np_proton(v))


def brems_pp_proton_factor(v):
    """Dimensionless pairing-reduction factor; v = Delta(T)/(k_B T)."""
    c = 0.1747 + sqrt(0.8253**2 + (0.07933 * v) ** 2)
    d = 0.7333 + sqrt(0.2667**2 + (0.1678 * v) ** 2)
    z1 = 4.228 - sqrt(4.228**2 + 4.0 * v**2)
    z2 = 7.762 - sqrt(7.762**2 + 9.0 * v**2)
    return (exp(2.0 * np.log(c) + z1) + exp(7.5 * np.log(d) + z2)) / 2.0


def brems_np_neutron_singlet_factor(v):
    """Dimensionless pairing-reduction factor; v = Delta(T)/(k_B T)."""
    return brems_np_proton_factor(v)


def brems_np_neutron_triplet_factor(v):
    """Dimensionless pairing-reduction factor; v = Delta(T)/(k_B T)."""
    return brems_np_neutron_singlet_factor(
        v
    )  # Same expression as the singlet/proton channel (Y01 Eq. 220), as in the model. The v_p -> 0 limit of the combined Eq. (229) is the direct-Urca type-B factor, 0.39 times this at v_n = 5 (audit); the measured effect on L is below 0.01 % and it is kept (D-035).


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
    log_ratio = _log_direct_urca_combined(vn, vp) - _log_direct_urca_singlet(vp)
    return exp(log_ratio + _log_brems_np_proton(vp))


def pbf_singlet_function(v):
    """Dimensionless PBF (pair breaking and formation) function; v = Delta(T)/(k_B T)."""
    v2 = v**2
    a = 0.602 + 0.5942 * v2 + 0.288 * v2**2
    b = 0.5547 + sqrt(0.4453**2 + 0.0113 * v2)
    log_v2 = np.log(np.where(v2 > 0.0, v2, 1.0))
    log_factor = log_v2 + np.log(a) + 0.5 * np.log(b) + 2.245 - np.hypot(2.245, 2.0 * v)
    return np.where(v2 == 0.0, 0.0, exp(log_factor))


def pbf_triplet_m0_function(v):
    """Dimensionless PBF (pair breaking and formation) function; v = Delta(T)/(k_B T)."""
    v2 = v**2
    a = 1.204 + 3.733 * v2 + 0.3191 * v2**2
    b = 0.7591 + sqrt(0.2409**2 + 0.3145 * v2)
    c = 1.0 + 0.3511 * v2
    log_v2 = np.log(np.where(v2 > 0.0, v2, 1.0))
    log_factor = log_v2 + np.log(a) + 2.0 * np.log(b) - np.log(c) + 0.4616 - np.hypot(0.4616, 2.0 * v)
    return np.where(v2 == 0.0, 0.0, exp(log_factor))
