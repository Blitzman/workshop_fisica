"""Volumetric heat capacity [erg cm^-3 K^-1].

Each contribution is a function ``term(state)`` registered in
``HEAT_CAPACITY_TERMS``; the total is the sum over the registry. ``state`` is a
``state.LocalState`` (temperature and background at the points of interest). A
term returns an array of the same shape as ``state.T``.

    @HEAT_CAPACITY_TERMS.register("my_term")
    def my_term(state):
        return 1e18 * state.T * (state.xh == 0)

Contributions
    electrons, muons   degenerate relativistic lepton gases (full dependence on k_F / m c)
    neutrons, protons  Landau quasiparticles with the effective masses of the background,
                       reduced by superfluidity (singlet crust neutrons, triplet core
                       neutrons, singlet core protons)
    ions               gas, liquid and Debye-solid regimes of the crust nuclei
"""

import numpy as np

from constants import HBAR, HBARC_MEV_FM, K_B, M_NEUTRON, M_PROTON, MEV_TO_ERG, MUON_MASS_MEV
from numerical_validation import checked
from registry import Registry
from superfluid import dimensionless_pairing_gap

HEAT_CAPACITY_TERMS = Registry("heat-capacity term")


def heat_capacity_pairing_factor(T_K, Tc_K, pairing_type):
    """Superfluid/normal heat-capacity ratio.

    At T_K = Tc_K the superfluid branch is used. For Tc_K = 0 or T_K > Tc_K the
    factor is 1 (no division by Tc_K). The jump at Tc_K is physical in this
    prescription: the factors are deliberately not normalized to R(0) = 1.
    """
    T_K, Tc_K = np.broadcast_arrays(checked(T_K, "T", True), checked(Tc_K, "Tc"))
    R = np.ones_like(T_K)
    mask = (Tc_K > 0) & (T_K <= Tc_K)
    v = dimensionless_pairing_gap(T_K[mask] / Tc_K[mask], pairing_type)
    if pairing_type == 1:
        a = 0.4186 + np.hypot(1.007, 0.5010 * v)
        z = 1.456 - np.hypot(1.456, v)
        R[mask] = a**2.5 * np.exp(np.maximum(z, -150.0))
    elif pairing_type == 2:
        a = 0.6893 + np.hypot(0.790, 0.2824 * v)
        z = 1.934 - np.hypot(1.934, v)
        R[mask] = a**2 * np.exp(np.maximum(z, -150.0))
    else:
        raise ValueError("pairing_type must be 1 or 2")
    return R


def debye_heat_capacity_factor(x):
    """Piecewise fit of the Debye function (x = T/T_D); its small jumps are not smoothed."""
    x = checked(x, "T/T_D")
    f = np.empty_like(x)
    lo, hi = x <= 0.15, x >= 0.4
    mid = ~(lo | hi)
    f[lo] = 77.9273 * x[lo] ** 3
    f[hi] = 1.0 - (0.05 / x[hi]) / x[hi]
    f[mid] = 1.69798 * x[mid] + 0.0083073
    return f


def landau_heat_capacity(kf, mass_ratio, particle_mass, T):
    """Degenerate nucleons as Landau quasiparticles, C = k_B^2 T m* k_F / (3 hbar^2) [erg cm^-3 K^-1].

    ``kf`` is the Fermi momentum [fm^-1], ``mass_ratio`` the Landau effective mass m*/m and
    ``particle_mass`` the bare mass [g]. The density of states at the Fermi surface is
    m* k_F / (pi^2 hbar^2) when m* is the Landau mass, k_F / v_F: the relativistic and the
    interaction corrections are contained in m*, so no further relativistic factor appears.
    """
    return K_B**2 / (3.0 * HBAR**2) * mass_ratio * particle_mass * (kf * 1e13) * T


def degenerate_lepton_heat_capacity(kf, mass_mev, T):
    """Degenerate relativistic free lepton gas, C = k_B^2 T E_F k_F / (3 (hbar c)^2) [erg cm^-3 K^-1].

    ``kf`` is the Fermi momentum [fm^-1], ``mass_mev`` the rest energy m c^2 [MeV] and
    E_F = sqrt((hbar c k_F)^2 + m^2 c^4) the Fermi energy including the rest mass.
    The ultra-relativistic limit is k_B^2 T k_F^2 / (3 hbar c), the non-relativistic one
    k_B^2 T m k_F / (3 hbar^2).
    """
    energy = np.sqrt((HBARC_MEV_FM * kf) ** 2 + mass_mev**2) * MEV_TO_ERG
    hbar_c = HBARC_MEV_FM * 1e-13 * MEV_TO_ERG  # erg cm
    return K_B**2 * T * energy * (kf * 1e13) / (3.0 * hbar_c**2)


@HEAT_CAPACITY_TERMS.register("electrons")
def electron_heat_capacity(state):
    """Degenerate electrons."""
    cv = np.zeros_like(state.T)
    mask = state.ye > 0
    xe = 0.010067 * (state.ye[mask] * state.rho[mask]) ** (1.0 / 3.0)
    cv[mask] = (
        state.nb[mask] * state.ye[mask] * 2.298976e-25 * state.T[mask] * np.hypot(xe, 1.0) / xe**2
    )
    return cv


@HEAT_CAPACITY_TERMS.register("muons")
def muon_heat_capacity(state):
    """Degenerate muons of the core (zero where there are none), with their full relativistic dependence."""
    return degenerate_lepton_heat_capacity(state.kf_mu, MUON_MASS_MEV, state.T)


@HEAT_CAPACITY_TERMS.register("neutrons")
def neutron_heat_capacity(state):
    """Free neutrons, reduced by crust (type 1) or core (type 2) pairing when superfluid."""
    cv = landau_heat_capacity(state.kf_n, state.mstar_n, M_NEUTRON, state.T)
    if state.superfluid:
        # Explicit regional selection, also in cells of mixed composition.
        cv *= heat_capacity_pairing_factor(state.T, np.where(state.crust, state.tc_crust, 0.0), 1)
        cv *= heat_capacity_pairing_factor(state.T, np.where(state.core, state.tc_n, 0.0), 2)
    return cv


@HEAT_CAPACITY_TERMS.register("protons")
def proton_heat_capacity(state):
    """Free protons, reduced by pairing when superfluid."""
    cv = landau_heat_capacity(state.kf_p, state.mstar_p, M_PROTON, state.T)
    if state.superfluid:
        cv *= heat_capacity_pairing_factor(state.T, state.tc_p, 1)
    return cv


@HEAT_CAPACITY_TERMS.register("ions")
def ion_heat_capacity(state):
    """Ions: gas, liquid and Debye solid regimes depending on the Coulomb parameter."""
    cv = np.zeros_like(state.T)
    ions = state.crust
    # Nothing (0/A, T/T_D, log(gamma)) is evaluated where there are no ions.
    ti, ri, ai, zi, xi = (
        state.T[ions],
        state.rho[ions],
        state.aa[ions],
        state.zz[ions],
        state.xh[ions],
    )
    ni = ri / 1.66e-24 * xi / ai
    T_D_K = 3.48e3 * np.sqrt(ri) * zi / ai
    gamma = 2.275e5 * zi**2 * (ri * xi / ai) ** (1.0 / 3.0) / ti
    ci = np.empty_like(ni)
    weak, solid = gamma <= 1.0, gamma >= 150.0
    liquid = ~(weak | solid)
    ci[weak] = ni[weak] * 4.143e-16 / 2.0
    ci[solid] = ni[solid] * 4.143e-16 * debye_heat_capacity_factor(ti[solid] / T_D_K[solid])
    ci[liquid] = ni[liquid] * 4.143e-16 / 2.0 * (1.0 + np.log10(gamma[liquid]) / np.log10(150.0))
    cv[ions] = ci
    return cv


def heat_capacity(state, *, terms=None):
    """Heat capacity of every registered term, plus 'total' [erg cm^-3 K^-1].

    ``terms`` optionally supplies an ordered snapshot of (name, callable) pairs.
    The default terms are electrons, muons, free neutrons, free protons and ions.
    T = 0 is not allowed: that limit must be treated separately, outside the thermal solver.
    """
    selected = HEAT_CAPACITY_TERMS.items() if terms is None else terms
    result = {name: term(state) for name, term in selected}
    result["total"] = sum(result.values())
    if any(not np.isfinite(v).all() or np.any(v < 0) for v in result.values()):
        raise FloatingPointError("Heat capacity not finite or negative")
    return result


def heat_capacity_internal(state):
    """Total heat capacity in internal units (1e40 erg km^-3 per 1e8 K), one value per point."""
    return heat_capacity(state)["total"] * 1e23 * 1e-40
