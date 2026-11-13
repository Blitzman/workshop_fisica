"""Neutrino emissivities of the core from Bottaro, Caputo and Fiorillo (arXiv:2406.18640, JCAP 2024).

Nucleon bremsstrahlung and modified Urca re-derived in the framework of Friman and Maxwell
(one-pion plus rho-meson exchange, Landau parameters for the short-range part). The closed
forms of the paper are implemented, in natural units (hbar = c = k_B = 1, energies in MeV):

    nn and pp bremsstrahlung   Eq. 3.11 with J of Eqs. 3.7-3.9      (identical nucleons, p_F of the species)
    np bremsstrahlung          Eq. 3.15 with the leading J of Eq. 3.16
    modified Urca, neutron     Eq. 3.28 (triangular approximation), electrons and muons
    modified Urca, proton      Eq. 3.39 (approximation of the maximum exchanged momenta), electrons
                               and muons, each lepton with its own threshold p_l + 3 p_p > p_n

Differences with the standard (Friman-Maxwell) rates implemented in ``emissivity.py``: the
rho-meson suppression of the tensor and spin-spin interactions, the exchange (u-channel)
and interference terms, the Landau-parameter structure of the amplitude, the proton Fermi
momentum (not the electron one) in the modified-Urca phase space, and the lepton factor
v_l = p_l / mu_l that makes the muon channel small where the muons are not relativistic.

Not included: the non-triangular modified-Urca integral (Eqs. 3.34 and 3.37, a four-dimensional
numerical integral per point, which the paper finds larger than the triangular result by
20-30% for the neutron branch), the in-medium effects the paper leaves aside, and the crust
(the rates apply where ``state.xh == 0``; the crust keeps the standard neutron bremsstrahlung).
The paper's emissivities are for one neutrino flavor; the bremsstrahlung is multiplied by 3
flavors, the modified Urca emits the flavor of its lepton only.

The short-range suppression factors of Friman and Maxwell (beta = 0.56 nn, 0.66 np, 0.70 pp,
0.68 modified Urca) are kept, as in the paper's comparisons (``short_range_suppression``).
Nucleon effective masses enter as overall powers (m*_n/m_n)^a (m*_p/m_p)^b; the Landau parameters
use the bare nucleon mass.

Use ``updated_rates()`` (context manager) or ``use_updated_rates()`` to replace the standard
processes in the registry.
"""

from contextlib import contextmanager

import numpy as np

from constants import (
    ELECTRON_MASS_MEV,
    FERMI_COUPLING_MEV,
    HBAR,
    HBARC_MEV_FM,
    K_B,
    MEV_TO_ERG,
    MUON_MASS_MEV,
    NEUTRON_MASS_MEV,
    PROTON_MASS_MEV,
)
import emissivity as standard
from emissivity import NEUTRINO_PROCESSES, process_factors

# Couplings and meson parameters of the paper (Sec. 2).
G_A = 1.27  # axial coupling
COS_CABIBBO = 0.974
M_PI = 138.0  # pion mass [MeV]
M_RHO = 769.0  # rho-meson mass [MeV]
C_RHO = 1.4  # rho-meson coupling strength
F_PI = 1.0  # pi-nucleon coupling (f_pi^2 / 4 pi = 0.08 is f_pi ~ 1)
F0_PRIME, G0 = 0.7, 1.1  # Landau parameters F0' and G0 (F0 and G0' drop out of the rates)
N_FLAVORS = 3  # neutrino flavors of the bremsstrahlung
BETA = dict(nn=0.56, np=0.66, pp=0.70, murca=0.68)  # short-range suppression of Friman and Maxwell

_G = FERMI_COUPLING_MEV  # G_F
_G_CHARGED = FERMI_COUPLING_MEV * COS_CABIBBO  # G = G_F cos(theta_C)
# Q [MeV^5] -> erg cm^-3 s^-1: (MeV/hbar c)^3 * (MeV / hbar) * erg per MeV
MEV5_TO_CGS = MEV_TO_ERG * 1e39 / HBARC_MEV_FM**3 * MEV_TO_ERG / HBAR


# ------------------------------------------------------------------ the functions of the paper
def phi(a):
    """phi(alpha), Eq. 3.8: the integral of x^4 / (x^2 + alpha^2)^2 over [0, 1]."""
    return 1.0 + a**2 / (2.0 * (1.0 + a**2)) - 1.5 * a * np.arctan(1.0 / a)


def psi(a):
    """psi(alpha), Eq. 3.17: the integral of x^2 / (x^2 + alpha^2) over [0, 1]."""
    return 1.0 - a * np.arctan(1.0 / a)


def cap_phi(a, b):
    """Phi(alpha, beta), Eq. 3.9: the integral of x^4 / ((x^2 + alpha^2)(x^2 + beta^2)) over [0, 1]."""
    return 1.0 - (a**3 * np.arctan(1.0 / a) - b**3 * np.arctan(1.0 / b)) / (a**2 - b**2)


def cap_psi(a, b):
    """Psi(alpha, beta), below Eq. 3.9 (symmetric in its arguments)."""
    r = np.sqrt(1.0 + a * a + b * b)
    return (
        1.0 - a * np.arctan(1.0 / a) - b * np.arctan(1.0 / b) + a * b / r * np.arctan(r / (a * b))
    )


def eta(mass, p):
    """eta(m) = p^2 / (m^2 + p^2) (Eq. 3.16 and below Eq. 3.28)."""
    return p**2 / (mass**2 + p**2)


def landau_f_minus_g(p_n):
    """f' - g [MeV^-2] from the Landau parameters: pi^2 (F0' - G0) / (2 m_N p_F(n))."""
    return np.pi**2 * (F0_PRIME - G0) / (2.0 * NEUTRON_MASS_MEV * p_n)


# ------------------------------------------------------------------ emissivities in natural units
def identical_nucleon_J(p_f):
    """J of Eqs. 3.6-3.9 for nn (p_F neutron) or pp (p_F proton) bremsstrahlung [MeV^-3]."""
    a, b = M_PI / (2.0 * p_f), M_RHO / (2.0 * p_f)
    bracket = (
        2.0 * (phi(a) + C_RHO**2 * phi(b) - 2.0 * C_RHO * cap_phi(a, b))
        + cap_psi(a, a)
        - 2.0 * C_RHO * cap_psi(a, b)
        + C_RHO**2 * cap_psi(b, b)
    )
    return (F_PI / M_PI) ** 4 * np.pi * p_f * bracket


def identical_nucleon_bremsstrahlung(p_f, mass, T):
    """nn or pp bremsstrahlung, one neutrino flavor [MeV^5] (Eq. 3.11); ``mass`` and T in MeV."""
    return 41.0 * G_A**2 * _G**2 * mass**4 * T**8 * identical_nucleon_J(p_f) / (56700.0 * np.pi**2)


def neutron_proton_J(p_n, p_p):
    """Leading J^np of Eq. 3.16 [MeV^-3] (p_p << p_n)."""
    a, b = M_PI / (2.0 * p_p), M_RHO / (2.0 * p_p)
    ep, er = eta(M_PI, p_n), eta(M_RHO, p_n)
    ratio = (F_PI / M_PI) ** 2
    g_minus_f = -landau_f_minus_g(p_n)
    meson = (F_PI / M_PI) ** 4 * (
        2.0 * ep**2
        + 4.0 * C_RHO**2 * er**2
        + phi(a)
        + 3.0 * C_RHO**2 * phi(b)
        + 2.0 * C_RHO * cap_phi(a, b)
        - 2.0 * (ep + C_RHO * er) * (psi(a) + C_RHO * psi(b))
        - 4.0 * C_RHO**2 * er * psi(b)
    )
    landau = 6.0 * g_minus_f**2 + 4.0 * ratio * g_minus_f * (
        psi(a) + 2.0 * C_RHO * psi(b) - ep - 2.0 * C_RHO * er
    )
    return np.pi * p_p * (meson + landau)


def neutron_proton_bremsstrahlung(p_n, p_p, m_n, m_p, T):
    """np bremsstrahlung, one neutrino flavor [MeV^5] (Eq. 3.15)."""
    return (
        41.0
        * G_A**2
        * _G**2
        * m_n**2
        * m_p**2
        * T**8
        * neutron_proton_J(p_n, p_p)
        / (14175.0 * np.pi**2)
    )


def lepton_velocity(p, mass):
    """v = p_F / mu for a lepton of rest energy ``mass`` (all in MeV)."""
    return p / np.sqrt(p**2 + mass**2)


def murca_neutron_branch(p_n, p_p, lepton_velocities, m_n, m_p, T):
    """Modified Urca, neutron branch, triangular approximation [MeV^5] (Eq. 3.28).

    ``lepton_velocities`` is the sum of p_l / mu_l over the leptons (electrons, muons). The factor 2
    of the inverse reaction is included.
    """
    fpg = landau_f_minus_g(p_n)
    pion_rho = eta(M_PI, p_n) - C_RHO * eta(M_RHO, p_n)
    reduced = 0.5 * 21.0 / 16.0 * (32.0 / 7.0 * (M_PI / F_PI) ** 4 * fpg**2 + 2.0 * pion_rho**2)
    return (
        2.0
        * lepton_velocities
        * 11513.0
        / (120960.0 * np.pi)
        * _G_CHARGED**2
        * G_A**2
        * (F_PI / M_PI) ** 4
        * m_n**3
        * m_p
        * p_p
        * T**8
        * reduced
    )


def murca_proton_branch(p_n, p_p, p_lepton, mu_lepton, m_n, m_p, T):
    """Modified Urca, proton branch, one lepton species [MeV^5] (Eq. 3.39).

    Open where p_l + 3 p_p > p_n, a threshold of each lepton. The factor 2 of the inverse reaction
    is included.
    """
    q = p_n - p_p  # the maximum momentum exchanged
    fpg = landau_f_minus_g(p_n)
    bracket = 2.0 * fpg**2 + (F_PI / M_PI) ** 4 * (eta(M_PI, q) - C_RHO * eta(M_RHO, q)) ** 2
    threshold = np.maximum(p_lepton + 3.0 * p_p - p_n, 0.0)
    return (
        2.0
        * threshold**2
        / mu_lepton
        * 11513.0
        / (645120.0 * np.pi)
        * _G_CHARGED**2
        * G_A**2
        * m_n
        * m_p**3
        * T**8
        * bracket
    )


# ------------------------------------------------------------------ processes of the registry
def _setup(state):
    """Core mask, momenta [MeV], masses [MeV] and temperature [MeV] of the points where the rates apply."""
    core = (state.xh == 0) & (state.yn > 0) & (state.yp > 0) & (state.kf_n > 0) & (state.kf_p > 0)
    pick = lambda a: a[core]  # noqa: E731
    return dict(
        core=core,
        p_n=HBARC_MEV_FM * pick(state.kf_n),
        p_p=HBARC_MEV_FM * pick(state.kf_p),
        p_e=HBARC_MEV_FM * pick(state.kf_e),
        p_mu=HBARC_MEV_FM * pick(state.kf_mu),
        mn=NEUTRON_MASS_MEV * pick(state.mstar_n),
        mp=PROTON_MASS_MEV * pick(state.mstar_p),
        T=K_B * pick(state.T) / MEV_TO_ERG,
    )


def make_processes(short_range_suppression=True):
    """The four updated processes as ``{registry name: function(state)}``."""
    beta = BETA if short_range_suppression else dict(nn=1.0, np=1.0, pp=1.0, murca=1.0)

    def modified_urca(state):
        out = np.zeros_like(state.T)
        k = _setup(state)
        f = process_factors(state)
        core = k["core"]
        ve = lepton_velocity(k["p_e"], ELECTRON_MASS_MEV)
        vmu = np.where(k["p_mu"] > 0, lepton_velocity(k["p_mu"], MUON_MASS_MEV), 0.0)
        neutron = murca_neutron_branch(k["p_n"], k["p_p"], ve + vmu, k["mn"], k["mp"], k["T"])
        proton = murca_proton_branch(
            k["p_n"],
            k["p_p"],
            k["p_e"],
            np.sqrt(k["p_e"] ** 2 + ELECTRON_MASS_MEV**2),
            k["mn"],
            k["mp"],
            k["T"],
        )
        has_mu = k["p_mu"] > 0
        proton = proton + np.where(
            has_mu,
            murca_proton_branch(
                k["p_n"],
                k["p_p"],
                k["p_mu"],
                np.sqrt(k["p_mu"] ** 2 + MUON_MASS_MEV**2),
                k["mn"],
                k["mp"],
                k["T"],
            ),
            0.0,
        )
        out[core] = beta["murca"] * (neutron * f.Rmn[core] + proton * f.Rmp[core]) * MEV5_TO_CGS
        return out

    def brems_nn(state):
        out = standard.neutron_neutron_bremsstrahlung(state)  # the crust keeps the standard rate
        k = _setup(state)
        f = process_factors(state)
        q = N_FLAVORS * beta["nn"] * identical_nucleon_bremsstrahlung(k["p_n"], k["mn"], k["T"])
        out[k["core"]] = q * f.Rnn[k["core"]] * MEV5_TO_CGS
        return out

    def brems_np(state):
        out = np.zeros_like(state.T)
        k = _setup(state)
        f = process_factors(state)
        q = (
            N_FLAVORS
            * beta["np"]
            * neutron_proton_bremsstrahlung(k["p_n"], k["p_p"], k["mn"], k["mp"], k["T"])
        )
        out[k["core"]] = q * f.Rnp[k["core"]] * MEV5_TO_CGS
        return out

    def brems_pp(state):
        out = np.zeros_like(state.T)
        k = _setup(state)
        f = process_factors(state)
        q = N_FLAVORS * beta["pp"] * identical_nucleon_bremsstrahlung(k["p_p"], k["mp"], k["T"])
        out[k["core"]] = q * f.Rpp[k["core"]] * MEV5_TO_CGS
        return out

    return dict(
        modified_urca=modified_urca, brems_nn=brems_nn, brems_np=brems_np, brems_pp=brems_pp
    )


_STANDARD = {
    "modified_urca": standard.modified_urca_emissivity,
    "brems_nn": standard.neutron_neutron_bremsstrahlung,
    "brems_np": standard.neutron_proton_bremsstrahlung,
    "brems_pp": standard.proton_proton_bremsstrahlung,
}


def use_updated_rates(short_range_suppression=True):
    """Replace modified Urca and the nn, np and pp bremsstrahlung of the registry by the updated rates."""
    for name, function in make_processes(short_range_suppression).items():
        NEUTRINO_PROCESSES.register(name, replace=True)(function)


def use_standard_rates():
    """Restore the standard (Friman-Maxwell) modified Urca and nucleon bremsstrahlung."""
    for name, function in _STANDARD.items():
        NEUTRINO_PROCESSES.register(name, replace=True)(function)


@contextmanager
def updated_rates(short_range_suppression=True):
    """Use updated rates temporarily; restore the exact previous registry even on error."""
    previous = NEUTRINO_PROCESSES.items()
    use_updated_rates(short_range_suppression)
    try:
        yield
    finally:
        NEUTRINO_PROCESSES.restore(previous)
