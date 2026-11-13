"""Thermal conductivity of the core: electrons, muons and neutrons. Local T in K, rho in g/cm^3, kappa in CGS.

Composition per baryon: Y_p = Y_e + Y_mu (charge neutrality), Y_n = 1 - Y_p. The nucleon-lepton and
lepton-lepton collision frequencies, the screening and the neutron conductivity are those of the
electron-and-neutron expressions of the model (see README, REFERENCES).

Lepton conduction, ``lepton_conductivity``: kappa_l = (pi^2/3) k_B^2 T n_l c^2 / (mu_l nu_l), written for
ultra-relativistic leptons as C T / nu_l (n_l/n_0)^(2/3) (n_0 = 0.16 fm^-3), times the velocity factor
v_l = p_F / mu_l that carries the finite mass (n_l / mu_l = n_l v_l / p_F). The total collision
frequency nu_l = nu_ll + nu_ln + nu_lp adds lepton-lepton, lepton-neutron and lepton-proton
collisions through a common screening q.

Electrons: v = 1 (m_e << p_F). Muons: the same expression with the muon density and
v = p_F / sqrt(p_F^2 + m_mu^2). Approximations: each lepton is treated as the only lepton
species in its collision frequencies and screening (no electron-muon collisions and no coupled
electron-muon kinetic equations), and the collision frequencies of the muon are those of the
electron expression for its own density. Flux-tube scattering is not included.
"""

import numpy as np

from constants import HBARC_MEV_FM, MUON_MASS_MEV
from effective_mass import effective_mass_ratios
from superfluid import dimensionless_pairing_gap


def _pairing_factors(tau, pairing):
    """Superfluid reduction factors (R, Z, R_C) of the collision frequencies and of the neutron conduction."""
    mask = tau <= 1
    y = dimensionless_pairing_gap(tau[mask], pairing)
    y2 = y * y
    R = (
        0.7694
        + np.hypot(0.2306, 0.07207 * y)
        + y2 * (27.0 + 0.1476 * y2) * np.exp(-np.hypot(4.273, y))
        + 0.5051 * np.expm1(4.273 - np.hypot(4.273, y))
    ) * np.exp(1.187 - np.hypot(1.187, y))
    Z = np.sqrt(0.9443 + np.hypot(0.0557, 0.1886 * y)) * np.exp(1.753 - np.hypot(1.753, y))
    RC = (0.647 + np.sqrt(0.353**2 + 0.109 * y2)) ** 1.5 * np.exp(1.39 - np.hypot(1.39, y))
    outputs = [np.ones_like(tau) for _ in range(3)]
    for out, value in zip(outputs, (R, Z, RC)):
        out[mask] = value
    return outputs


def lepton_conductivity(T8, n_lepton, n_n, n_p, mn, mp, Rn, Rp, Zn, Zp, velocity=1.0):
    """Conductivity of one lepton species and its collision frequencies [s^-1 / 1e15 units of the fit].

    Densities in fm^-3; ``mn`` and ``mp`` the nucleon effective-mass ratios; ``Rn``, ``Rp``, ``Zn``,
    ``Zp`` the superfluid factors; ``velocity`` the factor p_F / mu of the lepton (1 for electrons).
    Returns ``(kappa, nu_ll, nu_ln, nu_lp)``.
    """
    sfn = mn * (n_n * 0.16 / n_lepton**2) ** (1.0 / 3.0) * Zn
    sfp = mp * (n_p * 0.16 / n_lepton**2) ** (1.0 / 3.0) * Zp
    q = np.sqrt(0.00929 * (1.0 + 2.83 * (sfn + sfp)))
    nu_ll = 3.58e11 / q**3 / (n_lepton / 0.16) ** (1.0 / 3.0) * T8**2
    nu_ln = 1.15e12 / q**3 * mn**2 / (n_lepton / 0.16) * Rn * T8**2
    nu_lp = 1.15e12 / q**3 * mp**2 / (n_lepton / 0.16) * Rp * T8**2
    kappa = (
        1.7e24
        * 1.2
        * T8
        * 1.0e15
        / (nu_ll + nu_ln + nu_lp)
        * (n_lepton / 0.16) ** (2.0 / 3.0)
        * velocity
    )
    return kappa, nu_ll, nu_ln, nu_lp


def muon_velocity(n_muon):
    """p_F / mu of the muons of number density ``n_muon`` [fm^-3] (relativistic Fermi gas)."""
    p = HBARC_MEV_FM * (3.0 * np.pi**2 * n_muon) ** (1.0 / 3.0)
    return p / np.sqrt(p**2 + MUON_MASS_MEV**2)


def core_conductivity(T_K, rho, Y_e, tau_n=10.0, tau_p=10.0, *, effective_masses=None, Y_mu=0.0):
    """Electron, muon and neutron conductivity components.

    Inputs are scalars or broadcastable arrays. tau = T/Tc; use tau > 1 when there is no
    superfluidity. ``Y_mu`` is the muon fraction per baryon (0: npe matter). Flux-tube
    scattering is not included.
    """
    T_K, rho, Y_e, tau_n, tau_p, Y_mu = np.broadcast_arrays(
        *[np.asarray(v, float) for v in (T_K, rho, Y_e, tau_n, tau_p, Y_mu)]
    )
    values = np.array([T_K, rho, Y_e, tau_n, tau_p])
    if (
        not np.isfinite(values).all()
        or not np.isfinite(Y_mu).all()
        or np.any(values <= 0)
        or np.any(Y_mu < 0)
        or np.any(Y_e + Y_mu >= 1)
    ):
        raise ValueError(
            "T, rho, tau must be positive and finite; 0 < Y_e, Y_mu >= 0 and Y_e + Y_mu < 1 required"
        )
    T8 = T_K / 1.0e8
    Y_p = Y_e + Y_mu
    Y_n = 1.0 - Y_p
    nb = rho / 1.66e15
    ne = nb * Y_e
    n_mu = nb * Y_mu
    n_p = nb * Y_p
    n_n = nb * Y_n
    mn, mp = (
        effective_mass_ratios(rho / 1.66e-24 * Y_n, rho / 1.66e-24 * Y_p)
        if effective_masses is None
        else effective_masses
    )
    Rn, Zn, RC = _pairing_factors(tau_n, 2)
    Rp, Zp, _ = _pairing_factors(tau_p, 1)
    electrons, nuee, nuen, nuep = lepton_conductivity(T8, ne, n_n, n_p, mn, mp, Rn, Rp, Zn, Zp)
    muons = np.zeros_like(electrons)
    nu_mm, nu_mn, nu_mp = (
        np.zeros_like(electrons),
        np.zeros_like(electrons),
        np.zeros_like(electrons),
    )
    present = n_mu > 0
    if present.any():
        n_safe = np.where(present, n_mu, 1.0)
        kappa, a, b, c = lepton_conductivity(
            T8, n_safe, n_n, n_p, mn, mp, Rn, Rp, Zn, Zp, velocity=muon_velocity(n_safe)
        )
        muons = np.where(present, kappa, 0.0)
        nu_mm, nu_mn, nu_mp = (np.where(present, v, 0.0) for v in (a, b, c))
    # Interaction/superfluid corrections to the neutron collision frequencies are set to 1.
    nunn = 3.48e15 * mn**3 * T8**2
    nunp = 3.48e15 * mn * mp**2 * T8**2
    neutrons = 7.2e23 * T8 / mn * RC**2 * 1.0e15 / (nunn + nunp) * (nb / 0.16 * Y_n)
    return dict(
        electrons=electrons,
        muons=muons,
        neutrons=neutrons,
        total=electrons + muons + neutrons,
        nu_ee=nuee,
        nu_en=nuen,
        nu_ep=nuep,
        nu_mm=nu_mm,
        nu_mn=nu_mn,
        nu_mp=nu_mp,
        nu_nn=nunn,
        nu_np=nunp,
    )
