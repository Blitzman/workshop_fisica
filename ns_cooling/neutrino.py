"""Neutrino emissivities in npe-mu matter. Yakovlev, Kaminker, Gnedin & Haensel 2001 (arXiv:astro-ph/0012122).

Densities in fm^-3, T9 = T / 1e9 K, emissivities in erg cm^-3 s^-1. ms_n, ms_p are m*/m.
"""
import numpy as np

N0 = 0.16                                   # fm^-3
ALPHA_N, BETA_N = 1.13, 0.68                # after Eq. (140)
A_NN, A_NP, A_PP = 0.59, 1.06, 0.11         # after Eq. (167)
B_NN, B_NP, B_PP = 0.56, 0.66, 0.7
N_NU = 3


def _kf(n):
    return (3 * np.pi**2 * np.maximum(n, 0.0))**(1 / 3)


def murca(nn, np_, ne, nmu, T9, ms_n, ms_p):
    """Modified Urca: neutron (140) and proton (142) branches, electrons and muons."""
    kn, kp = _kf(nn), _kf(np_)
    base = 8.1e21 * ms_n**3 * ms_p * (np_ / N0)**(1 / 3) * T9**8 * ALPHA_N * BETA_N
    out = {}
    for lep, nl in (('e', ne), ('mu', nmu)):
        kl = _kf(nl)
        v = 1.0 if lep == 'e' else np.where(ne > 0, (nmu / np.where(ne > 0, ne, 1.0))**(1 / 3), 0.0)
        ok = (nl > 0)
        klp = np.where(ok, kl, 1.0)
        qn = np.where(ok, base * v, 0.0)
        tri = kl + 3 * kp - kn
        qp = np.where(ok & (tri > 0), qn * (ms_p / ms_n)**2 * tri**2 / (8 * klp * kp), 0.0)
        out[f'n_{lep}'], out[f'p_{lep}'] = qn, qp
    return out


def brems(nn, np_, T9, ms_n, ms_p):
    """Nucleon-nucleon bremsstrahlung, Eqs. (165)-(167)."""
    return {'nn': 7.5e19 * ms_n**4 * (nn / N0)**(1 / 3) * A_NN * B_NN * N_NU * T9**8,
            'np': 1.5e20 * (ms_n * ms_p)**2 * (np_ / N0)**(1 / 3) * A_NP * B_NP * N_NU * T9**8,
            'pp': 7.5e19 * ms_p**4 * (np_ / N0)**(1 / 3) * A_PP * B_PP * N_NU * T9**8}
