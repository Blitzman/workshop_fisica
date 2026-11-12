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


HBARC, MN_MEV, ME_MEV, MMU_MEV = 197.3269804, 939.56542, 0.51099895, 105.6583755   # CODATA 2018 [MeV fm, MeV]


def r_dir(kn, kp, kl, ms_n, ms_p, m_l, lower_A=None):
    """Shternin, Baldo & Haensel 2018, Eqs. (8)-(9), closed form. Momenta in fm^-1; valid below the direct-Urca threshold."""
    pn, pp, pl = kn * HBARC, kp * HBARC, kl * HBARC
    mu = np.sqrt(m_l**2 + pl**2)
    lo = (pn - pl) if lower_A is None else lower_A * HBARC
    hi = pn + pl

    def prim_A(k):                      # primitive of 1 / (pp^2 - k^2)^2, for k > pp
        return k / (2 * pp**2 * (pp**2 - k**2)) + np.log((k + pp) / (k - pp)) / (4 * pp**3)
    RA = 2 * (ms_p * MN_MEV)**2 * mu**2 / pl * (prim_A(hi) - prim_A(lo))
    # primitive of k / (pn^2 - k^2)^2 is 1 / (2 (pn^2 - k^2))
    RB = 2 * (ms_n * MN_MEV)**2 * mu**2 / (pl * pp) * (1 / (2 * (pn**2 - (pp + pl)**2)) - 1 / (2 * (pn**2 - (pp - pl)**2)))
    return RA, RB


def alpha_u(kn, kp, kl, ms_n, ms_p, m_l):
    """Their practical expression, Eq. (14)."""
    RA, RB = r_dir(kn, kp, kl, ms_n, ms_p, m_l)
    return 0.6 * RA + 0.2 * RB


def murca(nn, np_, ne, nmu, T9, ms_n, ms_p, inmedium=False):
    """Modified Urca: neutron (140) and proton (142) branches, electrons and muons.

    inmedium=True: the neutron branch uses alpha_U of Shternin et al. 2018 in place of alpha_n beta_n."""
    kn, kp = _kf(nn), _kf(np_)
    base = 8.1e21 * ms_n**3 * ms_p * (np_ / N0)**(1 / 3) * T9**8 * ALPHA_N * BETA_N
    out = {}
    for lep, nl in (('e', ne), ('mu', nmu)):
        kl = _kf(nl)
        v = 1.0 if lep == 'e' else np.where(ne > 0, (nmu / np.where(ne > 0, ne, 1.0))**(1 / 3), 0.0)
        ok = (nl > 0)
        klp = np.where(ok, kl, 1.0)
        qn = np.where(ok, base * v, 0.0)
        qn_std = qn                                           # the proton branch (142) is scaled from the standard rate
        if inmedium:
            below = ok & (kn > kp + klp)                      # only below the direct-Urca threshold
            au = alpha_u(np.where(below, kn, 3.0), np.where(below, kp, 1.0), np.where(below, klp, 1.0), ms_n, ms_p,
                         ME_MEV if lep == 'e' else MMU_MEV)
            qn = np.where(below, qn * au / (ALPHA_N * BETA_N), qn)
        tri = kl + 3 * kp - kn
        qp = np.where(ok & (tri > 0), (qn_std if inmedium else qn) * (ms_p / ms_n)**2 * tri**2 / (8 * klp * kp), 0.0)
        out[f'n_{lep}'], out[f'p_{lep}'] = qn, qp
    return out


def brems(nn, np_, T9, ms_n, ms_p):
    """Nucleon-nucleon bremsstrahlung, Eqs. (165)-(167)."""
    return {'nn': 7.5e19 * ms_n**4 * (nn / N0)**(1 / 3) * A_NN * B_NN * N_NU * T9**8,
            'np': 1.5e20 * (ms_n * ms_p)**2 * (np_ / N0)**(1 / 3) * A_NP * B_NP * N_NU * T9**8,
            'pp': 7.5e19 * ms_p**4 * (np_ / N0)**(1 / 3) * A_PP * B_PP * N_NU * T9**8}
