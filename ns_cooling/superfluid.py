"""Nucleon superfluidity: gaps, reduction factors and Cooper-pair neutrino emission.

Sources (equation numbers in the comments):
  Y01  Yakovlev, Kaminker, Gnedin & Haensel 2001, arXiv:astro-ph/0012122 (pp. 91, 95-106, 108-110)
  YLS  Yakovlev, Levenfish & Shibanov 1999, arXiv:astro-ph/9906456, Eq. (18)
  Ho15 Ho, Elshamouty, Heinke & Potekhin 2015, arXiv:1412.7759, Eqs. (1)-(2), Table II
  PC18 Potekhin & Chabrier 2018, arXiv:1711.07662, Sect. 3 (update of the pair-breaking rate)
Type A: singlet (protons in the core); type B: triplet (neutrons in the core). v = gap / k_B T.
"""
import numpy as np
from scipy import constants as K

MEV_K = K.e * 1e6 / K.k                      # 1 MeV in kelvin
HBARC = K.hbar * K.c / (K.e * 1e6) * 1e15    # MeV fm
MN_MEV = K.m_n * K.c**2 / (K.e * 1e6)
N0 = 0.16

# Ho15 Table II: Delta0 [MeV], k0 [fm^-1], k1 [fm^-2], k2 [fm^-1], k3 [fm^-2]
GAPS_P = {'AO': (14, 0.15, 0.22, 1.05, 3.8), 'BCLL': (1.69, 0.05, 0.07, 1.05, 0.16), 'BS': (17, 0.0, 2.9, 0.8, 0.08),
          'CCDK': (102, 0.0, 9.0, 1.3, 1.5), 'CCYms': (35, 0.0, 5.0, 1.1, 0.5), 'CCYps': (34, 0.0, 5.0, 0.95, 0.3),
          'EEHO': (4.5, 0.0, 0.57, 1.2, 0.35), 'EEHOr': (61, 0.0, 6.0, 1.1, 0.6), 'T': (48, 0.15, 2.1, 1.2, 2.8)}
GAPS_NT = {'AO': (4.0, 1.2, 0.45, 3.3, 5.0), 'BEEHS': (0.45, 1.0, 0.40, 3.2, 0.25), 'EEHO': (0.48, 1.28, 0.1, 2.37, 0.02),
           'EEHOr': (0.23, 1.2, 0.026, 1.6, 0.0080), 'SYHHP': (1.0, 2.08, 0.04, 2.7, 0.013), 'T': (1.2, 1.55, 0.05, 2.35, 0.07),
           'TTav': (3.0, 1.1, 0.60, 2.92, 3.0), 'TToa': (2.1, 1.1, 0.60, 3.2, 2.4)}


def gap(kf, par):
    """Ho15 Eq. (2); zero outside [k0, k2]."""
    d0, k0, k1, k2, k3 = par
    kf = np.asarray(kf, float)
    d = d0 * (kf - k0)**2 / ((kf - k0)**2 + k1) * (kf - k2)**2 / ((kf - k2)**2 + k3)
    return np.where((kf > k0) & (kf < k2), d, 0.0)


def tc(kf, par, kind):
    """Critical temperature [K], Ho15 Eq. (1): k T_c = 0.5669 Delta (singlet), 0.1187 Delta (triplet)."""
    return {'A': 0.5669, 'B': 0.1187}[kind] * gap(kf, par) * MEV_K


def v_of_tau(tau, kind):
    """Y01 Eq. (188). Zero for tau >= 1 (and where there is no superfluidity, tau = inf)."""
    tau = np.asarray(tau, float)
    t = np.clip(tau, 1e-6, 1.0)
    if kind == 'A':
        v = np.sqrt(1 - t) * (1.456 - 0.157 / np.sqrt(t) + 1.764 / t)
    else:
        v = np.sqrt(1 - t) * (0.7893 + 1.188 / t)
    return np.where(tau < 1, v, 0.0)


def _s(a, b, v):
    return np.sqrt(a**2 + (b * v)**2)


# --- heat capacity, YLS Eq. (18)
def c_A(v):
    return (0.4186 + _s(1.007, 0.5010, v))**2.5 * np.exp(1.456 - _s(1.456, 1, v))


def c_B(v):
    return (0.6893 + _s(0.790, 0.2824, v))**2 * np.exp(1.934 - _s(1.934, 1, v))


def heat_factor(v, kind):
    f = c_A(v) if kind == 'A' else c_B(v)
    return np.where(v > 0, f, 1.0)             # the factor jumps at T_c; it is 1 in normal matter


# --- direct Urca, Y01 Eq. (199), (204), (205)
def d_A(v):
    return (0.2312 + _s(0.7688, 0.1438, v))**5.5 * np.exp(3.427 - _s(3.427, 1, v))


def d_B(v):
    return (0.2546 + _s(0.7454, 0.1284, v))**5 * np.exp(2.701 - _s(2.701, 1, v))


def d_BA(vn, vp):
    """Neutrons B (v1 = vn), protons A (v2 = vp). Fit (205) for sqrt(vn^2 + vp^2) <= 5, rule (204) beyond."""
    vn, vp = np.broadcast_arrays(np.asarray(vn, float), np.asarray(vp, float))
    fit = (1e4 - 2.839 * vp**4 - 5.022 * vn**4) / (1e4 + 757.0 * vp**2 + 1494 * vn**2 + 211.1 * vp**2 * vn**2 + 0.4832 * vp**4 * vn**4)
    return np.where(np.hypot(vn, vp) <= 5, fit, np.minimum(d_B(vn), d_A(vp)))


# --- modified Urca, Y01 Eqs. (211), (212), (215), (225)-(227)
def mn_pA(v):
    a, b = 0.1477 + _s(0.8523, 0.1175, v), 0.1477 + _s(0.8523, 0.1297, v)
    return (a**7.5 + b**5.5) / 2 * np.exp(3.4370 - _s(3.4370, 1, v))


def mp_pA(v):
    return (0.2414 + _s(0.7586, 0.1318, v))**7 * np.exp(5.339 - _s(5.339, 2, v))


def mp_nB(v):
    a, b = 0.1612 + _s(0.8388, 0.1117, v), 0.1612 + _s(0.8388, 0.1274, v)
    return (a**7 + b**5) / 2 * np.exp(2.398 - _s(2.398, 1, v))


def murca_n(vn, vp):
    """Neutron branch. (227) for neutrons only; (226) when protons are superfluid."""
    both = d_BA(2 * vn, vp) / d_A(vp) * mn_pA(vp)
    return np.where(vp > 0, both, mp_pA(vn))


def murca_p(vn, vp):
    """Proton branch. (212) for protons only; (225) when neutrons are superfluid."""
    both = d_BA(vn, 2 * vp) / d_B(vn) * mp_nB(vn)
    return np.where(vn > 0, both, mp_pA(vp))


# --- bremsstrahlung, Y01 Eqs. (220), (221), (228), (229)
def np_pA(v):
    a, b = 0.9982 + _s(0.0018, 0.3815, v), 0.3949 + _s(0.6051, 0.2666, v)
    return (a * np.exp(1.306 - _s(1.306, 1, v)) + 1.732 * b**7 * np.exp(3.303 - _s(3.303, 2, v))) / 2.732


def pp_pA(v):
    c, d = 0.1747 + _s(0.8253, 0.07933, v), 0.7333 + _s(0.2667, 0.1678, v)
    return (c**2 * np.exp(4.228 - _s(4.228, 2, v)) + d**7.5 * np.exp(7.762 - _s(7.762, 3, v))) / 2


def brems_nn(vn):
    return pp_pA(vn)                                           # (228)


def brems_np(vn, vp):
    return d_BA(vn, vp) / d_A(vp) * np_pA(vp)                  # (229)


def brems_pp(vp):
    return pp_pA(vp)


# --- Cooper-pair emission, Y01 Eqs. (236), (238), (241)
def F_A(v):
    return ((0.602 * v**2 + 0.5942 * v**4 + 0.288 * v**6) * (0.5547 + np.sqrt(0.4453**2 + 0.0113 * v**2))**0.5
            * np.exp(-np.sqrt(4 * v**2 + 2.245**2) + 2.245))


def F_B(v):
    return ((1.204 * v**2 + 3.733 * v**4 + 0.3191 * v**6) / (1 + 0.3511 * v**2)
            * (0.7591 + np.sqrt(0.2409**2 + 0.3145 * v**2))**2 * np.exp(-np.sqrt(4 * v**2 + 0.4616**2) + 0.4616))


def pbf(n, T9, v, ms, species, pbf='Y01', n_nu=3):
    """Emissivity [erg cm^-3 s^-1]. species 'n' (triplet, a = 4.17) or 'p' (singlet, a from (238))."""
    pf_mc = 0.353 * (np.maximum(n, 0.0) / N0)**(1 / 3)         # p_F / (m_N c), text after (237)
    vf = pf_mc / ms                                            # v_F / c
    if species == 'n':
        a, F = 4.17 * (0.19 if pbf == 'updated' else 1.0), F_B(v)       # PC18: triplet rate times 0.19
    else:
        vec = 0.0064 * (vf**2 if pbf == 'updated' else 1.0)             # PC18: vector part suppressed ~ (p_F/m*c)^2
        a, F = vec + 1.59 * vf**2 * (ms**2 + 11 / 42), F_A(v)
    return 1.170e21 * ms * pf_mc * T9**7 * n_nu * a * F
