"""Ecuación de estado BSk26 (Goriely et al. 2013; Pearson et al. 2018) en forma analítica.

Ajustes de Potekhin et al. (bskfit18.f, http://www.ioffe.ru/astro/NSG/BSk/), los mismos que usa
la literatura para BSk22-26:
  - BSKfit: lg P(lg rho), con rho la densidad de masa-energía (g/cm^3);
  - BSkEofN: energía por barión en función de n_b, de donde rho(n_b).
Con las dos se construye una tabla (n_b, eps, p) y se interpola en logaritmos.

Unidades de salida: MeV/fm^3 para eps y p, fm^-3 para n_b.
"""
import numpy as np
from scipy.interpolate import CubicSpline

MEV_FM3_TO_GCC = 1.782662e12        # MeV/fm^3 -> g/cm^3 (constante de bskfit18.f)
MEV_FM3_TO_DYN = 1.602176634e33     # MeV/fm^3 -> dyn/cm^2
M_N_MEV = 939.5654                  # masa del neutrón (bskfit18.f)
GROUND = -9.1536                    # m(Fe56)/56 - m_n (bskfit18.f)

# BSKfit, columna BSk26 de bskfit18.f
_A = [3.672, 7.844, .00876, .22604, 3.129, 11.939, 13.738, 1.3389,
      3.112, -23.031, 1.6264, 4.83, 14.272, 23.28, -1.542, 2.10,
      15.31, .083, 11.66, 6.16, -.042, 14.18, 14.8]
# BSkEofN, columna BSk26
_B = [1.3995e8, 1.066e9, 7.472e8, 3.599e7, 2.906e6,
      .3051, -.54068, 388.5, 776.5, 15.435,
      18.66, .3029, 2.2483, .569]


def _fermi(x):
    return 1.0 / (np.exp(np.clip(x, -40, 40)) + 1.0)


def lgP_of_lgrho(x):
    """lg P [dyn/cm^2] en función de lg rho [g/cm^3] (BSKfit)."""
    a = [None] + _A
    p01 = (a[1] + a[2] * x + a[3] * x**3) / (1 + a[4] * x)
    p = (p01 * _fermi(a[5] * (x - a[6])) + (a[7] + a[8] * x) * _fermi(a[9] * (a[6] - x))
         + (a[10] + a[11] * x) * _fermi(a[12] * (a[13] - x)) + (a[14] + a[15] * x) * _fermi(a[16] * (a[17] - x))
         + a[18] / (1 + ((x - a[19]) * a[20])**2) + a[21] / (1 + ((x - a[22]) * a[23])**2))
    return p


def eps_of_n(n):
    """Densidad de energía total [MeV/fm^3] en función de n_b [fm^-3] (BSkEofN)."""
    b = [None] + _B
    EH = (b[10] * n)**b[13] / (1 + b[12] * n)
    WH = 1 / (1 + (b[11] * n)**b[14])
    EM = b[6] * n**b[7] * (1 + b[8] * n)
    WM = 1 / (1 + b[9] * n)
    EL = (b[1] * n)**1.16667 / (1 + np.sqrt(b[2] * n))
    EL = EL / (1 + np.sqrt(b[3] * n)) * (1 + np.sqrt(b[4] * n)) / (1 + np.sqrt(b[5] * n))
    E = EL * WM + EM * (1 - WM) * WH + EH * (1 - WH)
    return n * (M_N_MEV + E + GROUND)


class EOS:
    """Tabla monótona en lg p: da eps(p), n_b(p) y c_s^2 = dp/deps."""

    def __init__(self, n_min=1e-8, n_max=1.6, N=4000):
        n = np.logspace(np.log10(n_min), np.log10(n_max), N)
        e = eps_of_n(n)
        p = 10**lgP_of_lgrho(np.log10(e * MEV_FM3_TO_GCC)) / MEV_FM3_TO_DYN
        ok = np.concatenate([[True], np.diff(p) > 0])     # la tabla debe ser monótona en p
        self.n, self.e, self.p = n[ok], e[ok], p[ok]
        lp = np.log(self.p)
        self._le = CubicSpline(lp, np.log(self.e))
        self._ln = CubicSpline(lp, np.log(self.n))
        # c_s^2 = dp/deps a lo largo del ajuste P(rho), coherente con eps(p)
        self.p_min = self.p[0]

    def eps(self, p):
        return np.exp(self._le(np.log(p)))

    def nb(self, p):
        return np.exp(self._ln(np.log(p)))

    def cs2(self, p):
        return p / self.eps(p) / self._le(np.log(p), 1)
