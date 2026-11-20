"""Lindblom-Detweiler (LD) equations as typed in four sources, plus the corrections found by step 7.

Moved from axion_ns/fullgr.py (lines of the "Lindblom-Detweiler (sin axion)" block, unchanged except that
`NR.x` is now `x` from the notation contract). Sources:
  ZL  Zhao and Lattimer, arXiv:2204.03037 Sec. II (sign of 2 w^2 W in X' differs from the corrected one)
  BU  Burgio, Ferrari, Gualtieri and Schulze, arXiv:1106.2736 appendix A (typos in five of six equations)
  KU  Kunjipurayil, Zhao, Kumar, Agrawal and Prakash, arXiv:2205.02081 (correct X')
  ZH  Zhang, Parisi, Vasquez Flores, Lenzi and Chen, arXiv:2601.21911 appendix A (w^3 in V)
  METFIX  ZL with Q read from the metric and the corrected sign (sw = -1): the base of the extended system (D2)
These are REFERENCES (some knowingly wrong as printed): nothing in a derivation module may import this file.
"""
import sympy as sp

from common.notation import (r, w, l, G, Phi0, lam0, n0, H0, H1, K, W, Vt, Xf, d, eps, p, cs2, ep, psi, eL, eh, nn, pp,
                             GM, b_TOV, Q_TOV, Q_met)

def lindblom_detweiler(Q, b, sw=1):
    """Ecuaciones de LD en la forma de Zhao y Lattimer (arXiv:2204.03037, Sec. II). `sw` = signo del
    término 2 w^2 W en X' (+1 tal como está impreso)."""
    return {
        "H1'": r * d(H1) + (l + 1 + 2 * b * eL + 4 * sp.pi * G * r**2 * eL * (p - eps)) * H1
        - eL * (H0 + K - 16 * sp.pi * G * ep * Vt),
        "K'": r * d(K) - (H0 + (nn + 1) * H1 + (eL * Q - l - 1) * K - 8 * sp.pi * G * ep * eh * W),
        "W'": r * d(W) + (l + 1) * (W + l * eh * Vt) - r**2 * eh * (Xf * sp.exp(-Phi0) / (ep * cs2) + H0 / 2 + K),
        "X'": r * d(Xf) + l * Xf - ep * sp.exp(Phi0) / 2 * (
            (3 * eL * Q - 1) * K - 4 * (nn + 1) * eL * Q / r**2 * Vt + (1 - eL * Q) * H0
            + (r**2 * w**2 * sp.exp(-psi) + nn + 1) * H1
            + (sw * 2 * w**2 * sp.exp(lam0 - psi) - 8 * sp.pi * G * ep * eh
               + r**2 * d(sp.exp(-lam0) * d(psi) / r**2)) * W),
        'H0': H0 * (2 * b + nn + Q) - (8 * sp.pi * G * r**2 * sp.exp(-Phi0) * Xf
                                      - ((nn + 1) * Q - w**2 * r**2 * sp.exp(-psi - 2 * lam0)) * H1
                                      + (nn - w**2 * r**2 * sp.exp(-psi) - Q * (Q * eL - 1)) * K),
        'V': Vt - (Xf / ep - Q / r**2 * sp.exp(Phi0 + lam0) * W - sp.exp(Phi0) * H0 / 2) * sp.exp(Phi0) / w**2,
    }


LD_ZL = lindblom_detweiler(Q_TOV, b_TOV)
LD_ZL_FIX = {"X' con -2w^2": lindblom_detweiler(Q_TOV, b_TOV, sw=-1)["X'"]}

# Burgio, Ferrari, Gualtieri y Schulze (arXiv:1106.2736, apéndice A), tal como está en el LaTeX
_M = GM
_x1 = lambda h0c: d(Xf) - (-l / r * Xf + ep * sp.exp(psi / 2) / 2 * (
    h0c + (r * w**2 * sp.exp(-psi) + l * (l + 1) / (2 * r)) * H1
    + (sp.Rational(3, 2) * d(psi) - 1 / r) * K - l * (l + 1) / r**2 * d(psi) * Vt
    - 2 / r * (4 * sp.pi * G * ep * eh + w**2 * sp.exp(lam0 - psi) - r**2 / 2 * d(sp.exp(-lam0) * d(psi) / r**2)) * W))
_h0 = lambda s: (3 * _M + nn * r + 4 * sp.pi * G * r**3 * p) * H0 - 8 * sp.pi * G * r**3 * sp.exp(-psi / 2) * Xf \
    + (l * (l + 1) / 2 * (_M + 4 * sp.pi * G * r**3 * p) - w**2 * r**3 * sp.exp(-2 * lam0 - psi)) * H1 \
    - (nn * r - w**2 * r**3 * sp.exp(-psi) + s * eL / r * (_M + 4 * sp.pi * G * r**3 * p) * (3 * _M - r + 4 * sp.pi * G * r**3 * p)) * K
_h1 = lambda h1: d(H1) - (-(l + 1 + 2 * _M * eL / r + 4 * sp.pi * G * r**2 * eL * (p - eps)) / r * h1
                          + eL / r * (H0 + K - 16 * sp.pi * G * ep * Vt))
_k = lambda s: d(K) - (H0 / r + l * (l + 1) / (2 * r) * H1 - ((l + 1) / r + s * d(psi) / 2) * K
                       - 8 * sp.pi * G * ep * eh / r * W)
_v = lambda base, k=2: Xf - (w**k * ep * sp.exp(-psi / 2) * Vt - pp / r * base**((psi - 2 * lam0) / 2) * W
                              + sp.exp(psi / 2) / 2 * ep * H0)
LD_BU = {
    "H1'": _h1(1), "K'": _k(1),
    "W'": d(W) - (-(l + 1) / r * W + r * eh * (sp.exp(-psi / 2) / (ep * cs2) * Xf - l * (l + 1) / r**2 * Vt + H0 / 2 + K)),
    "X'": _x1(1 / r + d(psi) / 2), 'H0': _h0(1), 'V': _v(r),
}
# la corrección mínima que hace coincidir cada una
LD_BU_FIX = {
    "H1'": ('falta el factor H_1 en el primer término', _h1(H1)),
    "K'": ('signo de psi\'/2 en el coeficiente de K', _k(-1)),
    "X'": ('falta H_0 y el signo de psi\'/2 es el contrario: (1/r - psi\'/2) H_0', _x1((1 / r - d(psi) / 2) * H0)),
    'H0': ('signo del último término de K', _h0(-1)),
    'V': ('r^{(psi-lambda)/2} donde va e^{(psi-lambda)/2}', _v(sp.E)),
}

# Dos transcripciones más, posteriores, para ver si las discrepancias son de las fuentes o nuestras.
# Kunjipurayil, Zhao, Kumar, Agrawal y Prakash (arXiv:2205.02081, Ecs. 9-14): forma de Zhao y Lattimer,
# con el término de W de X' escrito como -2[w^2 e^{lambda/2-nu} + 4 pi (eps+p) e^{lambda/2} - r^2 (e^{lambda/2} Q/r^3)'].
LD_KU = dict(LD_ZL)
LD_KU["X'"] = r * d(Xf) + l * Xf - ep * sp.exp(Phi0) / 2 * (
    (3 * eL * Q_TOV - 1) * K - 4 * (nn + 1) * eL * Q_TOV / r**2 * Vt + (1 - eL * Q_TOV) * H0
    + (r**2 * w**2 * sp.exp(-psi) + nn + 1) * H1
    - 2 * (w**2 * eh * sp.exp(-psi) + 4 * sp.pi * G * ep * eh - r**2 * d(eh * Q_TOV / r**3)) * W)
# Zhang, Parisi, Vásquez Flores, Lenzi y Chen (arXiv:2601.21911, apéndice A): forma de Burgio et al.
LD_ZH = {"H1'": _h1(H1), "K'": _k(-1), "W'": LD_BU["W'"], "X'": _x1((1 / r - d(psi) / 2) * H0),
         'H0': _h0(-1), 'V': _v(sp.E, 3)}


LD_METFIX = lindblom_detweiler(Q_met, b_TOV, sw=-1)   # fullgr._LDm
