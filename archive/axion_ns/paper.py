"""Ecuaciones de arXiv:2510.16245 transcritas tal cual (forma 'lhs - rhs').

Variables del paper -> variables del modelo:
  p, eps         -> p_of(n0), eps_of(n0, a0)
  M(r)           -> r (1 - e^{-2 lambda0}) / (2G)      [lambda = -1/2 ln(1 - 2GM/r)]
  delta p        -> (dp/dn) delta n                    [Ec. (3.2b) del paper]
  delta eps      -> (eps+p)/n delta n + n m*' delta a  [primera ley con eps(n, a)]
  c_s^2          -> (dp/dn) n / (eps + p)
"""
import sympy as sp

from .model import r, G, w, V, mstar, eps_of, p_of, Phi0, lam0, a0, n0, dPhi, dlam, da1, dn1, xi1

d = lambda f, k=1: sp.diff(f, r, k)
pi = sp.pi

eps, p = eps_of(n0, a0), p_of(n0)
M = r * (1 - sp.exp(-2 * lam0)) / (2 * G)
Va = sp.diff(V(a0), a0)
Vaa = sp.diff(V(a0), a0, 2)
ma = sp.diff(mstar(a0), a0)
maa = sp.diff(mstar(a0), a0, 2)
nb = n0
x_ = sp.Dummy('x')
dpdn = sp.diff(p_of(x_), x_).subs(x_, n0)
cs2 = dpdn * n0 / (eps + p)
dp = dpdn * dn1
de = (eps + p) / n0 * dn1 + n0 * ma * da1

# ---- fondo, Ec. (2.6a)-(2.6d)
BACKGROUND = {
    '(2.6a) p\'': d(p) + (eps + p) * (d(Phi0) + nb / (eps + p) * ma * d(a0)),
    '(2.6b) Phi\'': d(Phi0) - 1 / r**2 / (1 - 2 * G * M / r) * (
        G * M + 4 * pi * G * r**3 * (p - V(a0) + d(a0)**2 / 2 * (1 - 2 * G * M / r))),
    '(2.6c) M\'': d(M) - 4 * pi * r**2 * (eps + V(a0) + d(a0)**2 / 2 * (1 - 2 * G * M / r)),
    '(2.6d) KG': d(a0, 2) * (1 - 2 * G * M / r) + 2 * d(a0) / r * (
        1 - G * M / r - 2 * pi * G * r**2 * (eps - p + 2 * V(a0))) - (Va + nb * ma),
}

# ---- perturbaciones radiales, Ec. (3.3a)-(3.3e)
e2l, e2P = sp.exp(-2 * lam0), sp.exp(-2 * Phi0)
PERTURBATIONS = {
    '(3.3a) KG': e2l * d(da1, 2) + e2l * d(da1) * (d(Phi0) - d(lam0) + 2 / r)
    - 2 * dlam * (Va + nb * ma) + e2l * d(a0) * (d(dPhi) - d(dlam))
    - (Vaa + nb * maa - w**2 * e2P) * da1 - dn1 * ma,
    '(3.3b) rr': d(dPhi) - dlam * (1 / r + 2 * d(Phi0))
    - 4 * pi * G * r * sp.exp(2 * lam0) * (dp - Va * da1 - dlam * e2l * d(a0)**2 + e2l * d(a0) * d(da1)),
    '(3.3c) tr': dlam + 4 * pi * G * r * sp.exp(2 * lam0) * ((eps + p) * xi1 - e2l * d(a0) * da1),
    '(3.3d) Euler': d(dp) + (de + dp) * d(Phi0) + (eps + p) * d(dPhi)
    - w**2 * (eps + p) * sp.exp(2 * (lam0 - Phi0)) * xi1
    + ma * (dn1 * d(a0) + nb * d(da1)) + nb * maa * d(a0) * da1,
    '(3.3e) continuidad': d((eps + p) * xi1) + (eps + p) * xi1 * (2 / r + d(lam0) + d(Phi0))
    + (dp / cs2 + (eps + p) * dlam),
}
