"""Integración numérica de las ecuaciones derivadas, sin axión, para un politropo.

Sirve para validar con números lo que el paso 7 valida con álgebra: el sistema extendido de la
página (fullgr.EXT) con el axión apagado, junto con las condiciones en el centro y el empalme con
Zerilli, tiene que dar las frecuencias publicadas. Lo mismo para Cowling con las (8a)-(8c) del paso 5.
Referencia: Yoshida y Kojima (1997, MNRAS 289, 117; gr-qc/9705081), Tablas 1-2: politropo n = 1,
M/R = 0.1 y 0.2, modos f, p1 y p2 con l = 2..5, cálculo exacto y Cowling, en unidades de sqrt(M/R^3).

Unidades G = c = 1 y K = 1 (los politropos son invariantes de escala).
"""
import numpy as np
import sympy as sp
from scipy.integrate import solve_ivp
from scipy.optimize import brentq, minimize_scalar

from . import fullgr as F
from . import nonradial as NR
from .model import r, Phi0, lam0, a0, n0, w, G, V, mstar, mN, eps0

H0, H1, K, W, Vt, dn, dth, Xf = F.H0, F.H1, F.K, F.W, F.Vt, F.dn, F.dth, F.Xf
d, l = F.d, NR.l

# Yoshida y Kojima (1997), Tablas 1 y 2 (n = 1): {M/R: {l: {modo: (exacto, Cowling)}}}
YK97 = {0.1: {2: {'f': (1.21, 1.52), 'p1': (3.11, 3.43), 'p2': (4.83, 5.09)},
              3: {'f': (1.62, 1.82), 'p1': (3.61, 3.85), 'p2': (5.36, 5.58)},
              4: {'f': (1.93, 2.07), 'p1': (4.03, 4.22), 'p2': (5.84, 6.01)},
              5: {'f': (2.18, 2.29), 'p1': (4.40, 4.55), 'p2': (6.27, 6.42)}},
        0.2: {2: {'f': (1.17, 1.35), 'p1': (2.70, 2.99), 'p2': (4.12, 4.38)},
              3: {'f': (1.52, 1.64), 'p1': (3.19, 3.39), 'p2': (4.64, 4.84)},
              4: {'f': (1.78, 1.87), 'p1': (3.59, 3.74), 'p2': (5.10, 5.26)},
              5: {'f': (2.01, 2.08), 'p1': (3.94, 4.06), 'p2': (5.52, 5.65)}}}

# ---------------------------------------------------------------- ecuación de estado
# p = K rho^2 con rho la densidad de masa en reposo ('rest') o la de energía ('energy').
_x = sp.Symbol('x', positive=True)
EOS = {'rest': _x + _x**2, 'energy': _x / (1 - _x)}     # eps(n) con m_b = 1 y K = 1


def _eos_num(kind):
    e = EOS[kind]
    p = sp.expand(_x * sp.diff(e, _x) - e)
    f = sp.lambdify(_x, [e, p, sp.diff(p, _x), sp.diff(p, _x, 2), sp.diff(e, _x)], 'numpy')
    n_of_p = (lambda q: np.sqrt(q)) if kind == 'rest' else (lambda q: np.sqrt(q) / (1 + np.sqrt(q)))
    return f, n_of_p


# ---------------------------------------------------------------- de sympy a una matriz numérica
_B = sp.symbols('bn bn1 bn2 bP bP1 bP2 bL bL1')     # n, n', n'', Phi, Phi', Phi'', lambda, lambda'


def _numeric(e, kind):
    """Sin axión, con la EOS concreta, y las funciones de fondo cambiadas por símbolos."""
    e = e.subs(dth, 0).subs(a0, 0).doit()
    e = e.replace(lambda x: isinstance(x, sp.Subs), lambda x: 0)
    e = e.subs({mstar(0): mN, V(0): 0}).subs(eps0, sp.Lambda(_x, EOS[kind])).doit()
    bn, bn1, bn2, bP, bP1, bP2, bL, bL1 = _B
    for k, v in [(d(Phi0, 2), bP2), (d(n0, 2), bn2), (d(Phi0), bP1), (d(lam0), bL1), (d(n0), bn1)]:
        e = e.subs(k, v)
    return e.subs({Phi0: bP, lam0: bL, n0: bn, G: 1})


def _matrix(eqs, ys, alg, kind):
    """y' = A(r, w, l, fondo) y para el sistema `eqs` (derivadas de `ys`, ligaduras `alg`)."""
    eqs = {k: _numeric(v, kind) for k, v in eqs.items()}
    sol = {}
    for name, var in alg:                       # ligaduras algebraicas, en orden
        sol[var] = sp.solve(eqs[name].subs(sol), var)[0]
    ders = [d(y) for y in ys]
    rest = [eqs[k].subs(sol) for k in eqs if k not in dict(alg)]
    rhs = sp.solve(rest, ders, dict=True)[0]
    A = sp.Matrix([[sp.factor(sp.cancel(sp.diff(rhs[dy], y))) for y in ys] for dy in ders])
    return sp.lambdify((r, w, l, *_B), A, 'numpy')


def full_system(kind, which='EXT'):
    """LD en variables (H1, K, W, X). 'EXT': el sistema extendido de la página sin axión.
    'ZL': Zhao y Lattimer tal como está impreso, que solo difiere en el signo de 2 w^2 e^{lambda-2Phi} W
    dentro de X' (paso 7); se construye sumando ese término, sin rehacer el álgebra."""
    A = _matrix({k: F.EXT[k] for k in ("H1'", "K'", "W'", "X'", 'H0', 'V')}, [H1, K, W, Xf],
                [('H0', H0), ('V', Vt)], kind)
    if which == 'EXT':
        return A
    eos = _eos_num(kind)[0]

    def A_zl(rr, omega, ll, n, n1, n2, P, P1, P2, L, L1):
        M = np.array(A(rr, omega, ll, n, n1, n2, P, P1, P2, L, L1), dtype=float)
        e, p = eos(n)[:2]
        M[3, 2] += 2 * omega**2 * (e + p) * np.exp(L - P) / rr      # r X' = ... + (e+p) e^Phi w^2 e^{lambda-2Phi} W
        return M
    return A_zl


def cowling_system(kind):
    """(8a)-(8c) del Paper II sin axión, en variables (W, X) con X = -e^Phi Delta p."""
    dn_of = sp.solve(Xf - F.X_def, dn)[0]
    eqs = {k: NR.PAPER[k].subs(dn, dn_of).doit() for k in ('(8a) momento radial', '(8b) bariones',
                                                         '(8c) momento angular')}
    return _matrix(eqs, [W, Xf], [('(8c) momento angular', Vt)], kind)


# ---------------------------------------------------------------- fondo (TOV)
class Star:
    def __init__(self, kind, compactness=None, mass=None):
        """Por compacidad M/R o por masa M (en unidades de K^{1/2}, rama estable)."""
        self.kind = kind
        self.eos, self.n_of_p = _eos_num(kind)
        lo, hi = 1e-4, 0.3 if kind == 'rest' else 0.2
        target = (lambda nc: self._solve(nc)[0] - compactness) if mass is None else (lambda nc: self._solve(nc)[1] - mass)
        self.nc = brentq(target, lo, hi if mass is None else 0.15, xtol=1e-14)
        self.MR, self.M, self.R, self.Phic = self._solve(self.nc)

    def bg(self, rr, m, p, Phi):
        """(n, n', n'', Phi, Phi', Phi'', lambda, lambda') y (m', p', Phi')."""
        n = self.n_of_p(max(p, 0.0))
        e, _, dpdn, d2p, dedn = self.eos(n)
        D = rr * (rr - 2 * m)
        Pp = (m + 4 * np.pi * rr**3 * p) / D
        mp, pp = 4 * np.pi * rr**2 * e, -(e + p) * Pp
        Npr = mp + 12 * np.pi * rr**2 * p + 4 * np.pi * rr**3 * pp
        Dp = 2 * rr - 2 * m - 2 * rr * mp
        Ppp = (Npr * D - (m + 4 * np.pi * rr**3 * p) * Dp) / D**2
        lam, lamp = -0.5 * np.log(1 - 2 * m / rr), (rr * mp - m) / D
        n1 = pp / dpdn if dpdn > 0 else 0.0
        ppp = -(dedn * n1 + pp) * Pp - (e + p) * Ppp
        n2 = (ppp - d2p * n1**2) / dpdn if dpdn > 0 else 0.0
        return (n, n1, n2, Phi, Pp, Ppp, lam, lamp), (mp, pp, Pp)

    def _start(self, nc, r0):
        e, p = self.eos(nc)[:2]
        return [4 * np.pi / 3 * e * r0**3, p - 2 * np.pi / 3 * (e + p) * (e + 3 * p) * r0**2, 0.0]

    def _solve(self, nc):
        r0 = 1e-6
        surf = lambda rr, y: y[1] - 1e-14 * self.eos(nc)[1]
        surf.terminal = True
        s = solve_ivp(lambda rr, y: self.bg(rr, *y)[1], (r0, 50), self._start(nc, r0), events=surf,
                      method='DOP853', rtol=1e-12, atol=1e-30)
        R = s.t[-1]
        M = s.y[0, -1]
        Phic = 0.5 * np.log(1 - 2 * M / R) - s.y[2, -1]
        return M / R, M, R, Phic

    def integrate(self, A, y0, omega, ll, r0=None, r1=None):
        """Integra fondo + perturbaciones desde r0 hasta r1 (cerca de la superficie)."""
        r0 = r0 or 1e-4 * self.R
        r1 = r1 or self.R * (1 - 1e-7)
        b0 = self._start(self.nc, r0)
        b0[2] = self.Phic

        def rhs(rr, y):
            B, dB = self.bg(rr, *y[:3])
            return np.concatenate([dB, np.asarray(A(rr, omega, ll, *B), dtype=float) @ y[3:]])
        s = solve_ivp(rhs, (r0, r1), np.concatenate([b0, y0]), method='DOP853', rtol=1e-11, atol=1e-14)
        return s.y[:, -1]

    def center(self, omega, ll, K0, W0):
        """Valores en el centro (LD): H1(0), K(0), W(0), X(0)."""
        e, p = self.eos(self.nc)[:2]
        X0 = (e + p) * np.exp(self.Phic) * ((4 * np.pi / 3 * (e + 3 * p)
                                             - omega**2 * np.exp(-2 * self.Phic) / ll) * W0 + K0 / 2)
        return np.array([(2 * ll * K0 + 16 * np.pi * (e + p) * W0) / (ll * (ll + 1)), K0, W0, X0])


# ---------------------------------------------------------------- frecuencias
def cowling_X(star, A, omega, ll):
    y0 = star.center(omega, ll, 0.0, 1.0)[2:]
    return star.integrate(A, y0, omega, ll)[-1]


def cowling_modes(star, A, ll, grid):
    vals = [cowling_X(star, A, o, ll) for o in grid]
    roots = []
    for a, b, fa_, fb in zip(grid[:-1], grid[1:], vals[:-1], vals[1:]):
        if fa_ * fb < 0:
            roots.append(brentq(lambda o: cowling_X(star, A, o, ll), a, b, xtol=1e-12))
    return roots


def outgoing_amplitude(star, A, omega, ll, n_far=40):
    """|A_out|^2 de la onda de Zerilli fuera, para la solución interior con X(R) = 0 y W(R) = 1.
    A frecuencia real el exterior es una onda estacionaria; su amplitud es mínima en la resonancia."""
    ends = []
    for K0, W0 in [(1.0, 0.0), (0.0, 1.0)]:
        ends.append(star.integrate(A, star.center(omega, ll, K0, W0), omega, ll))
    y1, y2 = ends[0][3:], ends[1][3:]
    c = -y1[3] / y2[3]                              # X(R) = 0
    y = y1 + c * y2
    y = y / y[2]                                     # W(R) = 1
    M, R = star.M, star.R
    rr = R * (1 - 1e-7)
    B, _ = star.bg(rr, *ends[0][:3])
    dy = np.asarray(A(rr, omega, ll, *B), dtype=float) @ y
    nn = (ll - 1) * (ll + 2) / 2
    f = lambda x: 1 - 2 * M / x
    h1, k, dh1, dk = y[0], y[1], dy[0], dy[1]
    pref = rr**(ll + 2) / (nn * rr + 3 * M)
    dpref = pref * ((ll + 2) / rr - nn / (nn * rr + 3 * M))
    Z = pref * (k - f(rr) * h1)
    dZ = dpref * (k - f(rr) * h1) + pref * (dk - 2 * M / rr**2 * h1 - f(rr) * dh1)
    VZ = lambda x: f(x) * (2 * nn**2 * (nn + 1) * x**3 + 6 * nn**2 * M * x**2 + 18 * nn * M**2 * x + 18 * M**3) \
        / (x**3 * (nn * x + 3 * M)**2)
    r_far = R + n_far * 2 * np.pi / omega
    s = solve_ivp(lambda x, u: [u[1] / f(x), -(omega**2 - VZ(x)) * u[0] / f(x)], (R, r_far), [Z, f(R) * dZ],
                  method='DOP853', rtol=1e-11, atol=1e-14)
    Zf, Pf = s.y[0, -1], s.y[1, -1]                 # Z y dZ/dr*
    L2 = ll * (ll + 1)
    rs = r_far + 2 * M * np.log(r_far / (2 * M) - 1)
    u = np.exp(1j * omega * rs) * (1 + 1j * L2 / (2 * omega * r_far))
    du = 1j * omega * u + f(r_far) * np.exp(1j * omega * rs) * (-1j * L2 / (2 * omega * r_far**2))
    # Z = 2 Re(a u), dZ/dr* = 2 Re(a du): dos ecuaciones reales para a complejo
    Mx = np.array([[2 * u.real, -2 * u.imag], [2 * du.real, -2 * du.imag]])
    ar, ai = np.linalg.solve(Mx, [Zf, Pf])
    return ar**2 + ai**2


def full_mode(star, A, ll, guess, width=0.04):
    """Parte real (mínimo de |A_out|^2) e imaginaria (ajuste parabólico |A|^2 ~ (w-w0)^2 + g^2)."""
    unit = np.sqrt(star.M / star.R**3)
    amp = lambda s: outgoing_amplitude(star, A, s * unit, ll)
    m = minimize_scalar(amp, bracket=(guess * (1 - width), guess, guess * (1 + width)), tol=1e-10)
    s0 = m.x
    h = 1e-3 * s0
    a0_, ap, am = amp(s0), amp(s0 + h), amp(s0 - h)
    curv = (ap + am - 2 * a0_) / (2 * h**2)
    gamma = np.sqrt(max(a0_ / curv, 0.0))
    return s0, gamma
