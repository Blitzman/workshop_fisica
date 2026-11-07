"""Perturbaciones no radiales polares con la métrica perturbada (más allá de Cowling).

Trabajo futuro que plantea el Paper II (arXiv:2610.37308, Sec. IV). Misma acción que model.py;
ahora también se perturba la métrica, en el gauge de Regge-Wheeler y con la convención de
Lindblom y Detweiler (LD), llevada a signatura (+,-,-,-) y dependencia e^{-i w t}:

    ds^2 = e^{2 Phi}(1 + H0 rY) dt^2 - 2 i w r H1 rY dt dr - e^{2 lambda}(1 - H2 rY) dr^2
           - r^2 (1 - K rY) dOmega^2,          rY = (r/R)^l Y(th) e^{-i w t}

El fluido y el axión siguen la convención del Paper II (nonradial.py), que coincide con la de LD
para W y V. Con H0 = H1 = H2 = K = 0 se recupera exactamente el caso de Cowling.

La geometría se linealiza a mano (Gamma, Ricci y Einstein a primer orden) para no calcular el
tensor de Einstein completo de una métrica que depende de (t, r, th).
"""
import sympy as sp

from .model import r, th, t, V, eps_of, p_of, Phi0, lam0, a0, n0, w, G
from .nonradial import X, l, Rs, fa, A, N, XR, XT, dth, dn, W, Vt, Y, _legendre, _divide

H0, H1, H2, K = [sp.Function(s)(r) for s in ('H0', 'H1', 'H2', 'K')]
eta = sp.Symbol('eta')
E = sp.exp(-sp.I * w * t)
rl = (r / Rs)**l

g0 = sp.diag(sp.exp(2 * Phi0), -sp.exp(2 * lam0), -r**2, -r**2 * sp.sin(th)**2)
g0inv = sp.diag(sp.exp(-2 * Phi0), -sp.exp(-2 * lam0), -1 / r**2, -1 / (r**2 * sp.sin(th)**2))
h = sp.zeros(4, 4)
h[0, 0] = sp.exp(2 * Phi0) * H0
h[0, 1] = h[1, 0] = -sp.I * w * r * H1
h[1, 1] = sp.exp(2 * lam0) * H2
h[2, 2] = r**2 * K
h[3, 3] = r**2 * sp.sin(th)**2 * K
h = h * rl * Y * E
g1inv = -g0inv * h * g0inv                      # g^{mu nu} a primer orden
sqrtg0 = sp.exp(Phi0 + lam0) * r**2 * sp.sin(th)
sqrtg1 = sqrtg0 * sum(g0inv[i, i] * h[i, i] for i in range(4)) / 2


def _gamma(gup, gdn_list):
    """Gamma^a_bc = 1/2 g^{ad}(d_b g_dc + d_c g_db - d_d g_bc), con g^{ad} diagonal o no."""
    Gam = [[[0] * 4 for _ in range(4)] for _ in range(4)]
    for a in range(4):
        for b in range(4):
            for c in range(b, 4):
                s = 0
                for gu, gd in gdn_list(gup):
                    s += sum(gu[a, d_] * (sp.diff(gd[d_, c], X[b]) + sp.diff(gd[d_, b], X[c])
                                          - sp.diff(gd[b, c], X[d_])) for d_ in range(4)) / 2
                Gam[a][b][c] = Gam[a][c][b] = sp.expand(s)
    return Gam


def _ricci(G0, G1=None):
    """R_bd = d_a Gam^a_db - d_d Gam^a_ab + Gam^a_ae Gam^e_db - Gam^a_de Gam^e_ab (parte lineal si G1)."""
    Ric = sp.zeros(4, 4)
    for b in range(4):
        for d_ in range(b, 4):
            if G1 is None:
                s = sum(sp.diff(G0[a][d_][b], X[a]) - sp.diff(G0[a][a][b], X[d_]) for a in range(4))
                s += sum(G0[a][a][e] * G0[e][d_][b] - G0[a][d_][e] * G0[e][a][b]
                         for a in range(4) for e in range(4))
            else:
                s = sum(sp.diff(G1[a][d_][b], X[a]) - sp.diff(G1[a][a][b], X[d_]) for a in range(4))
                s += sum(G0[a][a][e] * G1[e][d_][b] + G1[a][a][e] * G0[e][d_][b]
                         - G0[a][d_][e] * G1[e][a][b] - G1[a][d_][e] * G0[e][a][b]
                         for a in range(4) for e in range(4))
            Ric[b, d_] = Ric[d_, b] = sp.expand(s)
    return Ric


def linear_einstein():
    """G_{mu nu} a primer orden: (G0, G1)."""
    Gam0 = _gamma(None, lambda _: [(g0inv, g0)])
    Gam1 = _gamma(None, lambda _: [(g0inv, h), (g1inv, g0)])
    Ric0, Ric1 = _ricci(Gam0), _ricci(Gam0, Gam1)
    R0 = sum(g0inv[i, i] * Ric0[i, i] for i in range(4))
    R1 = sum(g0inv[i, i] * Ric1[i, i] for i in range(4)) + \
        sum(g1inv[i, j] * Ric0[i, j] for i in range(4) for j in range(4))
    G0 = Ric0 - R0 * g0 / 2
    G1 = Ric1 - (R1 * g0 + R0 * h) / 2
    return G0, G1, Gam0, Gam1


def derive():
    """Ecuaciones de campo con eta marcando el orden en la métrica; los campos materiales son
    generales (t, r, th) y se linealizan después con `linearize`."""
    G0, G1, Gam0, Gam1 = linear_einstein()
    g = g0 + eta * h
    ginv = g0inv + eta * g1inv
    sqrtg = sqrtg0 + eta * sqrtg1
    Gam = [[[Gam0[a][b][c] + eta * Gam1[a][b][c] for c in range(4)] for b in range(4)] for a in range(4)]
    dA = [sp.diff(A, x) for x in X]

    # Klein-Gordon por Euler-Lagrange (n fijo); normalización box a + V' + d eps/da
    L_a = sp.Rational(1, 2) * sum(ginv[i, j] * dA[i] * dA[j] for i in range(4) for j in range(4)) \
        - V(A) - eps_of(N, A)
    KG = -sp.euler_equations(sqrtg * L_a, A, [t, r, th])[0].lhs / sqrtg0

    # T_{mu nu} del axión por variación respecto a g^{mu nu}
    hh = sp.Matrix(4, 4, lambda i, j: sp.Symbol(f'h{i}{j}'))
    L_h = sp.Rational(1, 2) * sum(hh[i, j] * dA[i] * dA[j] for i in range(4) for j in range(4)) - V(A)
    L_val = L_h.subs({hh[i, j]: ginv[i, j] for i in range(4) for j in range(4)})
    T_ax = sp.Matrix(4, 4, lambda i, j: 2 * sp.diff(L_h, hh[i, j]) - g[i, j] * L_val)

    # fluido perfecto: u^mu = u^t (1, d_t xi^r, d_t xi^th, 0), normalizada con la métrica perturbada
    v = sp.Matrix([1, sp.diff(XR, t), sp.diff(XT, t), 0])
    ut = 1 / sp.sqrt((v.T * g * v)[0, 0])
    u_up = ut * v
    u_dn = g * u_up
    e, p = eps_of(N, A), p_of(N)
    T_fl = sp.Matrix(4, 4, lambda i, j: (e + p) * u_dn[i] * u_dn[j] - p * g[i, j])
    T = T_fl + T_ax

    out = {'KG': KG}
    for (i, j), name in [((0, 0), 'Ett'), ((0, 1), 'Etr'), ((0, 2), 'Etth'), ((1, 1), 'Err'),
                         ((1, 2), 'Erth'), ((2, 2), 'Ethth'), ((3, 3), 'Ephph')]:
        out[name] = (G0[i, j] + eta * G1[i, j]) - 8 * sp.pi * G * T[i, j]

    def div(nu):  # nabla_mu T^mu_nu con la conexión a primer orden
        M = ginv * T
        s = sum(sp.diff(sqrtg * M[m, nu], X[m]) for m in range(4)) / sqrtg
        return s - sum(Gam[a][m][nu] * M[m, a] for m in range(4) for a in range(4))

    for nu, name in [(1, 'FLUID_r'), (2, 'FLUID_th')]:
        out[name] = div(nu) - KG * sp.diff(A, X[nu])
    out['BAR'] = sum(sp.diff(sqrtg * N * u_up[m], X[m]) for m in range(4)) / (r**2 * sp.sin(th))
    return out


# factores angulares tras eliminar Y'' con Legendre
Yp = sp.diff(Y, th)
ANG = {'Y': Y, 'dY': Yp, 'traceless': -2 * sp.cos(th) / sp.sin(th) * Yp - l * (l + 1) * Y}


def linearize(expr, angular):
    """Parte lineal (métrica y materia), dividida por e^{-i w t} (r/R)^l y el factor angular."""
    sub = {A: a0 + eta * fa * rl * dth * Y * E,
           N: n0 + eta * rl * dn * Y * E,
           XR: eta * rl * sp.exp(-lam0) * W / r * Y * E,
           XT: -eta * rl * Vt / r**2 * Yp * E}
    lin = sp.diff(expr.subs(sub).doit(), eta).subs(eta, 0)
    lin = _legendre(sp.expand(lin / E / rl))
    lin = sp.expand(sp.powdenest(sp.powsimp(sp.expand(_divide(lin, ANG[angular])), force=True), force=True))
    if lin.has(th):
        raise ValueError('la dependencia angular no se factoriza')
    return lin


def linear_system():
    """Todas las ecuaciones linealizadas, con su factor angular."""
    D = derive()
    out = {}
    for k, ang in [('Ett', 'Y'), ('Etr', 'Y'), ('Err', 'Y'), ('Etth', 'dY'), ('Erth', 'dY'),
                   ('KG', 'Y'), ('FLUID_r', 'Y'), ('FLUID_th', 'dY'), ('BAR', 'Y')]:
        out[k] = linearize(D[k], ang)
    s = sp.sin(th)**2
    out['Eang'] = linearize(D['Ethth'] + D['Ephph'] / s, 'Y')        # traza angular
    out['Etl'] = linearize(D['Ethth'] - D['Ephph'] / s, 'traceless')  # parte sin traza
    return out


# ---------------------------------------------------------------- Lindblom-Detweiler (sin axión)
# Transcritas del LaTeX de dos fuentes, con G explícita y en nuestras variables:
# nu = 2 Phi, lambda_LD = 2 lambda, X = -e^{Phi} Delta p, Delta p = delta p + xi^r p'.
from . import nonradial as NR  # noqa: E402

d = lambda f, k=1: sp.diff(f, r, k)
Xf = sp.Function('X')(r)
ep, p, cs2 = NR.eps + NR.p, NR.p, NR.cs2
psi, eL, eh = 2 * Phi0, sp.exp(2 * lam0), sp.exp(lam0)
nn = (l - 1) * (l + 2) / 2
pp = NR.dpdn * d(n0)
X_def = -sp.exp(Phi0) * (NR.dp + sp.exp(-lam0) * W / r * pp)
GM = r * (1 - sp.exp(-2 * lam0)) / 2
b_TOV = GM / r
Q_TOV = b_TOV + 4 * sp.pi * G * r**2 * p          # Q de LD (fluido)
Q_met = r * sp.exp(-2 * lam0) * d(Phi0)            # el mismo Q leído de la métrica


def lindblom_detweiler(Q, b, sw=1):
    """Ecuaciones de LD en la forma de Zhao y Lattimer (arXiv:2204.03037, Sec. II). `sw` = signo del
    término 2 w^2 W en X' (+1 tal como está impreso)."""
    return {
        "H1'": r * d(H1) + (l + 1 + 2 * b * eL + 4 * sp.pi * G * r**2 * eL * (p - NR.eps)) * H1
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
_h1 = lambda h1: d(H1) - (-(l + 1 + 2 * _M * eL / r + 4 * sp.pi * G * r**2 * eL * (p - NR.eps)) / r * h1
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

# ---------------------------------------------------------------- sistema extendido con axión
# LD (corregida) con b y Q leídos de la métrica, más los términos nuevos. Todos los términos
# nuevos se anulan sin axión (theta_0' = 0, dm*/dtheta = 0, V = 0, delta theta = 0).
_Pp = d(Phi0)
thp = d(a0) / fa
dm, dmm = NR.dm, NR.dmm
Dtheta = dth + sp.exp(-lam0) * thp * W / r          # perturbación lagrangiana de theta
_LDm = lindblom_detweiler(Q_met, b_TOV, sw=-1)
NEW = {   # términos nuevos, en la forma  ecuación_LD = NEW
    "H1'": 8 * sp.pi * G * r**2 * eL * V(a0) * H1,
    "K'": 8 * sp.pi * G * r * fa**2 * thp * dth,
    "W'": sp.S(0),
    "X'": sp.exp(Phi0) * (n0 * dm * r * thp * (H0 / 2 + K) - l * (l + 1) * n0 * dm * thp * Vt / r
                          + n0 * dm * r * d(dth)
                          + (n0 * dm * (l + r * _Pp) + n0 * dmm * r * thp - 4 * sp.pi * G * fa**2 * r * thp * ep) * dth
                          + n0 * sp.exp(-lam0) * (dm**2 * eL * n0 / fa**2 - 4 * dm * thp / r + dmm * thp**2) * W),
    'H0': -4 * sp.pi * G * r**2 * sp.exp(-2 * lam0) * d(a0)**2 * H0
    - 8 * sp.pi * G * r * sp.exp(-2 * lam0) * fa * d(a0) * (r * d(dth) + (l + 2 + r * _Pp) * dth)
    - 8 * sp.pi * G * r * sp.exp(-lam0) * n0 * dm * thp * W,
    'V': -n0 * dm * sp.exp(2 * Phi0) * Dtheta / (w**2 * ep),
    'KG': -(n0 * dm / fa**2 - 2 * sp.exp(-2 * lam0) * _Pp * thp) * H0
    + 16 * sp.pi * G * sp.exp(-2 * lam0) * fa**2 * thp**2 * dth,
}
EXT = {k: (NR.PAPER['(8d) KG interior'] if k == 'KG' else _LDm[k]) - v for k, v in NEW.items()}
