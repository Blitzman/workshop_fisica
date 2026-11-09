"""Etapa B: el sistema extendido (fullgr.EXT, Lindblom-Detweiler + axión) integrado con el axión encendido.

El sistema no se reescribe a mano: las ecuaciones de fullgr.EXT, derivadas desde la acción en el paso 7,
se convierten en una matriz numérica y' = A y sustituyendo el fondo por símbolos:

  y = (H1, K, W, X, dtheta, dtheta')            dentro;   H0 y V algebraicos
  y = (H1, K, dtheta, dtheta')                  fuera (sin fluido; la cola del axión sigue acoplada)

Unidades geométricas en km. La densidad bariónica va en unidades de n_c (x = n_b/n_c), de modo que
n_b dm*/dtheta = x E0 g'(theta) y n_b (m* - m_N) = -x E0 (1 - g), con E0 = eps m_pi^2 f_pi^2.
"""
from functools import lru_cache

import numpy as np
import sympy as sp

from . import fullgr as F
from . import nonradial as NR
from .model import r, Phi0, lam0, a0, n0, w, G, V, mstar, mN, eps0

H0, H1, K, W, Vt, dn, dth, Xf = F.H0, F.H1, F.K, F.W, F.Vt, F.dn, F.dth, F.Xf
d, l = F.d, NR.l
fa = NR.fa

# símbolos del fondo, en este orden en las funciones numéricas
BG = sp.symbols('bP bP1 bP2 bL bL1 x x1 th th1 e0 pr cs0 fa Mst M1 M2 Vv V2')
(bP, bP1, bP2, bL, bL1, x, x1, th, th1, e0, pr, cs0, faS, Mst, M1, M2, Vv, V2) = BG


def to_numeric(e):
    """Cambia las funciones del fondo por símbolos. eps0'(n) = (eps0 + p)/n (primera ley) y
    eps0''(n) = (dp/dn)/n con dp/dn = (eps0 + p) c_s0^2 / n; c_s0^2 = dp/deps0 de la EOS sin axión."""
    e = e.doit()
    dpdn = (e0 + pr) * cs0 / x
    rep = [(sp.Derivative(eps0(n0), (n0, 2)), dpdn / x), (sp.Derivative(eps0(n0), n0), (e0 + pr) / x),
           (eps0(n0), e0),
           (sp.Derivative(mstar(a0), (a0, 2)), M2 / faS**2), (sp.Derivative(mstar(a0), a0), M1 / faS),
           (mstar(a0), mN + Mst),
           (sp.Derivative(V(a0), (a0, 2)), V2 / faS**2), (V(a0), Vv),
           (sp.Derivative(Phi0, (r, 2)), bP2), (sp.Derivative(Phi0, r), bP1), (sp.Derivative(lam0, r), bL1),
           (sp.Derivative(a0, r), faS * th1), (sp.Derivative(n0, r), x1)]
    for k, v in rep:
        e = e.subs(k, v)
    e = e.subs({Phi0: bP, lam0: bL, a0: faS * th, n0: x, G: 1, fa: faS})
    return e


def _linear(eqs, ys, dys, alg):
    """Matriz A con dys = A ys a partir de ecuaciones lineales; `alg` = [(eq, var)] ligaduras en orden."""
    sol = {}
    for eq, var in alg:
        sol[var] = sp.solve(eq.subs(sol), var)[0]
    rows = []
    for eq, dy in zip(eqs, dys):
        if isinstance(eq, int):                  # fila identidad: (y_j)' = y_eq
            rows.append([1 if k == eq else 0 for k in range(len(ys))])
            continue
        eq = sp.expand(eq.subs(sol))
        c = eq.coeff(dy)
        rest = sp.expand(eq - c * dy)
        rows.append([sp.expand(-rest.coeff(y) / c) for y in ys])
        left = sp.expand(rest - sum(rest.coeff(y) * y for y in ys))
        if left != 0:
            raise ValueError(f'término no lineal o variable sin resolver: {left.free_symbols}')
    return sp.Matrix(rows)


def _states(funcs):
    """Funciones -> símbolos, para poder extraer coeficientes."""
    return {f: sp.Symbol(f'S{i}') for i, f in enumerate(funcs)}


@lru_cache
def interior_equations():
    dn_of = sp.solve(Xf - F.X_def, dn)[0]
    E = {k: to_numeric(v.subs(dn, dn_of)) for k, v in F.EXT.items()}
    S = _states([d(dth, 2), d(H1), d(K), d(W), d(Xf), d(dth), H0, Vt, H1, K, W, Xf, dth])
    return {k: v.subs(S) for k, v in E.items()}, S


@lru_cache
def interior_matrix():
    """A(r, w, l, fondo) del sistema extendido dentro de la estrella (6 x 6)."""
    E, S = interior_equations()
    ys = [S[H1], S[K], S[W], S[Xf], S[dth], S[d(dth)]]
    dys = [S[d(H1)], S[d(K)], S[d(W)], S[d(Xf)], S[d(dth)], S[d(dth, 2)]]
    eqs = [E["H1'"], E["K'"], E["W'"], E["X'"], 5, E['KG']]
    A = _linear(eqs, ys, dys, [(E['H0'], S[H0]), (E['V'], S[Vt])])
    return sp.lambdify((r, w, l, *BG), A, 'numpy', cse=True)


@lru_cache
def exterior_matrix():
    """Fuera: sin fluido (x = e0 = p = 0, W = X = V = 0); quedan la métrica y el axión (4 x 4)."""
    zero = {x: 0, e0: 0, pr: 0, x1: 0}
    E = {}
    for k in ("H1'", "K'", 'KG', 'H0'):
        e = to_numeric(F.EXT[k].subs({W: 0, Xf: 0, Vt: 0, dn: 0}).doit())
        E[k] = sp.expand(e.subs(zero))
    S = _states([d(dth, 2), d(H1), d(K), d(dth), H0, H1, K, dth])
    E = {k: v.subs(S) for k, v in E.items()}
    ys = [S[H1], S[K], S[dth], S[d(dth)]]
    dys = [S[d(H1)], S[d(K)], S[d(dth)], S[d(dth, 2)]]
    eqs = [E["H1'"], E["K'"], 3, E['KG']]
    A = _linear(eqs, ys, dys, [(E['H0'], S[H0])])
    return sp.lambdify((r, w, l, *BG), A, 'numpy', cse=True)


def substitution_checks():
    """Comprobaciones de las sustituciones del axión: cantidades del paper II, convertidas, frente a lo esperado."""
    eps = e0 + x * Mst
    cases = [('n dm*/dtheta', n0 * NR.dm, x * M1),
             ('n d2m*/dtheta2', n0 * NR.dmm, x * M2),
             ('m_eff^2', NR.meff2, (V2 + x * M2) / faS**2),
             ('eps(n, theta)', NR.eps, eps), ('p(n)', NR.p, pr),
             ('c_s^2', NR.cs2, (e0 + pr) * cs0 / (eps + pr)),
             ("theta_0'", F.thp, th1), ('V(a_0)', V(a0), Vv)]
    return [(name, sp.simplify(to_numeric(e) - want) == 0) for name, e, want in cases]


# ---------------------------------------------------------------- fondo de axstar como valores numéricos
def bgvals(st, rr, h=1e-5):
    """Valores de BG en rr para una axstar.Star (dentro o fuera)."""
    from . import axstar as AX
    m, p, Phi, t, tp = st.sol.sol(rr)
    der = st.rhs(rr, [m, p, Phi, t, tp])
    dm, dp, dPhi = der[0], der[1], der[2]
    y2, y1 = st.sol.sol(rr + h), st.sol.sol(rr - h)
    P2 = (st.rhs(rr + h, list(y2))[2] - st.rhs(rr - h, list(y1))[2]) / (2 * h)
    lam = -0.5 * np.log(1 - 2 * m / rr)
    dlam = (dm * rr - m) / (rr * (rr - 2 * m))
    pm = p / AX.KGEO
    inside = pm > AX._eos.p_min * 1.0001
    if inside:
        e0v, xv, cs0v = AX._eos.eps(pm) * AX.KGEO, AX._eos.nb(pm) / st.nc, AX._eos.cs2(pm)
        x1v = dp / ((e0v + p) * cs0v / xv)
    else:
        e0v, xv, cs0v, x1v, p = 0.0, 0.0, 1.0, 0.0, 0.0
    if st.nu_a:
        E0, fav = st.E0, np.sqrt(st.fa2)
        G, G1, G2 = AX.g(t), AX.g1(t), AX.g2(t)
        ax = (fav, -E0 * (1 - G), E0 * G1, E0 * G2, E0 * (1 - G), -E0 * G2)
    else:
        t = tp = 0.0
        ax = (1.0, 0.0, 0.0, 0.0, 0.0, 0.0)
    return (Phi + st.Phi_shift, dPhi, P2, lam, dlam, xv, x1v, t, tp, e0v, p, cs0v, *ax)


# ---------------------------------------------------------------- modos con métrica dinámica y axión
def center(st, omega, ll, K0, W0, T0, r0):
    """Condiciones en el centro del sistema extendido (tercera parte de la página): LD más las
    correcciones axiónicas en X(0) y delta theta''(0) (a orden dominante: delta theta(r0) = delta theta(0))."""
    from . import axstar as AX
    b = bgvals(st, r0)
    P, xv, t, e0v, p, M1v, Vvv = b[0], b[5], b[7], b[9], b[10], b[14], b[16]
    eps = e0v + xv * b[13]
    tpp0 = st._start(r0, st.thc)[4] / r0 if st.nu_a else 0.0
    H1 = (2 * ll * K0 + 16 * np.pi * (eps + p) * W0) / (ll * (ll + 1))
    X = (eps + p) * np.exp(P) * ((4 * np.pi / 3 * (eps + 3 * p) - omega**2 * np.exp(-2 * P) / ll) * W0 + K0 / 2)
    if st.nu_a:
        X += np.exp(P) * (xv * M1v * (T0 + tpp0 * W0) - 8 * np.pi / 3 * (eps + p) * Vvv * W0)
    return np.array([H1, K0, W0, X, T0, 0.0], dtype=complex)


def _rhs(st, A, omega, ll, keep=None):
    """keep: índices de las variables que evolucionan (las demás, delta theta, quedan en 0)."""
    if keep is None:
        return lambda rr, y: np.asarray(A(rr, omega, ll, *bgvals(st, rr)), dtype=complex) @ y

    def f(rr, y):
        M = np.asarray(A(rr, omega, ll, *bgvals(st, rr)), dtype=complex)
        out = np.zeros_like(y)
        out[keep] = M[np.ix_(keep, keep)] @ y[keep]
        return out
    return f


ATOL0 = 1e-15


def basis(st, omega, ll, dyn=True, r0=1e-3, r1=0.5, rtol=1e-10, rm=None):
    """Soluciones regulares desde el centro hasta R, y su continuación fuera hasta r_m.
    dyn=False: delta theta apagado (solo el fondo con axión, como la Tabla I del Paper II)."""
    from scipy.integrate import solve_ivp
    Ai, Ae = interior_matrix(), exterior_matrix()
    dyn = dyn and bool(st.nu_a)
    seeds = [(1.0, 0.0, 0.0), (0.0, 1.0, 0.0)] + ([(0.0, 0.0, 1.0)] if dyn else [])
    ki, ke = (None, None) if dyn else ([0, 1, 2, 3], [0, 1])
    R, rm = st.R, rm or r_match(st)
    out = []
    for K0, W0, T0 in seeds:
        f = _rhs(st, Ai, omega, ll, ki)
        # H29: BDF fallaba en su primer paso y su resultado se usaba sin mirarlo. Cerca del centro hay componentes
        # que arrancan en 0: con atol = 1e-40 cualquier integrador se ahoga ahí; ATOL0 es absoluto sobre semillas O(1)
        s = _ok(solve_ivp(f, (r0, r1), center(st, omega, ll, K0, W0, T0, r0), method='DOP853', rtol=rtol, atol=ATOL0))
        # Mismo mecanismo en 0,5 km -> R: con atol = 1e-40 la semilla axiónica a nu_a >= 13 kHz tardaba de minutos a
        # horas por evaluación. Con atol relativo a la escala de la solución las raíces del paso 10 a 6, 10 y 12 kHz
        # cambian <= 7e-9 relativo (cuaderno, H29: segundo atasco)
        y1 = s.y[:, -1]
        s = _ok(solve_ivp(f, (r1, R), y1, method='DOP853', rtol=rtol, atol=ATOL0 * np.max(np.abs(y1))))
        yR = s.y[:, -1]
        ye = np.array([yR[0], yR[1], yR[4], yR[5]])
        if rm > R:
            s = _ok(solve_ivp(_rhs(st, Ae, omega, ll, ke), (R, rm), ye, method='DOP853', rtol=rtol, atol=1e-40))
            ye = s.y[:, -1]
        out.append((yR[3], ye))
    return out


def _ok(sol):
    """H29: si solve_ivp no llega al final, se lanza (antes se tomaba s.y[:, -1] sin mirar status)."""
    if sol.status < 0:
        raise RuntimeError(f'solve_ivp: {sol.message}')
    return sol


def r_match(st):
    """Radio donde se empalma con el exterior desacoplado: donde la cola de theta_0 ya es despreciable."""
    if not st.nu_a:
        return st.R
    if not hasattr(st, '_rm'):
        rr = np.linspace(st.R, st.sol.t[-1], 2000)
        th = np.abs(st.sol.sol(rr)[3])
        i = np.argmax(th < 1e-7)
        st._rm = rr[i] if i else rr[np.argmin(th)]
    return st._rm


def zerilli_amplitudes(st, omega, ll, ye, rm, n_far=40):
    """(A_out, A_in) de la onda de Zerilli construida con (H1, K) en rm (fuera, sin axión de fondo)."""
    from scipy.integrate import solve_ivp
    Ae = exterior_matrix()
    dy = np.asarray(Ae(rm, omega, ll, *bgvals(st, rm)), dtype=complex) @ ye
    M = st.sol.sol(rm)[0] if rm < st.sol.t[-1] else st.sol.y[0, -1]
    nn = (ll - 1) * (ll + 2) / 2
    f = lambda x_: 1 - 2 * M / x_
    h1, k, dh1, dk = ye[0], ye[1], dy[0], dy[1]
    pref = rm**(ll + 2) / (nn * rm + 3 * M)
    dpref = pref * ((ll + 2) / rm - nn / (nn * rm + 3 * M))
    Z = pref * (k - f(rm) * h1)
    dZ = dpref * (k - f(rm) * h1) + pref * (dk - 2 * M / rm**2 * h1 - f(rm) * dh1)
    VZ = lambda x_: f(x_) * (2 * nn**2 * (nn + 1) * x_**3 + 6 * nn**2 * M * x_**2 + 18 * nn * M**2 * x_ + 18 * M**3) \
        / (x_**3 * (nn * x_ + 3 * M)**2)
    r_far = rm + n_far * 2 * np.pi / abs(omega.real)
    s = _ok(solve_ivp(lambda x_, u: [u[1] / f(x_), -(omega**2 - VZ(x_)) * u[0] / f(x_)], (rm, r_far),
                      np.array([Z, f(rm) * dZ], dtype=complex), method='DOP853', rtol=1e-11, atol=1e-60))   # H32
    Zf, Pf = s.y[0, -1], s.y[1, -1]
    L2 = ll * (ll + 1)
    rs = r_far + 2 * M * np.log(r_far / (2 * M) - 1)
    out = []
    for sg in (1, -1):        # u = saliente, v = entrante (continuaciones analíticas en omega)
        e = np.exp(sg * 1j * omega * rs)
        u = e * (1 + sg * 1j * L2 / (2 * omega * r_far))
        du = sg * 1j * omega * u + f(r_far) * e * (-sg * 1j * L2 / (2 * omega * r_far**2))
        out.append((u, du))
    (u, du), (v, dv) = out
    a_out, a_in = np.linalg.solve(np.array([[u, v], [du, dv]]), np.array([Zf, Pf]))
    return a_out, a_in


def incoming(st, omega, ll, dyn=True, rtol=1e-10, rm=None, n_far=40):
    """Amplitud entrante de la onda gravitacional para la solución que cumple Delta p(R) = 0 y la
    condición del axión fuera. Es analítica en omega: sus ceros son los modos cuasinormales."""
    from . import axstar as AX
    rm = rm or r_match(st)
    B = basis(st, omega, ll, dyn=dyn, rtol=rtol, rm=rm)
    c1 = np.array([b[0] for b in B])                              # X(R)
    if len(B) == 3:
        L = AX.exterior_L(st, ll, omega, r_stop=rm) - ll / rm
        c2 = np.array([b[1][3] - L * b[1][2] for b in B])         # delta theta' - L delta theta en rm
        c = np.cross(c1, c2)
    else:
        c = np.array([c1[1], -c1[0]])
    ye = sum(ci * b[1] for ci, b in zip(c, B))
    return zerilli_amplitudes(st, omega, ll, ye, rm, n_far)[1]
