"""Modos radiales (l = 0) en relatividad general completa con el condensado de axiones (hipótesis H11).

Nada se transcribe a mano: se parte de las ecuaciones linealizadas desde la acción (pipeline.linear_equations,
las mismas que el paso 2 compara con la Ec. (3.3) del Paper I) y se convierten en una matriz numérica con
axgr.to_numeric / axgr._linear.

Variables:  y = (xi, DP, dth, dth')     con  DP = Delta p = delta p + xi p0'(r)  (presión lagrangiana)
            delta n = (DP - xi T) / (dp/dn),  T = p0'(r) escrito con la ecuación de fondo (TOV con axión)
Usar T de la TOV (en vez de dp/dn n0') evita derivadas terceras de la EOS: p0'' solo necesita Phi0'' y a0''.
Ligaduras algebraicas: delta lambda (G_tr) y delta Phi' (G_rr).  delta a = f_a dth.
"""
from functools import lru_cache

import numpy as np
import sympy as sp

from . import axgr as AG
from .model import r, w, G, n0, a0, dPhi, dlam, da1, dn1, xi1, p_of, eps0, Phi0, V
from .nonradial import fa
from .pipeline import linear_equations, background_equations, background_rules

d = lambda f, k=1: sp.diff(f, r, k)
DP = sp.Function('DeltaP')(r)
dth = sp.Function('dtheta')(r)
_nn = sp.Symbol('nn', positive=True)


@lru_cache
def reduced_equations():
    """Ecuaciones lineales en (xi, DP, dth) con el fondo ya reducido. Devuelve (eqs, pn, T)."""
    L = linear_equations()
    B = background_equations()
    pn = sp.diff(p_of(_nn), _nn).subs(_nn, n0)                      # dp/dn en n0
    # T = p0'(r) desde la ecuación de fondo del fluido (resuelta para n0')
    n0p = sp.solve(B['FLUID_r'], d(n0))[0]
    T = sp.expand(pn * n0p)
    a0pp = sp.solve(B['KG'], d(a0, 2))[0]
    dn_of = (DP - xi1 * T) / pn
    sub = {da1: fa * dth, dn1: dn_of}
    out = {}
    for k in ('Gtr', 'Grr', 'KG', 'FLUID_r', 'BAR'):
        e = L[k].subs(sub).doit()
        e = e.subs(d(a0, 3), d(a0pp)).subs(d(a0, 2), a0pp).doit()
        e = sp.expand(e.subs(d(a0, 2), a0pp))                        # la que reintroduce d(a0pp)
        # delta Phi sin derivar solo aparece en la conservación bariónica y su coeficiente es cero
        # módulo el fondo (el delta Phi de sqrt(-g) cancela el de u^t). Se comprueba y se elimina.
        c = e.coeff(dPhi)
        if c != 0:
            if sp.simplify(background_rules().reduce(c)) != 0:
                raise ValueError(f'{k}: delta Phi sin derivar con coeficiente no nulo')
            e = sp.expand(e - c * dPhi)
        out[k] = e
    return out, pn, T


def leftover_derivatives(e):
    """Derivativas que to_numeric no sabe sustituir (diagnóstico; debe estar vacío)."""
    e = AG.to_numeric(e)
    return sorted({str(x) for x in e.atoms(sp.Derivative)})


@lru_cache
def first_order_system():
    """Sistema de primer orden (xi', DP', dth'') en función de (xi, DP, dth, dth'), todo simbólico:
    delta lambda y su derivada desde G_tr, delta Phi' desde G_rr, xi' desde la conservación bariónica."""
    E, _, _ = reduced_equations()
    a0pp = sp.solve(background_equations()['KG'], d(a0, 2))[0]
    fix = lambda e: sp.expand(e.doit().subs(d(a0, 2), a0pp))
    dlam_s = sp.solve(E['Gtr'], dlam)[0]
    rep1 = lambda e: fix(e.subs(d(dlam), d(dlam_s)).subs(dlam, dlam_s))
    dPhip_s = sp.solve(rep1(E['Grr']), d(dPhi))[0]
    xip_s = sp.solve(rep1(E['BAR']), d(xi1))[0]
    rep2 = lambda e: sp.expand(rep1(e).subs(d(dPhi), dPhip_s).subs(d(xi1), xip_s))
    DPp_s = sp.solve(rep2(E['FLUID_r']), d(DP))[0]
    thpp_s = sp.solve(rep2(E['KG']), d(dth, 2))[0]
    return {'xi': xip_s, 'DP': DPp_s, 'thpp': thpp_s}


def _to_num(e):
    # V'(a0) = -M1/f_a en este modelo: V = E0 (1 - g) y n (m* - m_N) = -x E0 (1 - g)  (step15 lo comprueba)
    # Primero V'' (si no, subs de V' dentro de V'' = d(V')/da daría 0 al ser M1 un símbolo)
    e = e.doit().subs(sp.Derivative(V(a0), (a0, 2)), AG.V2 / AG.faS**2)
    e = e.subs(sp.Derivative(V(a0), a0), -AG.M1 / AG.faS)
    return sp.expand(AG.to_numeric(e))


@lru_cache
def interior_matrix():
    """A(r, w, fondo) con y' = A y dentro de la estrella, y = (xi, DP, dth, dth')."""
    F = first_order_system()
    S = AG._states([d(dth), xi1, DP, dth])
    ys = [S[xi1], S[DP], S[dth], S[d(dth)]]
    rows = []
    for key in ('xi', 'DP', None, 'thpp'):
        if key is None:
            rows.append([0, 0, 0, 1])
            continue
        e = _to_num(F[key]).subs(S)
        left = sp.expand(e - sum(e.coeff(y) * y for y in ys))
        if left != 0 or e.atoms(sp.Derivative):
            raise ValueError(f'{key}: sobra {left.free_symbols} / {e.atoms(sp.Derivative)}')
        rows.append([sp.expand(e.coeff(y)) for y in ys])
    return sp.lambdify((r, w, *AG.BG), sp.Matrix(rows), 'numpy', cse=True)


@lru_cache
def exterior_matrix():
    """Fuera de la estrella (sin fluido: xi = DP = 0, x = e0 = p = 0) pero con la cola del condensado:
    delta lambda y delta Phi' siguen acoplados al axión. y = (dth, dth')."""
    E, _, _ = reduced_equations()
    a0pp = sp.solve(background_equations()['KG'], d(a0, 2))[0]
    fix = lambda e: sp.expand(e.doit().subs({xi1: 0, DP: 0}).doit().subs(d(a0, 2), a0pp))
    Ex = {k: fix(v) for k, v in E.items()}
    dlam_s = sp.solve(Ex['Gtr'], dlam)[0]
    rep1 = lambda e: sp.expand(e.subs(d(dlam), d(dlam_s)).subs(dlam, dlam_s).doit().subs(d(a0, 2), a0pp))
    dPhip_s = sp.solve(rep1(Ex['Grr']), d(dPhi))[0]
    thpp = sp.solve(sp.expand(rep1(Ex['KG']).subs(d(dPhi), dPhip_s)), d(dth, 2))[0]
    zero = {AG.x: 0, AG.e0: 0, AG.pr: 0, AG.x1: 0}
    S = AG._states([d(dth), dth])
    e = sp.expand(_to_num(thpp).subs(zero).subs(S))
    ys = [S[dth], S[d(dth)]]
    left = sp.expand(e - sum(e.coeff(y) * y for y in ys))
    if left != 0 or e.atoms(sp.Derivative):
        raise ValueError(f'exterior: sobra {left.free_symbols} / {e.atoms(sp.Derivative)}')
    M = sp.Matrix([[0, 1], [sp.expand(e.coeff(ys[0])), sp.expand(e.coeff(ys[1]))]])
    return sp.lambdify((r, w, *AG.BG), M, 'numpy', cse=True)


def exterior_psi(omega, M, r_m, ma, fa_v, E0, n_decay=40.0, rtol=1e-11, sheet=1, dphi=0.0):
    """(psi, psi') del axión en r_m con condición saliente (o ligada), integrando la ecuación exterior
    derivada (exterior_matrix con theta_0 = 0: Schwarzschild + KG masivo, = Ec. 3.7 del Paper I) a lo largo
    del rayo complejo r = r_m + s e^{i phi}, con phi = pi/2 - arg(k): ahí la onda saliente decae como e^{-|k| s}.
    En la recta real (lo que hace axstar.exterior_L) la integración hacia dentro es inestable para modos muy
    amortiguados: la onda saliente crece con r y cualquier error inicial mezcla la entrante amplificada."""
    from scipy.integrate import solve_ivp
    from . import axstar as AX
    Ae = exterior_matrix()
    k = np.sqrt(complex(omega**2 - ma**2))
    if k.imag < 0 and abs(k.real) < 1e-14:
        k = -k
    k = sheet * k                      # sheet = -1: onda entrante (control negativo)
    phi = np.pi / 2 - np.angle(k) + dphi    # dphi != 0: prueba de independencia del contorno
    e = np.exp(1j * phi)
    S = n_decay / (abs(k) * np.cos(dphi))
    eta = M * (2 * omega**2 - ma**2) / k
    G2 = AX.g2(0.0)

    def bg(rr):
        f = 1 - 2 * M / rr
        P1 = M / (rr * (rr - 2 * M))
        P2 = -2 * (rr - M) * M / (rr * (rr - 2 * M))**2
        return (0.5 * np.log(f), P1, P2, -0.5 * np.log(f), -P1, 0, 0, 0.0, 0.0, 0, 0, 1.0,
                fa_v, 0.0, 0.0, E0 * G2, 0.0, -E0 * G2)

    rinf = r_m + S * e
    psi = np.exp(1j * (k * rinf + eta * np.log(rinf))) / rinf
    dpsi = psi * (1j * (k + eta / rinf) - 1 / rinf)
    sc = psi
    rhs = lambda s, y: e * (np.asarray(Ae(r_m + s * e, omega, *bg(r_m + s * e)), dtype=complex) @ y)
    sol = solve_ivp(rhs, (S, 0.0), [psi / sc, dpsi / sc], method='DOP853', rtol=rtol, atol=1e-60)
    return sol.y[0, -1], sol.y[1, -1]


# ---------------------------------------------------------------- integración numérica y modos
def _bgvals(st, rr):
    return AG.bgvals(st, rr)


def interior(st, omega, seed, r0=1e-3, r1=0.5, rtol=1e-10, r_eval=None):
    """Solución regular desde el centro hasta R. seed = 'fluid' (xi ~ r) o 'axion' (dth(0) = 1).
    Las condiciones en r0 salen de la propia fila de xi' (regularidad: xi'(r0) = 1 o 0)."""
    from scipy.integrate import solve_ivp
    A = interior_matrix()
    M = np.asarray(A(r0, omega, *_bgvals(st, r0)), dtype=complex)
    xi, th, dxi = (r0, 0.0, 1.0) if seed == 'fluid' else (0.0, 1.0, 0.0)
    y0 = np.array([xi, (dxi - M[0, 0] * xi - M[0, 2] * th) / M[0, 1], th, 0.0], dtype=complex)
    f = lambda rr, y: np.asarray(A(rr, omega, *_bgvals(st, rr)), dtype=complex) @ y
    s = solve_ivp(f, (r0, r1), y0, method='BDF', rtol=rtol, atol=1e-40)
    s = solve_ivp(f, (r1, st.R), s.y[:, -1], method='DOP853', rtol=rtol, atol=1e-40, t_eval=r_eval)
    return s


def outward(st, omega, y_R, rm, rtol=1e-10, r_eval=None):
    """(dth, dth') de R a r_m con la cola del condensado (exterior_matrix)."""
    from scipy.integrate import solve_ivp
    Ae = exterior_matrix()
    return solve_ivp(lambda rr, u: np.asarray(Ae(rr, omega, *_bgvals(st, rr)), dtype=complex) @ u,
                     (st.R, rm), y_R, method='DOP853', rtol=rtol, atol=1e-40, t_eval=r_eval)


def det(st, omega, rm=None, sheet=1, rtol=1e-10, r0=1e-3, n_decay=40.0):
    """Determinante de Delta p(R) = 0 y del empalme con la onda exterior (saliente o ligada) en r_m.
    sheet = -1 es la mutación «onda entrante» (control negativo). rtol, r0, n_decay: estudio de robustez."""
    rm = rm or AG.r_match(st)
    cols = []
    for seed in ('fluid', 'axion'):
        s = interior(st, omega, seed, r0=r0, rtol=rtol)
        cols.append((s.y[1, -1], outward(st, omega, s.y[2:, -1], rm, rtol=rtol).y[:, -1]))
    ps, dps = exterior_psi(omega, st.sol.y[0, -1], rm, st.ma, np.sqrt(st.fa2), st.E0, sheet=sheet,
                           n_decay=n_decay, rtol=min(rtol, 1e-11))
    L = dps / ps
    (a, ya), (c, yc) = cols
    return a * (yc[1] - L * yc[0]) - c * (ya[1] - L * ya[0])


def fluid_fraction(st, omega, rm=None):
    """Fracción de la energía (pesos cinéticos sencillos) que lleva el fluido en el modo de frecuencia omega."""
    rm = rm or AG.r_match(st)
    ri, re = np.linspace(0.5, st.R, 600), np.linspace(st.R, rm, 1500)
    sols = []
    for seed in ('fluid', 'axion'):
        s = interior(st, omega, seed, r_eval=ri)
        sols.append((s.y, outward(st, omega, s.y[2:, -1], rm, r_eval=re).y))
    (yf, ef), (ya, ea) = sols
    a, c = yf[1, -1], ya[1, -1]
    y, e = c * yf - a * ya, c * ef - a * ea
    Bi = np.array([_bgvals(st, x) for x in ri])
    Be = np.array([_bgvals(st, x) for x in re])
    eps = Bi[:, 9] + Bi[:, 5] * Bi[:, 13]
    Ef = np.trapezoid((eps + Bi[:, 10]) * np.exp(3 * Bi[:, 3] - Bi[:, 0]) * abs(y[0])**2 * ri**2, ri)
    Ea = (np.trapezoid(st.fa2 * np.exp(Bi[:, 3] - Bi[:, 0]) * abs(y[2])**2 * ri**2, ri)
          + np.trapezoid(st.fa2 * np.exp(Be[:, 3] - Be[:, 0]) * abs(e[0])**2 * re**2, re))
    return float(Ef / (Ef + Ea))


def muller(f, z, tol=1e-9, maxit=30):
    z0, z1, z2 = z, z + 0.01, z - 0.002j
    f0, f1, f2 = f(z0), f(z1), f(z2)
    for _ in range(maxit):
        h1, h2 = z1 - z0, z2 - z1
        d1, d2 = (f1 - f0) / h1, (f2 - f1) / h2
        a = (d2 - d1) / (h2 + h1)
        b = a * h2 + d2
        D = np.sqrt(b * b - 4 * f2 * a)
        den = b + D if abs(b + D) > abs(b - D) else b - D
        dz = -2 * f2 / den
        z0, z1, z2 = z1, z2, z2 + dz
        f0, f1, f2 = f1, f2, f(z2)
        if abs(dz) < tol * abs(z2):
            break
    return z2
