"""Paso 7 (trabajo nuevo): más allá de Cowling, con la métrica perturbada.

Es el primer "natural next step" que plantea el Paper II (Sec. IV). No hay ecuaciones publicadas con
las que comparar, así que la derivación se valida por sus límites y por consistencia interna:

  * la parte sin traza de las ecuaciones angulares da H2 = H0 (y el axión no la toca);
  * con la métrica congelada se recuperan exactamente las ecuaciones de Cowling del paso 5;
  * Bianchi: G_tt y la traza angular son redundantes, y la ligadura algebraica de H0 se conserva;
  * sin axión salen las ecuaciones de Lindblom y Detweiler (LD), contrastadas con dos transcripciones;
  * en el exterior: Zerilli para la métrica y la Ec. (13) del Paper II para el axión, desacoplados;
  * en el centro: las condiciones de regularidad de LD y sus correcciones axiónicas.

Después se escribe el sistema extendido (LD + axión) en forma compacta y se comprueba ecuación a
ecuación que se sigue de la acción.
"""
from functools import lru_cache

import sympy as sp

from . import fullgr as F
from . import nonradial as NR
from .compare import is_zero
from .model import r, Phi0, lam0, a0, n0, w, G, V, mstar, mN, eps0
from .pipeline import background_rules, background_equations, report
from .step5_nonradial import derived as derived5

d = F.d
H0, H1, H2, K, W, Vt, dn, dth = F.H0, F.H1, F.H2, F.K, F.W, F.Vt, F.dn, F.dth
l, fa = NR.l, NR.fa
PERT = [H0, H1, K, W, Vt, dn, dth, d(dth), F.Xf]

# ---------------------------------------------------------------- reducción "módulo el sistema"
_eS, _pS, _cS, _Vs, _tp, _dm, _dmm, _eL = sp.symbols('varepsilon p c_s2 V_a theta0p dm dmm eL', real=True)


def _phys(e):
    """Cambio de variables a cantidades físicas del fondo (biyectivo): ayuda a `simplify`."""
    e = e.subs(sp.Derivative(eps0(n0), (n0, 2)), _cS * (_eS + _pS) / n0**2)
    e = e.subs(sp.Derivative(eps0(n0), n0), (_pS + _eS) / n0 - (mstar(a0) - mN))
    e = e.subs(eps0(n0), _eS - n0 * (mstar(a0) - mN))
    e = e.subs(sp.Derivative(mstar(a0), (a0, 2)), _dmm / fa**2).subs(sp.Derivative(mstar(a0), a0), _dm / fa)
    e = e.subs(sp.Derivative(a0, r), fa * _tp).subs(V(a0), _Vs)
    return e.subs(sp.exp(2 * lam0), _eL).subs(sp.exp(lam0), sp.sqrt(_eL))


def vanishes(e):
    """Cero coeficiente a coeficiente en las perturbaciones, tras el cambio de variables."""
    e = sp.expand(e)
    for v in PERT:
        c = e.coeff(v)
        if c != 0 and sp.simplify(_phys(c)) != 0 and not is_zero(c):
            return False
        e = e - c * v
    return sp.simplify(_phys(sp.expand(e))) == 0


@lru_cache
def system():
    """Ecuaciones linealizadas, reglas para las derivadas y las dos ligaduras algebraicas."""
    L = F.linear_system()
    traceless = sp.factor(sp.simplify(L['Etl']))
    L = {k: sp.expand(v.subs(H2, H0).doit()) for k, v in L.items()}
    R = background_rules()
    piv = [d(K), d(H1), d(H0), d(W), d(dn), d(dth, 2)]
    sol = sp.solve([L[k] for k in ('Etr', 'Etth', 'Erth', 'BAR', 'FLUID_r', 'KG')], piv, dict=True)[0]
    first = {k: R.reduce(sp.expand(v)) for k, v in sol.items()}

    def sub1(e):
        for _ in range(3):
            e = e.subs(first)
        return e

    second = {d(f, 2): R.reduce(sub1(d(first[d(f)]))) for f in (K, H1, H0, W, dn)}
    second[d(dth, 3)] = R.reduce(sub1(d(first[d(dth, 2)])))

    def red(e):
        return sp.expand(R.reduce(sub1(e.subs(second))))

    AV, AH = red(L['FLUID_th']), red(L['Err'])
    solH0 = sp.solve(AH, H0)[0]
    solV = sp.solve(AV.subs(H0, solH0), Vt)[0]

    def elim(e):
        e = red(e.subs(F.Xf, F.X_def).doit())
        return e.subs(Vt, solV).subs(H0, solH0)

    return dict(L=L, traceless=traceless, first=first, red=red, elim=elim, AH=AH, R=R)


def noaxion(e):
    e = e.subs(dth, 0).subs(a0, 0).doit()
    e = e.replace(lambda x: isinstance(x, sp.Subs), lambda x: 0)
    return e.subs({mstar(0): mN, V(0): 0})


def implied(e, axion=True):
    """¿Se sigue la ecuación `e` del sistema derivado (con o sin axión)?"""
    res = system()['elim'](e)
    return vanishes(res) if axion else is_zero(noaxion(res))


# ---------------------------------------------------------------- comprobaciones
def structure_checks():
    S = system()
    L = S['L']
    res = [('Parte sin traza de las ecuaciones angulares: H2 = H0, sin contribución del axión',
            sp.simplify(S['traceless'] + (H0 - H2) / 2) == 0, f'{S["traceless"]}')]
    _, L5 = derived5()
    zero = {H0: 0, H1: 0, K: 0}
    ok = all(sp.simplify(sp.powdenest(L[k].subs(zero).doit() - L5[name], force=True)) == 0 for k, name in
             [('FLUID_r', '(8a) momento radial'), ('BAR', '(8b) bariones'),
              ('FLUID_th', '(8c) momento angular'), ('KG', '(8d) KG interior')])
    res.append(('Métrica congelada (H0 = H1 = K = 0): las cuatro ecuaciones de materia son las de Cowling (paso 5)', ok, ''))
    for k, label in [('Ett', 'G_tt'), ('Eang', 'la traza angular')]:
        res.append((f'Bianchi: {label} linealizada es redundante (con axión)', vanishes(S['elim'](L[k])), ''))
    res.append(('La ligadura algebraica de H0 se conserva al integrar (su derivada se anula con el sistema)',
                vanishes(S['elim'](d(S['AH']))), ''))
    return res


def ld_checks():
    res = []
    ok = {k: implied(e, axion=False) for k, e in F.LD_ZL.items()}
    good = [k for k in ok if ok[k]]
    res.append(('Sin axión: H1\', K\', W\', H0 y V de Lindblom-Detweiler (Zhao y Lattimer, arXiv:2204.03037)',
                good == ["H1'", "K'", "W'", 'H0', 'V'], f'coinciden: {good}'))
    fix = implied(F.LD_ZL_FIX["X' con -2w^2"], axion=False)
    res.append(('X\' de Zhao y Lattimer: no coincide tal cual; con -2 w^2 en el término de W, sí (errata de signo)',
                (not ok["X'"]) and fix, ''))
    bu = {k: implied(e, axion=False) for k, e in F.LD_BU.items()}
    bu_fix = {k: implied(e, axion=False) for k, (_, e) in F.LD_BU_FIX.items()}
    bad = [k for k in bu if not bu[k]]
    res.append(('Burgio et al. (arXiv:1106.2736): W\' coincide; H1\', K\', X\', H0 y V fallan tal cual y coinciden '
                'con una corrección mínima cada una', bad == list(F.LD_BU_FIX) and all(bu_fix.values()),
                '; '.join(f'{k}: {F.LD_BU_FIX[k][0]}' for k in bad)))
    ku = [k for k, e in F.LD_KU.items() if implied(e, axion=False)]
    res.append(('Kunjipurayil et al. (arXiv:2205.02081, coautor T. Zhao): las seis coinciden, X\' con -2 w^2',
                ku == list(F.LD_KU), f'coinciden: {ku}'))
    zh = {k: implied(e, axion=False) for k, e in F.LD_ZH.items()}
    res.append(('Zhang et al. (arXiv:2601.21911): H1\', K\', W\', X\' y H0 coinciden, escritas ya con las cinco '
                'correcciones de Burgio et al.; en V pone w^3 donde va w^2',
                [k for k in zh if not zh[k]] == ['V'] and bu_fix['V'], ''))
    wrong = F.LD_ZL["H1'"] + 8 * sp.pi * G * F.ep * F.eL * Vt   # 16 pi -> 8 pi
    res.append(('Control negativo: 16 pi -> 8 pi en H1\' de LD se detecta', not implied(wrong, axion=False), ''))
    return res


def exterior_checks():
    L = system()['L']
    M = sp.Symbol('M', positive=True)
    f = 1 - 2 * G * M / r
    first_der = lambda e: isinstance(e, sp.Subs) and isinstance(e.expr, sp.Derivative) and e.expr.derivative_count == 1

    def ext(e, axion=False):
        if not axion:
            e = e.subs(dth, 0)
        e = e.subs({W: 0, Vt: 0, dn: 0}).subs(a0, 0).subs(n0, 0).doit()
        e = e.replace(first_der, lambda x: 0)          # V'(0) = m*'(0) = 0: theta = 0 es un mínimo
        if not axion:
            e = e.replace(lambda x: isinstance(x, sp.Subs), lambda x: 0)
        e = e.subs({mstar(0): mN, V(0): 0, eps0(0): 0})
        e = e.subs({Phi0: sp.log(f) / 2, lam0: -sp.log(f) / 2}).doit()
        return sp.expand(sp.simplify(e))

    tr, tth, rth, rr = [ext(L[k]) for k in ('Etr', 'Etth', 'Erth', 'Err')]
    sol = sp.solve([tr, tth, rth], [d(K), d(H1), d(H0)], dict=True)[0]
    solH0 = sp.solve(sp.simplify(rr.subs(sol)), H0)[0]
    D = lambda e: sp.simplify(d(e).subs(sol).subs(H0, solH0))
    nn = F.nn
    VZ = f * (2 * nn**2 * (nn + 1) * r**3 + 6 * nn**2 * G * M * r**2 + 18 * nn * G**2 * M**2 * r + 18 * G**3 * M**3) \
        / (r**3 * (nn * r + 3 * G * M)**2)
    Z = r**(l + 2) / (nn * r + 3 * G * M) * (K - f * H1)
    zer = f * D(f * D(Z)) + (w**2 - VZ) * Z
    res = [('Exterior sin axión: Z = r^{l+2}(K - e^psi H1)/(nr+3M) cumple la ecuación de Zerilli',
            sp.simplify(zer) == 0, '')]
    kg = ext(L['KG'], axion=True)
    u = sp.Function('u')(r)
    m2 = sp.Symbol('m_a2')                              # V''(0) = m_a^2 f_a^2 ... sólo queda esa Subs
    kgu = kg.subs(dth, (NR.Rs / r)**l * u).doit() * (r / NR.Rs)**l / fa
    kgu = kgu.replace(lambda e: isinstance(e, sp.Subs), lambda e: m2)
    eq13 = d(u, 2) + (2 / r + 2 * G * M / (r**2 * f)) * d(u) - (m2 + l * (l + 1) / r**2 - w**2 / f) * u / f
    # (la ecuación derivada es -(1 - 2GM/r) veces la (13))
    res.append(('Exterior con axión (theta_0 = 0): delta theta no ve la métrica y cumple la Ec. (13) del Paper II',
                not kg.has(H0, H1, K) and sp.simplify(kgu + f * eq13) == 0,
                ''))
    return res


def origin_checks():
    """Desarrollo en r = 0 con un fondo regular que cumple las ecuaciones de fondo."""
    L = system()['L']
    Pc, P2, L2, nc, N2, ac, A2 = sp.symbols('Phi_c Phi_2 lambda_2 n_c n_2 a_c a_2')
    e0, e1, e2, e3, m0, m1, m2, v0, v1, v2 = sp.symbols('e0:4 m0:3 v0:3')
    x = sp.Symbol('x')
    funcs = {eps0: sp.Lambda(x, e0 + e1 * (x - nc) + e2 * (x - nc)**2 / 2 + e3 * (x - nc)**3 / 6),
             mstar: sp.Lambda(x, m0 + m1 * (x - ac) + m2 * (x - ac)**2 / 2),
             V: sp.Lambda(x, v0 + v1 * (x - ac) + v2 * (x - ac)**2 / 2)}
    bg = {Phi0: Pc + P2 * r**2 / 2, lam0: L2 * r**2 / 2, n0: nc + N2 * r**2 / 2, a0: ac + A2 * r**2 / 2}
    concrete = lambda e: e.subs(funcs).doit().subs(bg).doit()
    B = background_equations()
    beqs = []
    for k in ('Gtt', 'Grr', 'FLUID_r', 'KG'):
        s = sp.series(concrete(B[k]) * r, r, 0, 3).removeO()
        beqs += [q for q in (sp.expand(s.coeff(r, i)) for i in range(3)) if q != 0]
    bs = sp.solve(beqs, [L2, P2, N2, A2], dict=True)[0]
    c = sp.symbols('h0 h2 j0 j2 k0 k2 W0 W2 U0 U2 N0 N2p T0 T2')
    h0, h2, j0, j2, k0, k2, W0, W2, U0, U2, N0, N2p, T0, T2 = c
    pert = {H0: h0 + h2 * r**2, H1: j0 + j2 * r**2, K: k0 + k2 * r**2, W: W0 + W2 * r**2, Vt: U0 + U2 * r**2,
            dn: N0 + N2p * r**2, dth: T0 + T2 * r**2}
    eqs = []
    for k in ('Etr', 'Etth', 'Erth', 'Err', 'BAR', 'FLUID_r', 'FLUID_th', 'KG'):
        e = concrete(L[k]).subs(pert).doit().subs(bs)
        s = sp.series(e * r**2, r, 0, 3).removeO()
        eqs += [(k, q) for q in (sp.simplify(s.coeff(r, i)) for i in range(3)) if q != 0]
    sol = sp.solve([q for _, q in eqs], [h0, j0, U0, N0, h2, j2, k2, W2, U2, N2p, T2], dict=True)
    res = []
    if len(sol) != 1:
        return [('Regularidad en el centro: solución única', False, f'{len(sol)} soluciones')]
    sol = sol[0]
    epc = nc * e1 + nc * (m0 - mN)                     # (eps + p) del fluido en el centro
    Ec, pc = e0 + nc * (m0 - mN), nc * e1 - e0
    res.append(('Centro: H0(0) = K(0), H1(0) = [2l K(0) + 16 pi G (eps+p) W(0)]/[l(l+1)] y V(0) = -W(0)/l, '
                'como en LD, sin corrección axiónica',
                sp.simplify(sol[h0] - k0) == 0 and sp.simplify(sol[U0] + W0 / l) == 0
                and sp.simplify(sol[j0] - (2 * l * k0 + 16 * sp.pi * G * epc * W0) / (l * (l + 1))) == 0, ''))
    X0 = -sp.exp(Pc) * (nc * e2 * sol[N0] + W0 * nc * e2 * bs[N2])
    X0_LD = epc * sp.exp(Pc) * ((sp.Rational(4, 3) * sp.pi * G * (Ec + 3 * pc) - w**2 * sp.exp(-2 * Pc) / l) * W0 + k0 / 2)
    dmc, thpp = fa * m1, bs[A2] / fa                   # dm*/dtheta y theta_0''(0)
    X0_ax = X0_LD + sp.exp(Pc) * (nc * dmc * (T0 + thpp * W0) - sp.Rational(8, 3) * sp.pi * G * epc * v0 * W0)
    res.append(('Centro: X(0) = X_LD(0) + e^{Phi_c}[n dm*/dtheta Delta theta(0) - (8 pi G/3)(eps+p) V(theta_c) W(0)]',
                sp.simplify(X0 - X0_ax) == 0, 'Delta theta(0) = delta theta(0) + theta_0\'\'(0) W(0)'))
    kg = [q for k, q in eqs if k == 'KG'][0]
    T2s = sp.solve(kg.subs({h0: k0, j0: sol[j0], U0: -W0 / l}), T2)[0]
    c16 = (T0 * (-w**2 * sp.exp(-2 * Pc) - l * bs[P2] + l * (l + 2) * bs[L2] + v2 + nc * m2) + N0 * m1 / fa) / (2 * l + 3)
    res.append(('Centro: delta theta\'\'(0) = Ec. (16c) del Paper II - 3 K(0) theta_0\'\'(0)/(2l+3)',
                sp.simplify(2 * T2s - c16 + 3 * k0 * thpp / (2 * l + 3)) == 0, ''))
    return res


def extended_checks():
    res = []
    for k, label in [("H1'", "H1'"), ("K'", "K'"), ("W'", "W' (sin cambios)"), ("X'", "X'"),
                     ('H0', 'relación algebraica de H0'), ('V', 'relación algebraica de V'),
                     ('KG', 'Klein-Gordon: Ec. (8d) + términos métricos')]:
        res.append((f'Sistema extendido, {label}: se sigue de la acción', implied(F.EXT[k]), ''))
    ok = all(is_zero(noaxion(F.NEW[k].subs(F.Xf, F.X_def).doit())) for k in F.NEW)
    res.append(('Todos los términos nuevos se anulan sin axión', ok, ''))
    wrong = F.EXT["K'"] + 8 * sp.pi * G * r * fa**2 * F.thp * dth     # quitar el término de flujo
    res.append(('Control negativo: quitar el flujo 8 pi G r f_a^2 theta_0\' delta theta de K\' se detecta',
                not implied(wrong), ''))
    return res


def run():
    out = []
    for title, fn in [('Estructura y límites', structure_checks), ('Lindblom-Detweiler', ld_checks),
                      ('Exterior', exterior_checks), ('Centro', origin_checks), ('Sistema extendido', extended_checks)]:
        out += report(f'Paso 7: {title}', fn())
    return out


if __name__ == '__main__':
    run()
