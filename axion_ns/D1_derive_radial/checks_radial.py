"""D1 checks for Paper I (old steps 1 and 2): derived equations vs the equations typed in the paper.

Logic and check names are those of axion_ns/step1_background.py and step2_perturbations.py. The derived
equations come from D0/D1 artifacts; the paper equations come from R's paper1_equations artifact.
"""
import sympy as sp

from common.modulo import check_modulo, is_zero
from common.notation import r, G, w, V, mstar, mN, Phi0, lam0, a0, n0, dPhi, dlam, da1, dn1, xi1

d = lambda f, k=1: sp.diff(f, r, k)

# derived equation that corresponds to each one of the paper, and the variable to eliminate
MATCH = {
    '(3.3a) KG': ('KG', d(da1, 2)),
    '(3.3b) rr': ('Grr', d(dPhi)),
    '(3.3c) tr': ('Gtr', dlam),
    '(3.3d) Euler': ('FLUID_r', xi1),
    '(3.3e) continuidad': ('BAR', dn1),
}


def background_checks(BG, rules, PB):
    """BG: static background equations (D0); PB: paper (2.6) equations keyed '(2.6a) p\\'' ..."""
    res = []
    for name, e in PB.items():
        res.append((f'Paper {name} == derivada de la acción', is_zero(rules.reduce(e)), ''))
    # limit without axion: V = 0, m* = m_N constant -> standard TOV
    x = sp.Symbol('x')
    std = {V: sp.Lambda(x, 0), mstar: sp.Lambda(x, mN)}
    tov = rules.reduce(PB["(2.6a) p'"]).subs(std).doit()
    res.append(('Límite sin axión reproduce TOV estándar (p\' = -(eps+p) Phi\')', is_zero(tov), ''))
    return res


def perturbation_checks(L, rules, PP):
    """L: linearized derived equations; PP: paper (3.3) equations."""
    R = rules
    res = []
    for name, (der, piv) in MATCH.items():
        ok, _ = check_modulo(PP[name], L[der], piv, R)
        res.append((f'Paper {name} == linealización de la acción', ok, ''))

    # negative controls: planted errors must be detected
    bad = {
        '(3.3b) con 8 pi en vez de 4 pi': ('(3.3b) rr', lambda e: e.subs(sp.pi, 2 * sp.pi)),
        '(3.3a) con signo cambiado en el término 2 delta lambda': ('(3.3a) KG', lambda e: e + 4 * dlam * sp.diff(V(a0), a0)),
        '(3.3e) sin el término delta lambda': ('(3.3e) continuidad', lambda e: e.subs(dlam, 0)),
    }
    for label, (name, mutate) in bad.items():
        der, piv = MATCH[name]
        ok, _ = check_modulo(mutate(PP[name]), L[der], piv, R)
        res.append((f'Control negativo: {label} -> se detecta', not ok, ''))

    # redundancy: the linearized tt component follows from tr (Bianchi identity)
    sol_dn = sp.solve(L['BAR'], dn1)[0]
    tt = L['Gtt'].subs({d(dn1): d(sol_dn)}).subs({dn1: sol_dn}).doit()
    sol_dlam = sp.solve(L['Gtr'], dlam)[0]
    tt = R.reduce(tt.subs({d(dlam): d(sol_dlam)}).subs({dlam: sol_dlam}).doit())
    res.append(('G_tt linealizada es redundante (Bianchi + conservación bariónica)', is_zero(tt), ''))

    # exterior: Schwarzschild + vacuum axion -> Eq. (3.7) of the paper
    Mc, ma = sp.symbols('M m_a', positive=True)
    x = sp.Symbol('x')
    f = 1 - 2 * G * Mc / r
    ext = PP['(3.3a) KG'].subs({V: sp.Lambda(x, ma**2 * x**2 / 2), mstar: sp.Lambda(x, 0)}).doit()
    ext = ext.subs({n0: 0, a0: 0, Phi0: sp.log(f) / 2, lam0: -sp.log(f) / 2}).doit()
    eq11 = f * d(da1, 2) + 2 / r * (1 - G * Mc / r) * d(da1) - (ma**2 - w**2 / f) * da1
    res.append(('Exterior: (3.3a) -> KG en Schwarzschild, Ec. (3.7)', sp.simplify(ext - eq11) == 0, ''))
    return res
