"""Paso 1: acción -> ecuaciones de equilibrio (TOV + Klein-Gordon) y comparación con Ec. (2.6)."""
import sympy as sp

from . import paper
from .compare import is_zero
from .model import eps0, mstar, mN, V
from .pipeline import field_equations, background_equations, background_rules, report


def run():
    res = [('Identidad nabla_mu T^mu_r(axion) = (box a + V\') d_r a', field_equations()['ident_ax'] == 0, '')]
    bg, R = background_equations(), background_rules()
    res.append(('Ec. G_tr del fondo estático se anula idénticamente', bg['Gtr'] == 0, ''))
    for name, e in paper.BACKGROUND.items():
        res.append((f'Paper {name} == derivada de la acción', is_zero(R.reduce(e)), ''))

    # límite sin axión: V = 0, m* = m_N constante -> TOV estándar
    x = sp.Symbol('x')
    std = {V: sp.Lambda(x, 0), mstar: sp.Lambda(x, mN)}
    tov = R.reduce(paper.BACKGROUND["(2.6a) p'"]).subs(std).doit()
    res.append(('Límite sin axión reproduce TOV estándar (p\' = -(eps+p) Phi\')', is_zero(tov), ''))

    # observación: p(n, a) no depende de a (exacto, no aproximado, con la Ec. (2.7) del paper)
    from .model import eps_of, n0, a0
    nn = sp.Symbol('n')
    p_exact = (nn * sp.diff(eps_of(nn, a0), nn) - eps_of(nn, a0))
    res.append(('p = n d(eps)/dn - eps no depende de a (la Ec. (2.8) del paper es exacta)',
                sp.diff(p_exact, a0) == 0, 'el paper la presenta como aproximación (≈)'))
    return report('Paso 1: fondo (TOV + KG)', res)


if __name__ == '__main__':
    run()
