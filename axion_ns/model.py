"""Modelo físico: estrella de neutrones + axión ligero acoplado a los nucleones.

Referencia: A. Gómez-Bañón, P. Pnigouras, J. A. Pons, arXiv:2510.16245 (2025).

Punto de partida (signatura (+,-,-,-), c = hbar = 1, G explícita):

    S = Int d^4x sqrt(-g) [ -R/(16 pi G) + 1/2 g^{mu nu} d_mu a d_nu a - V(a) - eps(n, a) ]

con el fluido perfecto descrito por su densidad de energía eps(n, a). El acoplo axión-nucleón
entra a través de la masa efectiva m*(a) (con n_s = n_b, como en el paper):

    eps(n, a) = eps0(n) + n [m*(a) - m_N]      =>   p = n d eps/dn - eps = n eps0' - eps0 = p0(n)

Todo lo que se deriva aquí sale de esta acción: Einstein (geometría + T_{mu nu} por variación
métrica), Klein-Gordon (Euler-Lagrange) y fluido (conservación de T total + bariones).
Los campos son funciones generales de (t, r); el fondo estático y las perturbaciones se
obtienen luego por sustitución.
"""
import sympy as sp

from .geometry import einstein, div_mixed

# ---------------------------------------------------------------- símbolos y campos
t, r, th, ph = sp.symbols('t r theta phi', real=True)
X = [t, r, th, ph]
G = sp.Symbol('G', positive=True)
mN = sp.Symbol('m_N', positive=True)

Phi = sp.Function('Phi')(t, r)
lam = sp.Function('lambda')(t, r)
a = sp.Function('a')(t, r)
n = sp.Function('n')(t, r)
xi = sp.Function('xi')(t, r)  # desplazamiento radial del fluido

V = sp.Function('V')              # potencial del axión en vacío
mstar = sp.Function('mstar')      # masa efectiva del nucleón m*(a)
eps0 = sp.Function('epsilon0')    # densidad de energía sin axión eps0(n)


def eps_of(nn, aa):
    return eps0(nn) + nn * (mstar(aa) - mN)


def p_of(nn):
    # p = n d eps/dn - eps  (la dependencia en a se cancela exactamente)
    x = sp.Dummy('x')
    e = eps0(x)
    return (x * sp.diff(e, x) - e).subs(x, nn)


# ---------------------------------------------------------------- geometría
g = sp.diag(sp.exp(2 * Phi), -sp.exp(2 * lam), -r**2, -r**2 * sp.sin(th)**2)
sqrtg = sp.exp(Phi + lam) * r**2 * sp.sin(th)


def derive_field_equations():
    """Devuelve un dict con las ecuaciones de campo completas (= 0), dependientes de (t, r)."""
    Gmn, Gam, ginv = einstein(g, X)
    da = [sp.diff(a, x) for x in X]

    # --- Klein-Gordon por Euler-Lagrange de la densidad sqrt(-g) L_a (n fijo)
    L_a = sp.Rational(1, 2) * sum(ginv[i, j] * da[i] * da[j] for i in range(4) for j in range(4)) \
        - V(a) - eps_of(n, a)
    KG_el = sp.euler_equations(sqrtg * L_a, a, [t, r])[0].lhs
    # normalización: KG = box a + V'(a) + n m*'(a)
    KG = sp.simplify(-KG_el / sqrtg)

    # --- T_{mu nu} del axión por variación respecto a la métrica inversa:
    #     T_{mu nu} = 2 dL/dg^{mu nu} - g_{mu nu} L
    h = sp.Matrix(4, 4, lambda i, j: sp.Symbol(f'h{i}{j}'))
    L_h = sp.Rational(1, 2) * sum(h[i, j] * da[i] * da[j] for i in range(4) for j in range(4)) - V(a)
    L_val = L_h.subs({h[i, j]: ginv[i, j] for i in range(4) for j in range(4)})
    T_ax = sp.Matrix(4, 4, lambda i, j: 2 * sp.diff(L_h, h[i, j]) - g[i, j] * L_val)

    # --- fluido perfecto: u^r = u^t d_t xi, normalizada
    ut = 1 / sp.sqrt(sp.exp(2 * Phi) - sp.exp(2 * lam) * sp.diff(xi, t)**2)
    u_up = sp.Matrix([ut, ut * sp.diff(xi, t), 0, 0])
    u_dn = g * u_up
    e, p = eps_of(n, a), p_of(n)
    T_fl = sp.Matrix(4, 4, lambda i, j: (e + p) * u_dn[i] * u_dn[j] - p * g[i, j])

    T = T_fl + T_ax
    EIN = {(i, j): Gmn[i, j] - 8 * sp.pi * G * T[i, j] for (i, j) in [(0, 0), (0, 1), (1, 1)]}

    # --- conservación: nabla_mu T^mu_nu(total) = 0, eliminando KG
    mixed = lambda Tdn: ginv * Tdn  # T^mu_nu = g^{mu a} T_{a nu}
    divT_fl = div_mixed(mixed(T_fl), Gam, sqrtg, X, 1)
    divT_ax = div_mixed(mixed(T_ax), Gam, sqrtg, X, 1)
    # identidad: nabla T_ax^mu_r = (box a + V') d_r a  =>  con KG: = -n m*' d_r a
    ident_ax = sp.simplify(divT_ax - (KG - n * sp.diff(mstar(a), a)) * sp.diff(a, r))
    FLUID_r = divT_fl + divT_ax - KG * sp.diff(a, r)

    # --- conservación del número bariónico: d_mu (sqrt(-g) n u^mu) = 0
    BAR = sum(sp.diff(sqrtg * n * u_up[m], X[m]) for m in range(4)) / (r**2 * sp.sin(th))

    return dict(Gtt=EIN[(0, 0)], Gtr=EIN[(0, 1)], Grr=EIN[(1, 1)], KG=KG,
                FLUID_r=FLUID_r, BAR=BAR, ident_ax=ident_ax)


# ---------------------------------------------------------------- fondo y perturbaciones
Phi0, lam0, a0, n0 = [sp.Function(s)(r) for s in ('Phi0', 'lambda0', 'a0', 'n0')]
dPhi, dlam, da1, dn1, xi1 = [sp.Function(s)(r) for s in ('deltaPhi', 'deltalambda', 'deltaa', 'deltan', 'xi1')]
w = sp.Symbol('omega')
eta = sp.Symbol('eta')


def background(expr):
    """Configuración estática (orden eta^0)."""
    return sp.simplify(expr.subs({Phi: Phi0, lam: lam0, a: a0, n: n0, xi: 0}).doit())


def linearize(expr):
    """Parte lineal en las perturbaciones, con dependencia e^{-i w t} factorizada."""
    E = sp.exp(-sp.I * w * t)
    sub = {Phi: Phi0 + eta * dPhi * E, lam: lam0 + eta * dlam * E, a: a0 + eta * da1 * E,
           n: n0 + eta * dn1 * E, xi: eta * xi1 * E}
    lin = sp.diff(expr.subs(sub).doit(), eta).subs(eta, 0)
    return sp.expand(sp.powsimp(sp.expand(lin / E)))
