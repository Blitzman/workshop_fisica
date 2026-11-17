"""Notation contract shared by every symbolic module (D0, D1, D2, R).

This file is the ONLY place where the names of the symbols and undefined functions used in the
equation artifacts are defined. An artifact stores expressions in terms of exactly these objects
(symbol name + assumptions, function name), so any module that imports this file reads the same
expression back that the producer wrote. It holds definitions only: no derivation, no
paper equations.

Provenance: the contents are taken verbatim from axion_ns/model.py (symbols and fields),
axion_ns/nonradial.py (Paper II perturbation variables and shorthand) and axion_ns/fullgr.py
(metric-perturbation variables, Lindblom-Detweiler shorthand).

Conventions are those of CONVENTIONS.md: signature (+,-,-,-), c = hbar = 1, G explicit,
e^{-i w t} Y_lm, factor (r/R)^l extracted from the perturbations.
"""
import sympy as sp

# ------------------------------------------------------------------ coordinates, constants
t, r, th, ph = sp.symbols('t r theta phi', real=True)
X = [t, r, th, ph]
G = sp.Symbol('G', positive=True)
mN = sp.Symbol('m_N', positive=True)
w = sp.Symbol('omega')
eta = sp.Symbol('eta')

# ------------------------------------------------------------------ undefined functions of the model
V = sp.Function('V')                # axion potential in vacuum
mstar = sp.Function('mstar')        # effective nucleon mass m*(a)
eps0 = sp.Function('epsilon0')      # energy density without axion, eps0(n)


def eps_of(nn, aa):
    """eps(n, a) = eps0(n) + n [m*(a) - m_N]."""
    return eps0(nn) + nn * (mstar(aa) - mN)


def p_of(nn):
    """p = n d eps/dn - eps  (the dependence on a cancels exactly)."""
    x = sp.Dummy('x')
    e = eps0(x)
    return (x * sp.diff(e, x) - e).subs(x, nn)


# ------------------------------------------------------------------ general fields of (t, r)  [D0]
Phi = sp.Function('Phi')(t, r)
lam = sp.Function('lambda')(t, r)
a = sp.Function('a')(t, r)
n = sp.Function('n')(t, r)
xi = sp.Function('xi')(t, r)        # radial displacement of the fluid

# ------------------------------------------------------------------ static background and radial perturbations
Phi0, lam0, a0, n0 = [sp.Function(s)(r) for s in ('Phi0', 'lambda0', 'a0', 'n0')]
dPhi, dlam, da1, dn1, xi1 = [sp.Function(s)(r) for s in ('deltaPhi', 'deltalambda', 'deltaa', 'deltan', 'xi1')]

# ------------------------------------------------------------------ nonradial perturbations (Paper II convention)
l = sp.Symbol('l', positive=True, integer=True)
Rs = sp.Symbol('R', positive=True)
fa = sp.Symbol('f_a', positive=True)

# full fields of (t, r, theta)
A = sp.Function('a')(t, r, th)
N = sp.Function('n')(t, r, th)
XR = sp.Function('xi_r')(t, r, th)
XT = sp.Function('xi_th')(t, r, th)

# perturbations in the variables of the paper
dth, dn, W, Vt = [sp.Function(s)(r) for s in ('deltatheta', 'deltan', 'W', 'calV')]
Y = sp.Function('Y')(th)

# ------------------------------------------------------------------ metric perturbation (Lindblom-Detweiler 1983 convention)
H0, H1, H2, K = [sp.Function(s)(r) for s in ('H0', 'H1', 'H2', 'K')]
Xf = sp.Function('X')(r)            # Lagrangian pressure variable X = -e^Phi Delta p


# ------------------------------------------------------------------ shorthand (functions of the background)
def d(f, k=1):
    """k-th radial derivative."""
    return sp.diff(f, r, k)


eps, p = eps_of(n0, a0), p_of(n0)
_x = sp.Dummy('x')
dpdn = sp.diff(p_of(_x), _x).subs(_x, n0)
cs2 = dpdn * n0 / (eps + p)
dm = fa * sp.diff(mstar(a0), a0)            # dm*/dtheta
dmm = fa**2 * sp.diff(mstar(a0), a0, 2)     # d^2 m*/dtheta^2
thp = d(a0) / fa                            # theta_0'
meff2 = (fa**2 * sp.diff(V(a0), a0, 2) + n0 * dmm) / fa**2
dp = dpdn * dn                              # Eulerian perturbation of p
de = (eps + p) / n0 * dn + n0 * dm * dth    # first law with eps(n, theta)

# Lindblom-Detweiler shorthand (fullgr.py)
ep = eps + p
eL, eh = sp.exp(2 * lam0), sp.exp(lam0)
psi = 2 * Phi0
nn = (l - 1) * (l + 2) / 2
pp = dpdn * d(n0)
X_def = -sp.exp(Phi0) * (dp + sp.exp(-lam0) * W / r * pp)
GM = r * (1 - sp.exp(-2 * lam0)) / 2
b_TOV = GM / r
Q_TOV = b_TOV + 4 * sp.pi * G * r**2 * p            # Q of LD read from the fluid
Q_met = r * sp.exp(-2 * lam0) * d(Phi0)             # the same Q read from the metric
Pp = d(Phi0)
Dtheta = dth + sp.exp(-lam0) * thp * W / r          # Lagrangian perturbation of theta
