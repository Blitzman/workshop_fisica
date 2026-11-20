"""Ecuaciones de referencia para oscilaciones no barotrópicas (modos g y de interfaz axiónica).

Transcribe de forma independiente:
- Las ecuaciones relativistas no barotrópicas en la aproximación de Cowling de
  Reisenegger & Goldreich (1992, ApJ 395, 240) y Finn (1987, MNRAS 227, 265),
  extendidas con el acoplo axiónico de Gómez-Bañón, Pnigouras & Pons (Paper II, arXiv:2610.37308).
- El discriminante relativista de Schwarzschild con axión A(r) y la frecuencia de
  Brunt-Väisälä relativista N^2(r).
- Las relaciones eulerianas para delta eps y delta n con composición congelada (Delta Y_p = 0).
- La relación de dispersión WKB local para ondas internas de gravedad (límite de Boussinesq).

NOTA (Regla 1 de AGENTS.md): `axion_ns/gmodes.py` jamás importa este módulo.
"""
import sympy as sp

from .model import r, w, V, mstar, eps_of, p_of, Phi0, lam0, a0, n0
from .nonradial import l, fa, dth, W, Vt

d = lambda f, k=1: sp.diff(f, r, k)

# Variables termodinámicas de fondo y perturbación euleriana de presión
eps, p = eps_of(n0, a0), p_of(n0).doit()
dpdn_eq = sp.diff(p, n0)

# Velocidad del sonido de equilibrio c_e^2 y adiabática (composición congelada) c_s^2
ce2 = dpdn_eq * n0 / (eps + p)
cs2 = sp.Symbol('c_s2', positive=True)
dp = sp.Function('deltap')(r)
kr = sp.Symbol('k_r', real=True)

# Derivadas del acoplo axiónico respecto a theta = a / f_a y respecto a a
ma = sp.diff(mstar(a0), a0)
maa = sp.diff(mstar(a0), a0, 2)
dm = fa * ma
dmm = fa**2 * maa
thp = d(a0) / fa
meff2 = (fa**2 * sp.diff(V(a0), a0, 2) + n0 * dmm) / fa**2

# Desplazamiento radial xi^r = e^{-lambda_0} W / r
xi_r = sp.exp(-lam0) * W / r

# Gravedad efectiva con fuerza escalar del axión (del equilibrio hidrostático TOV+axión)
g_eff_p = -d(p) / (eps + p)
g_eff = d(Phi0) + n0 / (eps + p) * ma * d(a0)

# Discriminante relativista de Schwarzschild con axión A(r) (dos formas equivalentes)
A_def = (d(eps) - d(p) / cs2 - n0 * ma * d(a0)) / (eps + p)
A_ce = (1 / ce2 - 1 / cs2) * d(p) / (eps + p)

# Frecuencia de Brunt-Väisälä relativista al cuadrado N^2(r)
N2_def = -g_eff * A_def * sp.exp(2 * (Phi0 - lam0))
N2_ce = g_eff**2 * (1 / ce2 - 1 / cs2) * sp.exp(2 * (Phi0 - lam0))

# Perturbaciones eulerianas con composición congelada: Delta p = c_s^2 (Delta eps - n_b dm*/da Delta a)
de_euler = dp / cs2 - (eps + p) * A_def * xi_r + n0 * dm * dth
dn_euler = n0 / (eps + p) * dp / cs2 - n0 * A_def * xi_r

# Sistema no barotrópico de 4 ecuaciones en Cowling (Reisenegger & Goldreich 1992 + Paper II)
REF_EQS = {
    '(1) momento radial': d(dp) + (d(Phi0) + l / r) * dp - (
        (w**2 - N2_def) * sp.exp(2 * (lam0 - Phi0)) * sp.exp(-lam0) / r * (eps + p) * W
        - (dp / cs2 + n0 * dm * dth) * d(Phi0)
        - n0 * dm * (d(dth) + l / r * dth)
        - thp * (n0 / (eps + p) * dp / cs2 * dm + n0 * dmm * dth)
    ),
    '(2) bariones': d((eps + p) * W) + (eps + p) * W * (d(Phi0) + (l + 1) / r - A_def)
    + r * sp.exp(lam0) * (dp / cs2 + (eps + p) * Vt * l * (l + 1) / r**2),
    '(3) momento angular': (eps + p) * Vt + sp.exp(2 * Phi0) / w**2 * (dp + n0 * dm * dth),
    '(4) KG interior': sp.exp(-2 * lam0) * d(dth, 2)
    + sp.exp(-2 * lam0) * (d(Phi0) - d(lam0) + 2 * (l + 1) / r) * d(dth)
    + (w**2 * sp.exp(-2 * Phi0) - meff2 + l / r * sp.exp(-2 * lam0) * (d(Phi0) - d(lam0))
       + l * (l + 1) / r**2 * (sp.exp(-2 * lam0) - 1)) * dth
    - n0 / fa**2 * dm * (dp / (cs2 * (eps + p)) - A_def * xi_r),
}

# Relación de dispersión WKB local para ondas internas de gravedad (Boussinesq, omega^2 << L_l^2)
WKB_OMEGA2 = l * (l + 1) * N2_def / (kr**2 * r**2 * sp.exp(2 * lam0) + l * (l + 1))
