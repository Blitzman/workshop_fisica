"""Ecuaciones de referencia transcritas de la literatura para rotación lenta y modos axiales (r-modes).

Referencias independientes:
  1. J. B. Hartle (1967), "Slowly Rotating Relativistic Stars. I. Equations of Structure",
     Astrophys. J. 150, 1005, Ecs. (40)-(46).
  2. Y. Kojima (1998), "Quasi-toroidal oscillations in slowly rotating relativistic stars",
     Mon. Not. R. Astron. Soc. 293, 49, Ecs. (20) y (25).
  3. K. H. Lockitch & N. Andersson (1999), "r-modes in relativistic stars",
     Astrophys. J. 521, 764.

Este módulo NUNCA es importado por rotation.py; solo se usa en step12_rotation.py para verificar que las
expresiones deducidas desde la acción en rotation.py se reducen exactamente a 0 frente a la literatura.
"""
import sympy as sp

from .model import r, G, Phi0, lam0, a0, n0, eps_of, p_of
from .nonradial import l
from .rotation import Om, sig, m, M_sym, wLT, varpi, j_fn, h0, U, d

e, p = eps_of(n0, a0), p_of(n0)

# ---------------------------------------------------------------- 1. Hartle (1967), Ec. (43)
# Forma en términos de Phi0' + lam0' y densidad de entalpía del fluido (eps + p):
HARTLE_EQ = (
    d(varpi, 2)
    + (4 / r - d(Phi0) - d(lam0)) * d(varpi)
    - 16 * sp.pi * G * (e + p) * sp.exp(2 * lam0) * varpi
)

# Forma autoadjunta de Hartle con j(r) = exp(-(Phi0 + lam0)) y corrección axiónica de pared de dominio:
# (1/r^4) d/dr (r^4 j dvarpi/dr) + (4/r)(dj/dr) varpi = -16 pi G j(r) a0'(r)^2 varpi
HARTLE_J_EQ_AXION = (
    sp.diff(r**4 * j_fn * d(varpi), r) / r**4
    + 4 / r * d(j_fn) * varpi
    + 16 * sp.pi * G * j_fn * d(a0)**2 * varpi
)

# Forma sin axión de Hartle (1967, Ec. 43): (1/r^4) d/dr (r^4 j dvarpi/dr) + (4/r)(dj/dr) varpi = 0
HARTLE_J_EQ_NOAX = (
    sp.diff(r**4 * j_fn * d(varpi), r) / r**4
    + 4 / r * d(j_fn) * varpi
)

# Solución analítica exterior en vacío de Schwarzschild (j = 1): varpi_ext(r) = Omega - 2 G I / r^3
I_sym = sp.Symbol('I', positive=True)
VARPI_EXT = Om - 2 * G * I_sym / r**3

# Condición de regularidad en el origen r -> 0: varpi''(0) = (16 pi G / 5) (eps_c + p_c) varpi_c
wc, ec, pc = sp.symbols('varpi_c varepsilon_c p_c', real=True)
VARPI_DD0_REF = sp.Rational(16, 5) * sp.pi * G * (ec + pc) * wc

# ---------------------------------------------------------------- 2. Kojima (1998, Ecs. 20 y 25)
# Ecuación de vorticidad axial a orden O(Omega):
# l(l+1) sigma (U - h0) - 2 m varpi(r) U + 1/2 (l-1)(l+2) m (Omega - varpi(r)) h0 = 0
KOJIMA_VORTICITY = (
    l * (l + 1) * sig * (U - h0)
    - 2 * m * varpi * U
    + sp.Rational(1, 2) * (l - 1) * (l + 2) * m * (Om - varpi) * h0
)

# Relación de dispersión corrotante de los modos r con arrastre de marcos de Lense-Thirring:
SIGMA_RMODE_REF = 2 * m * varpi / (l * (l + 1))
KAPPA_LM_REF = wLT / Om

# Ecuación maestra radial de Kojima (1998, Ec. 25) antes de eliminar Phi0' + lam0':
KOJIMA_MASTER_RAW = (
    d(h0, 2)
    - (d(Phi0) + d(lam0)) * d(h0)
    - ((l * (l + 1) - 2) / r**2 * sp.exp(2 * lam0) + 2 / r**2 - 2 / r * (d(Phi0) + d(lam0))
       + 16 * sp.pi * G * (e + p) * sp.exp(2 * lam0)) * h0
    + 16 * sp.pi * G * (e + p) * sp.exp(2 * lam0) * U
)

# Ecuación maestra radial de Kojima tras usar las ecuaciones de equilibrio TOV + Klein-Gordon:
KOJIMA_MASTER_TOV = (
    d(h0, 2)
    - 4 * sp.pi * G * r * (sp.exp(2 * lam0) * (e + p) + d(a0)**2) * d(h0)
    - ((l * (l + 1) - 2) / r**2 * sp.exp(2 * lam0) + 2 / r**2
       + 8 * sp.pi * G * (e + p) * sp.exp(2 * lam0) - 8 * sp.pi * G * d(a0)**2) * h0
    + 16 * sp.pi * G * (e + p) * sp.exp(2 * lam0) * U
)

# En términos de la masa encerrada M(r) con e^{-2 lam0} = 1 - 2 G M / r:
KOJIMA_MASTER_M = (
    d(h0, 2)
    - 4 * sp.pi * G * r * ((e + p) / (1 - 2 * G * M_sym / r) + d(a0)**2) * d(h0)
    - ((l * (l + 1) / r**2 - 4 * G * M_sym / r**3 + 8 * sp.pi * G * (e + p)) / (1 - 2 * G * M_sym / r)
       - 8 * sp.pi * G * d(a0)**2) * h0
    + 16 * sp.pi * G * (e + p) / (1 - 2 * G * M_sym / r) * U
)
