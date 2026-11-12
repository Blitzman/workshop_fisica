"""T_s - T_b relations of Potekhin, Chabrier & Yakovlev 1997 (arXiv:astro-ph/9706148), Appendix A.3.

T_b [K] at the bottom of the heat-blanketing envelope, g14 = surface gravity / 1e14 cm s^-2. Returns T_s [K].
"""
import numpy as np


def ts_iron(Tb, g14):
    tb9 = Tb / 1e9
    tstar = np.sqrt(7 * tb9 * np.sqrt(g14))                    # (A6)
    zeta = tb9 - tstar / 1e3
    return 1e6 * (g14 * ((7 * zeta)**2.25 + (zeta / 3)**1.25))**0.25   # (A7)


def ts_accreted(Tb, g14):
    return 1e6 * (g14 * (18.1 * Tb / 1e9)**2.42)**0.25         # (A8)


def ts_partial(Tb, g14, eta):
    """eta = g14^2 dM/M, (A5); (A9)-(A10)."""
    a = (1.2 + (5.3e-6 / eta)**0.38) * (Tb / 1e9)**(5 / 3)
    return ((a * ts_iron(Tb, g14)**4 + ts_accreted(Tb, g14)**4) / (a + 1))**0.25
