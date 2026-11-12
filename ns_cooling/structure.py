"""Radial profile of a TOV star: r, m, P, Phi, rho, n_b (cgs; n_b in fm^-3)."""
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

from ns_cooling.tov import C, G, MSUN, star


def profile(tab, rho_c, rho_surf=1e6, rtol=1e-10, N=4000):
    Pc, Ps = tab.P_of_rho(rho_c), tab.P_of_rho(rho_surf)

    def rhs(r, y):
        P, m, _ = y
        rho = tab.rho_of_P(max(P, 0.5 * Ps))
        dP = -G * (rho + P / C**2) * (m + 4 * np.pi * r**3 * P / C**2) / (r * (r - 2 * G * m / C**2))
        return [dP, 4 * np.pi * r**2 * rho, -dP / (P + rho * C**2)]

    def surf(r, y):
        return y[0] - Ps
    surf.terminal = True
    r0 = 1.0
    sol = solve_ivp(rhs, (r0, 5e6), [Pc, 4 / 3 * np.pi * r0**3 * rho_c, 0.0], events=surf, rtol=rtol,
                    atol=[Ps * 1e-6, 1e10, 1e-12], method='DOP853', dense_output=True)
    assert sol.status == 1, (rho_c, sol.message)
    R = sol.t[-1]
    # uniform in r**2 near the centre is not needed; refine towards the surface, where the profile is steep
    r = R * (1 - np.linspace(1, 0, N)**2 * (1 - r0 / R))
    P, m, phi = sol.sol(r)
    M = m[-1]
    phi = phi - phi[-1] + 0.5 * np.log(1 - 2 * G * M / (R * C**2))      # match to Schwarzschild at the surface
    P = np.maximum(P, Ps)
    return {'r': r, 'm': m, 'P': P, 'phi': phi, 'rho': tab.rho_of_P(P), 'n': tab.n_of_P(P), 'M': M / MSUN, 'R': R}


def with_mass(tab, M, lo=14.4, hi=15.3, **kw):
    lr = brentq(lambda x: star(tab, 10**x)[0] - M, lo, hi, xtol=1e-10)
    return profile(tab, 10**lr, **kw)


def volume_integral(pr, f, mask):
    """Integral of f over proper volume, 4 pi r^2 dr / sqrt(1 - 2Gm/rc^2), restricted to mask."""
    w = 4 * np.pi * pr['r']**2 / np.sqrt(1 - 2 * G * pr['m'] / (pr['r'] * C**2))
    y = np.where(mask, f * w, 0.0)
    return float(np.sum(0.5 * (y[1:] + y[:-1]) * np.diff(pr['r'])))
