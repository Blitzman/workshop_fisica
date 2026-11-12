"""TOV integration for a barotropic table rho(P). cgs units."""
import numpy as np
from scipy.integrate import solve_ivp
from scipy.interpolate import CubicSpline

G, C, MSUN = 6.6743e-8, 2.99792458e10, 1.98841e33


class Table:
    def __init__(self, n, rho, P):
        ok = np.concatenate([[True], np.diff(P) > 0])
        self.n, self.rho, self.P = n[ok], rho[ok], P[ok]
        self._lr = CubicSpline(np.log(self.P), np.log(self.rho))
        self._ln = CubicSpline(np.log(self.P), np.log(self.n))
        self._lp = CubicSpline(np.log(self.rho), np.log(self.P))

    def rho_of_P(self, P):
        return np.exp(self._lr(np.log(P)))

    def n_of_P(self, P):
        return np.exp(self._ln(np.log(P)))

    def P_of_rho(self, rho):
        return np.exp(self._lp(np.log(rho)))


def star(tab, rho_c, rho_surf=1e6, rtol=1e-10, newtonian=False):
    """Mass [Msun] and radius [km] for central mass-energy density rho_c [g/cm^3]."""
    Pc, Ps = tab.P_of_rho(rho_c), tab.P_of_rho(rho_surf)

    def rhs(r, y):
        P, m = y
        rho = tab.rho_of_P(max(P, 0.5 * Ps))      # the last trial step may overshoot the surface
        if newtonian:                              # negative control only
            return [-G * rho * m / r**2, 4 * np.pi * r**2 * rho]
        dP = -G * (rho + P / C**2) * (m + 4 * np.pi * r**3 * P / C**2) / (r * (r - 2 * G * m / C**2))
        return [dP, 4 * np.pi * r**2 * rho]

    def surf(r, y):
        return y[0] - Ps
    surf.terminal = True
    r0 = 1.0
    sol = solve_ivp(rhs, (r0, 5e6), [Pc, 4 / 3 * np.pi * r0**3 * rho_c], events=surf, rtol=rtol, atol=[Ps * 1e-6, 1e10],
                    method='DOP853')
    assert sol.status == 1, (rho_c, sol.message)
    return sol.y[1, -1] / MSUN, sol.t[-1] / 1e5
