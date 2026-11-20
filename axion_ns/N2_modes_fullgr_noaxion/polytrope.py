"""N2: analytic polytropes p = K rho^2 and their TOV background (moved from axion_ns/fmodes.py: EOS, Star).

Units G = c = 1 and K = 1 (polytropes are scale invariant). `kind`: 'rest' (rho = rest-mass density) or 'energy' (rho = energy
density; the convention of Lindblom-Mendell-Ipser and Yoshida-Kojima). Differences with the original: `integrate` takes the
numerical system object of systems.py and explicit tolerances, and counts right-hand-side calls (`nfev`).
"""
import numpy as np
import sympy as sp
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

_x = sp.Symbol('x', positive=True)
EOS = {'rest': _x + _x**2, 'energy': _x / (1 - _x)}     # eps(n) with m_b = 1 and K = 1


def eos_num(kind):
    e = EOS[kind]
    p = sp.expand(_x * sp.diff(e, _x) - e)
    f = sp.lambdify(_x, [e, p, sp.diff(p, _x), sp.diff(p, _x, 2), sp.diff(e, _x)], 'numpy')
    n_of_p = (lambda q: np.sqrt(q)) if kind == 'rest' else (lambda q: np.sqrt(q) / (1 + np.sqrt(q)))
    return f, n_of_p


class Star:
    def __init__(self, kind, compactness=None, mass=None):
        """By compactness M/R or by mass M (units of K^{1/2}, stable branch)."""
        self.kind = kind
        self.eos, self.n_of_p = eos_num(kind)
        self.nfev = 0
        lo, hi = 1e-4, 0.3 if kind == 'rest' else 0.2
        target = (lambda nc: self._solve(nc)[0] - compactness) if mass is None else (lambda nc: self._solve(nc)[1] - mass)
        self.nc = brentq(target, lo, hi if mass is None else 0.15, xtol=1e-14)
        self.MR, self.M, self.R, self.Phic = self._solve(self.nc)

    def bg(self, rr, m, p, Phi):
        """(n, n', n'', Phi, Phi', Phi'', lambda, lambda') and (m', p', Phi')."""
        n = self.n_of_p(max(p, 0.0))
        e, _, dpdn, d2p, dedn = self.eos(n)
        D = rr * (rr - 2 * m)
        Pp = (m + 4 * np.pi * rr**3 * p) / D
        mp, pp = 4 * np.pi * rr**2 * e, -(e + p) * Pp
        Npr = mp + 12 * np.pi * rr**2 * p + 4 * np.pi * rr**3 * pp
        Dp = 2 * rr - 2 * m - 2 * rr * mp
        Ppp = (Npr * D - (m + 4 * np.pi * rr**3 * p) * Dp) / D**2
        lam, lamp = -0.5 * np.log(1 - 2 * m / rr), (rr * mp - m) / D
        n1 = pp / dpdn if dpdn > 0 else 0.0
        ppp = -(dedn * n1 + pp) * Pp - (e + p) * Ppp
        n2 = (ppp - d2p * n1**2) / dpdn if dpdn > 0 else 0.0
        return (n, n1, n2, Phi, Pp, Ppp, lam, lamp), (mp, pp, Pp)

    def _start(self, nc, r0):
        e, p = self.eos(nc)[:2]
        return [4 * np.pi / 3 * e * r0**3, p - 2 * np.pi / 3 * (e + p) * (e + 3 * p) * r0**2, 0.0]

    def _solve(self, nc):
        r0 = 1e-6
        surf = lambda rr, y: y[1] - 1e-14 * self.eos(nc)[1]
        surf.terminal = True
        s = solve_ivp(lambda rr, y: self.bg(rr, *y)[1], (r0, 50), self._start(nc, r0), events=surf,
                      method='DOP853', rtol=1e-12, atol=1e-30)
        if s.status < 0:                                # H32
            raise RuntimeError(f'polytrope: {s.message}')
        R = s.t[-1]
        M = s.y[0, -1]
        Phic = 0.5 * np.log(1 - 2 * M / R) - s.y[2, -1]
        return M / R, M, R, Phic

    def integrate(self, system, y0, omega, ll, r0=None, r1=None, rtol=1e-11, atol=1e-14):
        """Integrates background + perturbations from r0 to r1 (near the surface); `system` is a systems.Sys."""
        r0 = r0 or 1e-4 * self.R
        r1 = r1 or self.R * (1 - 1e-7)
        b0 = self._start(self.nc, r0)
        b0[2] = self.Phic
        ay = system.ay

        def rhs(rr, y):
            self.nfev += 1
            B, dB = self.bg(rr, y[0], y[1], y[2])
            if ay is not None:
                return [*dB, *ay(rr, omega, ll, *B, *y[3:])]
            return np.concatenate([dB, system.mat(rr, omega, ll, *B) @ y[3:]])
        # At the last point of the interval (r/R = 1 - 1e-7) the polytrope surface makes a coefficient divide by zero: the solver sees a NaN error
        # estimate, rejects that step and retries (identical in the legacy code; roots agree to 1e-15). Real failures still raise (status < 0).
        with np.errstate(divide='ignore', invalid='ignore'):
            s = solve_ivp(rhs, (r0, r1), np.concatenate([b0, y0]), method='DOP853', rtol=rtol, atol=atol)
        if s.status < 0:                                # H32
            raise RuntimeError(f'polytrope: {s.message}')
        return s.y[:, -1]

    def center(self, omega, ll, K0, W0):
        """Values at the centre (LD): H1(0), K(0), W(0), X(0)."""
        e, p = self.eos(self.nc)[:2]
        X0 = (e + p) * np.exp(self.Phic) * ((4 * np.pi / 3 * (e + 3 * p)
                                             - omega**2 * np.exp(-2 * self.Phic) / ll) * W0 + K0 / 2)
        return np.array([(2 * ll * K0 + 16 * np.pi * (e + p) * W0) / (ll * (ll + 1)), K0, W0, X0])
