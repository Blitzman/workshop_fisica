"""S1 (experimental): the star integrated with s = ln p as the independent variable.

The structure z = (r, m, Phi, theta, theta') is integrated against s = ln p from r0 = 10 m (series about the centre) until p has dropped to
`p_frac` times its central value (default 1e-10); then the (fluid-free) axion tail is integrated in r, with the same events and the same
shooting on theta_c as `Star`. Physics (matter, rhs) is inherited from `Star`; only the integration variable, the start and the stop change.

dr/ds = p / (dp/dr) is singular at the centre (dp/dr ~ r), hence the start at r0 > 0.
Beyond p = p_frac p_c the fluid (a ~9 m crust layer with eps ~ 1e-6 eps_c) is dropped: the vacuum tail starts there.

`state_at_r(r)` returns the exact state at arbitrary radii (vectorised): the centre series for r < r0, the ln p solution inverted by Newton
for r0 <= r <= R_stop, and the vacuum tail beyond. The tables of dense.py are built from it.
"""
import math

import numpy as np
from scipy.integrate import solve_ivp

from .star import BETA, KGEO, Star, _eos, g, g1, omega_of


def series_start(star, r0, thc):
    """State (m, p, Phi, theta, theta') at r0 from the expansion about the centre, to O(r^2) in p, Phi, theta and O(r^3) in m.

    Phi' = (4 pi r / 3)(e + 3p - 2V),  dp/dr = -(e + p) Phi' - q1 theta',  theta' = tpp0 r,  theta''(0) = src(0)/3. Phi(r0) = 0."""
    p = star.p0 * KGEO
    e, nb, x, cs2, q1, q2 = star.matter(p, thc)
    V = star.E0 * (1 - g(thc)) if star.nu_a else 0.0
    tpp0 = (4 * star.ma**2 / BETA * (x - 1) * g1(thc)) / 3 if star.nu_a else 0.0
    m = 4 * np.pi / 3 * (e + V) * r0**3
    dp2 = (e + p) * (4 * np.pi / 3) * (e + 3 * p - 2 * V) + q1 * tpp0           # p = p_c - dp2 r^2 / 2
    return [m, p - 0.5 * dp2 * r0**2, 0.0, thc + tpp0 * r0**2 / 2, tpp0 * r0]


class StarLogP(Star):
    def __init__(self, nu_a=None, p0=None, e0=None, rtol=1e-8, atol=1e-30, max_step=0.25, shoot_tol=1e-12,
                 r0=0.01, p_frac=1e-10, max_step_s=0.25):
        from .star import E0, NC, P0, SIGMA_N                      # noqa: F401
        p0 = P0 if p0 is None else p0
        e0 = E0 if e0 is None else e0
        self.rtol, self.atol, self.max_step, self.shoot_tol = rtol, atol, max_step, shoot_tol
        self.r0, self.p_frac, self.max_step_s = r0, p_frac, max_step_s
        self.nu_a, self.p0, self.E0 = nu_a, p0, e0 * KGEO
        self.nc = e0 / SIGMA_N
        self.ma = omega_of(nu_a) if nu_a else 0.0
        self.fa2 = BETA * self.E0 / (4 * self.ma**2) if nu_a else np.inf
        self.nfev = self.nint = 0
        if nu_a:
            self._shoot()
        else:
            self.thc = 0.0
            self._integrate(0.0, keep=True)

    def _rhs_s(self, s, z):
        r, m, Phi, t, tp = z
        self.nfev += 1
        p = self.p0 * KGEO * math.exp(s)
        dm, dp, dPhi, dt, dtp = self.rhs(r, [m, p, Phi, t, tp])
        if dp >= 0:
            raise RuntimeError(f'StarLogP: dp/dr >= 0 at r = {r:.6g} km (s = {s:.4g}): p is not a valid coordinate there')
        drds = p / dp
        return [drds, dm * drds, dPhi * drds, dt * drds, dtp * drds]

    def _rhs_r(self, r, y):
        self.nfev += 1
        return self.rhs(r, y)

    def _integrate(self, thc, keep=False, rmax=None):
        """+1 if theta shoots up, -1 if it crosses zero (as Star._integrate)."""
        self.nint += 1
        rmax = rmax or (60.0 + (40 / self.ma if self.nu_a else 0.0))
        m0, p_r0, Phi0, t0, tp0 = series_start(self, self.r0, thc)
        pc = self.p0 * KGEO
        s0, s1 = math.log(p_r0 / pc), math.log(self.p_frac)
        up = lambda s, z: (z[3] - np.pi - 1e-3) if self.nu_a else 1.0
        down = lambda s, z: z[3] + 1e-3 if self.nu_a else 1.0
        for ev in (up, down):
            ev.terminal = True
        s_in = solve_ivp(self._rhs_s, (s0, s1), [self.r0, m0, Phi0, t0, tp0], events=[up, down], method='DOP853', rtol=self.rtol,
                         atol=self.atol, dense_output=keep, max_step=self.max_step_s)
        if s_in.status < 0:
            raise RuntimeError(f'StarLogP._integrate (structure): {s_in.message}')
        if keep:
            self.sol_in = s_in
        if not self.nu_a:
            self.R_stop, self.M_stop = s_in.y[0, -1], s_in.y[1, -1]
            return 0
        if s_in.t_events[1].size:
            return -1
        if s_in.t_events[0].size:
            return +1
        # vacuum tail in r (the same events as Star, `turn` included)
        r_st, m, Phi, t, tp = s_in.y[:, -1]
        turn = lambda r, y: y[4] if (self.nu_a and y[1] / KGEO < _eos.p_min * 1.0001) else -1.0
        up_r = lambda r, y: y[3] - np.pi - 1e-3
        down_r = lambda r, y: y[3] + 1e-3
        for ev in (up_r, down_r):
            ev.terminal = True
        turn.terminal, turn.direction = True, 1
        s_out = solve_ivp(self._rhs_r, (r_st, rmax), [m, 0.0, Phi, t, tp], events=[up_r, down_r, turn], method='DOP853',
                          rtol=self.rtol, atol=self.atol, dense_output=keep, max_step=self.max_step)
        if s_out.status < 0:
            raise RuntimeError(f'StarLogP._integrate (tail): {s_out.message}')
        if keep:
            self.sol_out = s_out
            self.R_stop, self.M_stop = r_st, m
        return -1 if s_out.t_events[1].size else +1

    # ---- derived quantities (compare with Star.R, Star.M, Star.Phi_shift, defined at p = p_frac p_c here)
    # ---- exact state at arbitrary radii (needs the dense solution: keep=True, which _shoot/__init__ leave in place)
    def _dp_dr(self, r, m, p, t, tp):
        """dp/dr inside the star (vectorised): same formulas as Star.matter / Star.rhs."""
        pm = p / KGEO
        e0, nb = _eos.eps(pm) * KGEO, _eos.nb(pm)
        x, F = nb / self.nc, 1 - 2 * m / r
        if self.nu_a:
            gg, gg1 = g(t), g1(t)
            e, q1 = e0 - self.E0 * x * (1 - gg), self.E0 * x * gg1
            V, kin = self.E0 * (1 - gg), self.fa2 * tp**2 / 2 * F
        else:
            e, q1, V, kin = e0, 0.0, 0.0, 0.0
        dPhi = (m + 4 * np.pi * r**3 * (p - V + kin)) / (r**2 * F)
        return -(e + p) * dPhi - q1 * tp

    def _centre_series(self, r):
        """(m, p, Phi, theta, theta') for r < r0 from the expansion about the centre (Phi(r0) = 0)."""
        pc = self.p0 * KGEO
        thc = self.thc
        e, nb, x, cs2, q1, q2 = self.matter(pc, thc)
        V = self.E0 * (1 - g(thc)) if self.nu_a else 0.0
        tpp0 = (4 * self.ma**2 / BETA * (x - 1) * g1(thc)) / 3 if self.nu_a else 0.0
        dp2 = (e + pc) * (4 * np.pi / 3) * (e + 3 * pc - 2 * V) + q1 * tpp0
        return np.vstack([4 * np.pi / 3 * (e + V) * r**3, pc - 0.5 * dp2 * r**2, (2 * np.pi / 3) * (e + 3 * pc - 2 * V) * (r**2 - self.r0**2),
                          thc + tpp0 * r**2 / 2, tpp0 * r])

    def state_at_r(self, r):
        """Rows (m, p, Phi, theta, theta') at r (array or scalar; Phi without the shift)."""
        r = np.asarray(r, float)
        scalar = r.ndim == 0
        r = np.atleast_1d(r)
        out = np.empty((5, r.size))
        pc = self.p0 * KGEO
        if not hasattr(self, '_rfine'):
            self._sfine = np.linspace(self.sol_in.t[0], self.sol_in.t[-1], 40001)
            self._rfine = self.sol_in.sol(self._sfine)[0]
        cen, body, tail = r < self.r0, (r >= self.r0) & (r <= self.R_stop), r > self.R_stop
        if cen.any():
            out[:, cen] = self._centre_series(r[cen])
        if body.any():
            rr = r[body]
            # r^2 is nearly linear in s near the centre (r ~ sqrt(-s)): Newton on r^2, started from the interpolation of r^2(s)
            s = np.interp(rr**2, self._rfine**2, self._sfine)
            for _ in range(6):
                z = self.sol_in.sol(s)
                p = pc * np.exp(s)
                dr_ds = p / self._dp_dr(z[0], z[1], p, z[3], z[4])
                s = np.clip(s - (z[0]**2 - rr**2) / (2 * z[0] * dr_ds), self.sol_in.t[-1], self.sol_in.t[0])
            z = self.sol_in.sol(s)
            out[:, body] = np.vstack([z[1], pc * np.exp(s), z[2], z[3], z[4]])
        if tail.any():
            if self.nu_a:
                out[:, tail] = self.sol_out.sol(r[tail])
            else:                                               # Schwarzschild vacuum
                rt = r[tail]
                Phi_R = self.sol_in.y[2, -1]
                out[:, tail] = np.vstack([np.full_like(rt, self.M_stop), np.zeros_like(rt),
                                          Phi_R + 0.5 * np.log((1 - 2 * self.M_stop / rt) / (1 - 2 * self.M_stop / self.R_stop)),
                                          np.zeros_like(rt), np.zeros_like(rt)])
        return out[:, 0] if scalar else out

    @property
    def R(self):
        return self.R_stop

    @property
    def M(self):
        return self.M_stop

    @property
    def r_end(self):
        return self.sol_out.t[-1] if self.nu_a else self.R_stop

    @property
    def M_total(self):
        return self.sol_out.y[0, -1] if self.nu_a else self.M_stop

    @property
    def Phi_shift(self):
        if not hasattr(self, '_Ps'):
            if self.nu_a:
                rend, (m, _, Phi) = self.sol_out.t[-1], self.sol_out.y[:3, -1]
            else:
                rend, m, Phi = self.R_stop, self.M_stop, self.sol_in.y[2, -1]
            self._Ps = 0.5 * np.log(1 - 2 * m / rend) - Phi
        return self._Ps

    def at_s(self, s):
        """(r, m, Phi + shift, theta, theta') at s = ln(p/p_c) (structure part)."""
        r, m, Phi, t, tp = self.sol_in.sol(s)
        return r, m, Phi + self.Phi_shift, t, tp

    def at_r_tail(self, r):
        m, p, Phi, t, tp = self.sol_out.sol(r)
        return m, Phi + self.Phi_shift, t, tp
