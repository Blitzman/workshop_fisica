"""S1: the star integrated with s = ln p as the independent variable, and the axion shot onto the decaying exterior branch.

Interior. The structure is integrated against s = ln p from r0 = 10 m (series about the centre) until p has dropped to `p_frac` times its
central value (default 1e-10). With the condensate the scalar is carried as u = ln(pi - theta) and w = du/dr (log-offset variables): the
shooting parameter is u_c = ln(delta_c), delta_c = pi - theta_c, which stays resolvable when theta_c rounds to pi (delta_c ~ 1e-10 at 20 kHz,
~1e-27 at 50 kHz). Klein-Gordon in these variables:  w' = -w^2 - c w - m_a^2 (1 - x) (sin delta / delta) / (g F),  c F = (2/r)(1 - m/r
- 2 pi r^2 (e - p + 2V)),  g = sqrt(1 - beta cos^2(delta/2)).  dr/ds = p / (dp/dr) is singular at the centre, hence the start at r0 > 0.

Exterior. Beyond p = p_frac p_c the fluid (a ~9 m crust layer) is dropped. The nonlinear vacuum equations (TOV + Klein-Gordon with scalar
backreaction) are integrated in r up to the matching radius r_m = R + match / m_a (default 5 reduced Compton wavelengths). There theta is
~1e-3 and the field obeys the linear equation on Schwarzschild with the enclosed mass M(r_m); its decaying branch is fixed by the logarithmic
derivative S(r) = theta'/theta, obtained from the Riccati equation  S' = m_a^2 / F - 2 (1 - M/r) S / (r F) - S^2  integrated INWARD (stable
for the decaying branch) from r_m + 30 / m_a. The shooting condition is  theta'(r_m) = S(r_m) theta(r_m),  solved for u_c with Brent's
method after a sign bracket. Beyond r_m the tail is continued on that decaying branch (theta' = S theta) with m and Phi still integrated
from the full vacuum equations, up to r_end = R + tail / m_a (default 15). The whole tail [R, r_end] is therefore on the decaying solution;
the previous outward-only shooting stopped at a numerical turning point whose position depended on the bisection tolerance.

The legacy direct solver in r (`star.Star`, outward bisection) is kept as the physics base class and as a reference.

`state_at_r(r)` returns the exact state at arbitrary radii (vectorised): the centre series for r < r0, the ln p solution inverted by Newton
for r0 <= r <= R_stop, and the tail beyond. The tables of dense.py are built from it.
"""
import math

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

from .star import BETA, KGEO, Star, _eos, g, g1, omega_of

DEFAULT_SHOOT_TOL = 1e-12        # Brent xtol in u_c = ln(delta_c)
DEFAULT_MATCH = 5.0              # matching distance beyond R, in reduced Compton wavelengths 1/m_a
DEFAULT_TAIL = 15.0              # table extent beyond R, in 1/m_a
RICCATI_FAR = 30.0               # start of the inward Riccati integration beyond r_m, in 1/m_a
U_LO, U_HI = math.log(1e-300), math.log(math.pi - 0.3)   # search interval for u_c (theta_c from pi - 1e-300 to 0.3)
TURNED = object()                # _interior: the scalar turned back up inside the star (classified +1 by trial)


def series_start(star, r0, thc):
    """State (m, p, Phi, theta, theta') at r0 from the expansion about the centre, to O(r^2) in p, Phi, theta and O(r^3) in m.

    Phi' = (4 pi r / 3)(e + 3p - 2V),  dp/dr = -(e + p) Phi' - q1 theta',  theta' = tpp0 r,  theta''(0) = src(0)/3. Phi(r0) = 0."""
    p = star.p0 * KGEO
    shape = 1 - g(thc) if star.nu_a else 0.0
    slope = g1(thc) if star.nu_a else 0.0
    e, _, _, x, q1 = star._background_matter(p, shape, slope)
    V = star.E0 * shape
    tpp0 = (4 * star.ma**2 / BETA * (x - 1) * slope) / 3 if star.nu_a else 0.0
    m = 4 * np.pi / 3 * (e + V) * r0**3
    dp2 = (e + p) * (4 * np.pi / 3) * (e + 3 * p - 2 * V) + q1 * tpp0           # p = p_c - dp2 r^2 / 2
    return [m, p - 0.5 * dp2 * r0**2, 0.0, thc + tpp0 * r0**2 / 2, tpp0 * r0]


def _physical(z):
    """(r, m, Phi, u, w) -> (r, m, Phi, theta, theta') (arrays of rows or a single state)."""
    y = np.array(z, dtype=float)
    d = np.exp(y[3])
    y[3] = np.pi - d
    y[4] = -d * y[4]
    return y


class _InteriorSolution:
    """The ln p solution seen in physical variables (the interface of a solve_ivp result: t, y, sol)."""

    def __init__(self, raw):
        self.raw, self.t, self.y = raw, raw.t, _physical(raw.y)

    def sol(self, s):
        return _physical(self.raw.sol(s))


class _TailSolution:
    """Vacuum tail on [R_stop, r_end]: nonlinear solution up to r_m, decaying-branch continuation beyond. Rows (m, p = 0, Phi, theta, theta')."""

    def __init__(self, near, far, riccati):
        self.near, self.far, self.ric = near, far, riccati
        self.r_match = near.t[-1]
        self.t = near.t if far is None else np.concatenate([near.t, far.t[1:]])
        self.y = self.sol(self.t)

    def sol(self, r):
        r = np.asarray(r, float)
        scalar = r.ndim == 0
        r = np.atleast_1d(r)
        out = np.empty((5, r.size))
        a = r <= self.r_match
        if a.any():
            out[:, a] = self.near.sol(r[a])
        b = ~a
        if b.any():
            m, Phi, lt = self.far.sol(r[b])
            t = np.exp(lt)
            out[:, b] = np.vstack([m, np.zeros_like(m), Phi, t, self.ric.sol(r[b])[0] * t])
        return out[:, 0] if scalar else out


class StarLogP(Star):
    def __init__(self, nu_a=None, p0=None, e0=None, rtol=1e-8, atol=1e-30, max_step=0.25, shoot_tol=DEFAULT_SHOOT_TOL,
                 r0=0.01, p_frac=1e-10, max_step_s=0.25, match=DEFAULT_MATCH, tail=DEFAULT_TAIL, seed=None):
        """seed: a guess for u_c = ln(pi - theta_c) (e.g. extrapolated from neighbouring nu_a); it only narrows the bracket."""
        from .star import E0, NC, P0, SIGMA_N                      # noqa: F401
        p0 = P0 if p0 is None else p0
        e0 = E0 if e0 is None else e0
        self.rtol, self.atol, self.max_step = rtol, atol, max_step
        self.shoot_tol = DEFAULT_SHOOT_TOL if not shoot_tol else shoot_tol     # None / 0 (legacy 'original precision' spelling): the default
        self.r0, self.p_frac, self.max_step_s = r0, p_frac, max_step_s
        self.match, self.tail, self.seed = match, max(tail, match), seed
        self.nu_a, self.p0, self.E0 = nu_a, p0, e0 * KGEO
        self.nc = e0 / SIGMA_N
        self.ma = omega_of(nu_a) if nu_a else 0.0
        self.fa2 = BETA * self.E0 / (4 * self.ma**2) if nu_a else np.inf
        self.nfev = self.nint = 0
        if nu_a:
            self._shoot()
        else:
            self.thc = 0.0
            self._integrate_noaxion()

    # ------------------------------------------------------------------ right-hand sides
    def _rhs_s(self, s, z):
        """Structure in s = ln p, physical theta (no condensate, or any caller that wants theta directly)."""
        r, m, Phi, t, tp = z
        self.nfev += 1
        p = self.p0 * KGEO * math.exp(s)
        dm, dp, dPhi, dt, dtp = self.rhs(r, [m, p, Phi, t, tp])
        if dp >= 0:
            raise RuntimeError(f'StarLogP: dp/dr >= 0 at r = {r:.6g} km (s = {s:.4g}): p is not a valid coordinate there')
        drds = p / dp
        return [drds, dm * drds, dPhi * drds, dt * drds, dtp * drds]

    def _rhs_su(self, s, z):
        """Structure in s = ln p with the log-offset scalar (u, w). Same physics as Star.rhs; g and g1 are evaluated from delta itself."""
        r, m, Phi, u, w = z
        self.nfev += 1
        d = math.exp(u)
        t, tp = math.pi - d, -d * w
        p = self.p0 * KGEO * math.exp(s)
        gg = math.sqrt(1 - BETA * math.cos(d / 2)**2)              # g(theta) with theta = pi - delta
        shape = 1 - gg
        slope = -BETA * math.sin(d) / (4 * gg)                      # g1(theta): sin(theta) = sin(delta)
        F = 1 - 2 * m / r
        V = self.E0 * shape
        e, _, _, x, q1 = self._background_matter(p, shape, slope)
        kin = self.fa2 * tp**2 / 2 * F
        dPhi = (m + 4 * np.pi * r**3 * (p - V + kin)) / (r**2 * F)
        dm = 4 * np.pi * r**2 * (e + V + kin)
        dp = -(e + p) * dPhi - q1 * tp
        if dp >= 0:
            raise RuntimeError(f'StarLogP: dp/dr >= 0 at r = {r:.6g} km (s = {s:.4g}): p is not a valid coordinate there')
        c = 2 / r * (1 - m / r - 2 * np.pi * r * r * (e - p + 2 * V)) / F
        sinc = math.sin(d) / d if d else 1.0
        dw = -w * w - c * w - self.ma**2 * (1 - x) * sinc / (gg * F)
        drds = p / dp
        return [drds, dm * drds, dPhi * drds, w * drds, dw * drds]

    def _rhs_r(self, r, y):
        self.nfev += 1
        return self.rhs(r, y)

    # ------------------------------------------------------------------ integrations
    def _integrate_noaxion(self):
        self.nint += 1
        m0, p_r0, Phi0, t0, tp0 = series_start(self, self.r0, 0.0)
        s0, s1 = math.log(p_r0 / (self.p0 * KGEO)), math.log(self.p_frac)
        s_in = solve_ivp(self._rhs_s, (s0, s1), [self.r0, m0, Phi0, t0, tp0], method='DOP853', rtol=self.rtol, atol=self.atol,
                         dense_output=True, max_step=self.max_step_s)
        if s_in.status < 0:
            raise RuntimeError(f'StarLogP (structure): {s_in.message}')
        self.sol_in = s_in
        self.R_stop, self.M_stop = s_in.y[0, -1], s_in.y[1, -1]

    def _interior(self, u_c, keep=False):
        """Structure + log-offset scalar from r0 to p = p_frac p_c. None if theta crosses below -1e-3 inside the star; TURNED if theta turns
        back up inside the star (w = u' falls through zero: delta would then collapse to 0 and u, w -> -inf before the surface; this happens
        in the crust at large m_a, e.g. 100 kHz, where the crust is thicker than 1/m_a)."""
        thc = math.pi - math.exp(u_c)
        m, p, Phi, _, _ = series_start(self, self.r0, thc)
        d = math.exp(u_c)
        x = self._background_matter(self.p0 * KGEO, 0.0, 0.0)[3]
        gg = math.sqrt(1 - BETA * math.cos(d / 2)**2)
        w = -self.r0 * self.ma**2 * (1 - x) * (math.sin(d) / d if d else 1.0) / (3 * gg)    # u'(r) = u''(0) r
        u = u_c + w * self.r0 / 2

        def cross(s, z):                                              # theta = -1e-3  <=>  u = ln(pi + 1e-3)
            return math.log(math.pi + 1e-3) - z[3]

        def turn(s, z):                                               # theta' = 0 from below  <=>  w = 0, falling
            return z[4]
        cross.terminal = turn.terminal = True
        turn.direction = -1
        sol = solve_ivp(self._rhs_su, (math.log(p / (self.p0 * KGEO)), math.log(self.p_frac)), [self.r0, m, Phi, u, w], events=[cross, turn],
                        method='DOP853', rtol=self.rtol, atol=self.atol, max_step=self.max_step_s, dense_output=keep)
        if sol.status < 0:
            raise RuntimeError(f'StarLogP (structure): {sol.message}')
        if sol.status == 1:
            return None if sol.t_events[0].size else TURNED
        return sol

    def _riccati(self, mass, rm):
        """S = theta'/theta on the decaying branch of the linear exterior equation (Schwarzschild, mass M), from r_m + 30/m_a inward to r_m."""
        r_far = rm + max(RICCATI_FAR, self.tail - self.match + 5.0) / self.ma
        ma2 = self.ma**2

        def f(r, S):
            F = 1 - 2 * mass / r
            return ma2 / F - 2 * (1 - mass / r) * S / (r * F) - S * S
        sol = solve_ivp(f, (r_far, rm), [-self.ma - (1 + mass * self.ma) / r_far], method='DOP853', rtol=1e-10, atol=1e-12, dense_output=True)
        if not sol.success:
            raise RuntimeError(f'StarLogP (exterior Riccati): {sol.message}')
        return sol

    def trial(self, u_c, keep=False):
        """One shot. Returns (sign, residual): sign +1 when theta turns back up (u_c too small), -1 when it crosses zero; residual =
        theta' - S theta at r_m (None when the shot does not reach r_m, the sign alone then brackets)."""
        self.nint += 1
        interior = self._interior(u_c, keep)
        if interior is None:
            return -1, None
        if interior is TURNED:
            return 1, None
        R, m, Phi, t, tp = _physical(interior.y[:, -1])
        if tp >= 0 and t > 0:                                         # already leaving the decaying branch at the surface
            return 1, None
        rm = R + self.match / self.ma

        def cross(r, y):
            return y[3] + 1e-3

        def turn(r, y):
            return y[4]
        cross.terminal = turn.terminal = True
        cross.direction, turn.direction = -1, 1
        near = solve_ivp(self._rhs_r, (R, rm), [m, 0.0, Phi, t, tp], events=[cross, turn], method='DOP853', rtol=self.rtol, atol=self.atol,
                         max_step=self.max_step, dense_output=keep)
        if near.status < 0:
            raise RuntimeError(f'StarLogP (exterior): {near.message}')
        if near.status == 1:
            return (-1 if near.t_events[0].size else 1), None
        mass, _, Phi_m, t_m, tp_m = near.y[:, -1]
        ric = self._riccati(mass, rm)
        residual = tp_m - ric.y[0, -1] * t_m
        if keep:
            self.u_c, self.thc = u_c, thc_of(u_c)
            self.sol_in = _InteriorSolution(interior)
            self.R_stop, self.M_stop = R, m
            self.theta_match = t_m
            self.decay_mismatch = residual / (self.ma * abs(t_m))
            self.sol_out = self._continue_tail(near, ric)
        return int(np.sign(residual)), residual

    def _continue_tail(self, near, ric):
        """Tail beyond r_m on the decaying branch: d ln(theta)/dr = S(r); m and Phi from the full vacuum equations."""
        rm, rend = near.t[-1], self.R_stop + self.tail / self.ma
        if rend <= rm * (1 + 1e-12):
            return _TailSolution(near, None, ric)
        m, _, Phi, t, _ = near.y[:, -1]
        if t <= 0:
            raise RuntimeError('StarLogP: theta <= 0 at the matching radius')

        def f(r, z):
            m, Phi, lt = z
            S = ric.sol(r)[0]
            t = math.exp(lt)
            dm, _, dPhi, _, _ = self.rhs(r, [m, 0.0, Phi, t, S * t])
            return [dm, dPhi, S]
        far = solve_ivp(f, (rm, rend), [m, Phi, math.log(t)], method='DOP853', rtol=self.rtol, atol=self.atol, max_step=self.max_step,
                        dense_output=True)
        if not far.success:
            raise RuntimeError(f'StarLogP (tail continuation): {far.message}')
        return _TailSolution(near, far, ric)

    def _shoot(self):
        cache = {}

        def sample(u):
            if u not in cache:
                cache[u] = self.trial(u)
            return cache[u]
        lo, hi = U_LO, U_HI
        if self.seed is not None:
            for width in (0.25, 1.0, 4.0, 16.0):
                a, b = max(lo, self.seed - width), min(hi, self.seed + width)
                if a < b and sample(a)[0] > 0 > sample(b)[0]:
                    lo, hi = a, b
                    break
        if not sample(lo)[0] > 0 > sample(hi)[0]:
            # Both ends undershoot (or the delta = 1e-300 end overshoots): this does NOT distinguish "no condensate" from a condensate
            # with theta_c < 0.3, outside the search interval (or delta_c < 1e-300).
            raise ValueError('sin cambio de signo en theta_c en [0.3, pi - 1e-300]: o no hay condensado, o tiene theta_c < 0.3 '
                             '(fuera del intervalo de búsqueda), o delta_c < 1e-300')
        for _ in range(200):                       # bisection on the sign until both ends reach r_m (Brent needs finite residuals)
            if sample(lo)[1] is not None and sample(hi)[1] is not None:
                break
            mid = 0.5 * (lo + hi)
            if sample(mid)[0] > 0:
                lo = mid
            else:
                hi = mid
        else:
            raise RuntimeError('StarLogP: no finite matching residual on both sides of the bracket')
        self.bracket_trials = self.nint

        def residual(u):
            v = sample(u)[1]
            if v is None:
                raise RuntimeError('StarLogP: a Brent trial stopped before the matching radius')
            return v
        root = brentq(residual, lo, hi, xtol=self.shoot_tol, rtol=4 * np.finfo(float).eps)
        if self.trial(root, keep=True)[1] is None:
            raise RuntimeError('StarLogP: the final trial stopped before the matching radius')

    def _integrate(self, thc, keep=False):
        """Compatibility: one shot at theta_c (+1 / -1 as Star._integrate; 0 without condensate)."""
        if not self.nu_a:
            self._integrate_noaxion()
            return 0
        return self.trial(math.log(math.pi - thc), keep=keep)[0]

    # ---- exact state at arbitrary radii (needs the dense solution, which __init__ leaves in place)
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
        shape = 1 - g(thc) if self.nu_a else 0.0
        slope = g1(thc) if self.nu_a else 0.0
        e, _, _, x, q1 = self._background_matter(pc, shape, slope)
        V = self.E0 * shape
        tpp0 = (4 * self.ma**2 / BETA * (x - 1) * slope) / 3 if self.nu_a else 0.0
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
    def delta_c(self):
        return math.exp(self.u_c) if self.nu_a else math.pi

    @property
    def r_match(self):
        return self.sol_out.r_match if self.nu_a else self.R_stop

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


def thc_of(u_c):
    return math.pi - math.exp(u_c)


class SeedTracker:
    """Guesses u_c for the next nu_a from the stars already solved (linear extrapolation of the last two, or the last one)."""

    def __init__(self):
        self.done = []

    def add(self, nu, u_c):
        self.done.append((nu, u_c))

    def guess(self, nu):
        if not self.done:
            return None
        if len(self.done) == 1 or self.done[-1][0] == self.done[-2][0]:
            return self.done[-1][1]
        (n1, u1), (n2, u2) = self.done[-2:]
        return u2 + (u2 - u1) * (nu - n2) / (n2 - n1)
