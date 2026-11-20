"""N4: the background star as seen by the mode solver, rebuilt from S1's tables (no import from S1). Duplicate of N1/N3 star_table.py on purpose
(modules never import each other; tests/test_physics_crosscheck.py checks the formulas against S1).

Inputs (written by S1): eos_bsk26.csv and star_*.csv (+ .meta.json). The star table is the dense STATE table
(ln p, m, Phi, theta, theta' on ~7,000 nodes, finer in r near the surface); interpolated with ONE cubic spline lookup per radius (not-a-knot;
two splines, split at the node R_stop = r_k where the star hands over to the fluid-free tail; the tail has no ln p). This is the interpolation whose
error S1 measures and records. The EOS is rebuilt from its table exactly as bsk26.EOS builds it, but
ln eps, ln n and d ln eps / d ln p are evaluated together as one stacked polynomial. Everything else (energy density, sound
speed, couplings, Phi', lambda', m') is exact arithmetic from the state: `_matter` and `rhs` are copies of
axion_ns/axstar.py::Star.matter / Star.rhs, validated against S1's Star in tests/ (test_physics_crosscheck).

Per right-hand-side call this costs about 9 microseconds (the Hermite/CSV version of this class cost about 60), see MIGRATION.md.
"""
import math

import numpy as np
from scipy.interpolate import CubicSpline

from common import artifacts as A


class _FastSpline:
    """Cubic spline (not-a-knot) with scalar evaluation of all columns: one searchsorted + Horner."""

    def __init__(self, x, Y):
        cs = CubicSpline(x, Y, axis=0)
        self.x = np.ascontiguousarray(x)
        self.c = np.ascontiguousarray(cs.c.transpose(1, 0, 2))          # (n-1, 4, ncols)
        self.n = len(self.x)

    def __call__(self, r):
        i = int(np.searchsorted(self.x, r, 'right')) - 1
        i = 0 if i < 0 else (self.n - 2 if i > self.n - 2 else i)
        dx = r - self.x[i]
        c = self.c[i]
        return ((c[0] * dx + c[1]) * dx + c[2]) * dx + c[3]


class _EOS:
    """ln eps, ln n and d ln eps / d ln p as one stacked polynomial evaluation; same splines as bsk26.EOS."""

    def __init__(self, n, e, p):
        lp = np.log(p)
        le, ln = CubicSpline(lp, np.log(e)), CubicSpline(lp, np.log(n))
        dc = np.zeros_like(le.c)
        dc[1], dc[2], dc[3] = 3 * le.c[0], 2 * le.c[1], le.c[2]
        self.c = np.ascontiguousarray(np.stack([le.c, ln.c, dc], axis=-1).transpose(1, 0, 2))      # (n-1, 4, 3)
        self.x, self.nx = np.ascontiguousarray(lp), len(lp)
        self.p_min = p[0]

    def __call__(self, lp):
        i = int(np.searchsorted(self.x, lp, 'right')) - 1
        i = 0 if i < 0 else (self.nx - 2 if i > self.nx - 2 else i)
        dx = lp - self.x[i]
        c = self.c[i]
        return (((c[0] * dx + c[1]) * dx + c[2]) * dx + c[3]).tolist()


class _EOSVec:
    """Vectorised EOS (scipy splines, exactly as bsk26.EOS): eps(p), n(p), c_s^2(p) for arrays of p in MeV/fm^3.

    c_s^2 has three definitions (N4, step 23): 'eps' = dp/deps from the fits (the default of the whole project), 'n' = n dp/dn / (eps + p), and 'mid' = their geometric mean. The published fits
    eps(n) and P(rho) are independent and do not satisfy d eps/dn = (eps + p)/n exactly (see `mismatch`), so the two definitions differ by about 1e-3; the code uses 'eps'."""

    def __init__(self, n, e, p):
        lp = np.log(p)
        self._le, self._ln = CubicSpline(lp, np.log(e)), CubicSpline(lp, np.log(n))
        self.p_min = p[0]
        self.n, self.e, self.p = np.asarray(n), np.asarray(e), np.asarray(p)
        self.cs2_mode = 'eps'

    def eps(self, p):
        return np.exp(self._le(np.log(p)))

    def nb(self, p):
        return np.exp(self._ln(np.log(p)))

    def _cs2_eps(self, p):
        return p / self.eps(p) / self._le(np.log(p), 1)

    def _cs2_n(self, p):
        return p / (self.eps(p) + p) / self._ln(np.log(p), 1)

    def cs2(self, p):
        if self.cs2_mode == 'eps':
            return self._cs2_eps(p)
        if self.cs2_mode == 'n':
            return self._cs2_n(p)
        return np.sqrt(self._cs2_eps(p) * self._cs2_n(p))

    def mismatch(self, p):
        """(d ln eps / d ln n) eps / (eps + p) - 1: zero if the fits satisfy the first law d eps/dn = (eps + p)/n."""
        lp = np.log(p)
        e = self.eps(p)
        return float(self._le(lp, 1) / self._ln(lp, 1) * e / (e + p) - 1)


class StarTable:
    def __init__(self, eos_csv, star_csv):
        ec, _ = A.load_table(eos_csv)
        sc, sm = A.load_table(star_csv)
        par = sm['parameters']
        self.meta = sm
        self.eos = _EOS(ec['n_fm3'], ec['eps_MeV_fm3'], ec['p_MeV_fm3'])
        self.eosv = _EOSVec(ec['n_fm3'], ec['eps_MeV_fm3'], ec['p_MeV_fm3'])
        self.nu_a, self.p0, self.E0, self.nc = par['nu_a_kHz'], par['p0_MeV_fm3'], par['E0_geo_km2'], par['nc_fm3']
        self.ma = par['ma_per_km']
        self.fa2 = np.inf if par['fa2_geo'] is None else par['fa2_geo']
        self.beta, self.KGEO, self.thc = par['beta'], par['KGEO'], par['theta_c']
        self.MSUN = par['MSUN_km']
        self.R, self.M, self.Phi_shift = par['R_km'], par['M_km'], par['Phi_shift']
        self.rend, self.Mt = par['r_end_km'], par['M_total_km']
        self._thr = self.eos.p_min * 1.0001
        self._ax = bool(self.nu_a)
        r = sc['r_km']
        Y = np.column_stack([sc['ln_p'], sc['m_km'], sc['Phi'], sc['theta'], sc['thetap']])
        self.rk = par['r_k_km']
        i = int(np.searchsorted(r, self.rk))
        assert r[i] == self.rk, 'r_k must be a node of the table'
        # two splines, one on each side of r_k (S1 validates exactly this interpolation; one spline across the kink rings)
        # the tail (r > R_stop) has no fluid: no ln p column (stored as -inf), p = 0 there
        self._fi, self._fo = _FastSpline(r[:i + 1], Y[:i + 1]), _FastSpline(r[i:], Y[i:, 1:])
        self.r_nodes, self.i_stop, self.cols = r, i, sc          # raw nodes and columns (N3 builds its coefficient tables on them)

    # ---- background functions (same definitions as axion_ns/axstar.py)
    def g(self, t):
        return math.sqrt(1 - self.beta * math.sin(t / 2)**2)

    def g1(self, t):
        return -self.beta * math.sin(t) / (4 * self.g(t))

    def g2(self, t):
        gg = self.g(t)
        return -self.beta * math.cos(t) / (4 * gg) - self.beta**2 * math.sin(t)**2 / (16 * gg**3)

    def y_at(self, r):
        """(m, p, Phi, theta, theta') at r (Phi without the shift); replaces st.sol.sol(r)."""
        if r <= self.rk:
            lnp, m, Phi, t, tp = self._fi(r).tolist()
            return m, math.exp(lnp), Phi, t, tp
        m, Phi, t, tp = self._fo(r).tolist()
        return m, 0.0, Phi, t, tp

    def _matter(self, p, t):
        """(eps, n_b, x, c_s^2, q1, q2, g) in km^-2; zeros outside the star (same as Star.matter)."""
        pm = p / self.KGEO
        ax, b = self._ax, self.beta
        g = math.sqrt(1 - b * math.sin(t / 2)**2) if ax else 1.0
        if pm <= self._thr:
            return 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, g
        w = self.eos(math.log(pm))
        e0, nb, cs0 = math.exp(w[0]) * self.KGEO, math.exp(w[1]), pm / math.exp(w[0]) / w[2]
        x = nb / self.nc
        if ax:
            sn, cs = math.sin(t), math.cos(t)
            g1 = -b * sn / (4 * g)
            g2 = -b * cs / (4 * g) - b * b * sn * sn / (16 * g**3)
            e, q1, q2 = e0 - self.E0 * x * (1 - g), self.E0 * x * g1, self.E0 * x * g2
        else:
            e, q1, q2 = e0, 0.0, 0.0
        return e, nb, x, (e0 + p) * cs0 / (e + p), q1, q2, g

    def matter(self, p, t):
        """(eps, n_b, x, c_s^2 with the axion eps, q1 = n dm*/dtheta, q2 = n d2m*/dtheta2) in km^-2."""
        return self._matter(p, t)[:6]

    def _derived(self, r, p, m, t, tp):
        e, nb, x, cs2, q1, q2, g = self._matter(p, t)
        inside = nb != 0.0                      # nb > 0 exactly when p is above the EOS threshold
        F = 1 - 2 * m / r
        V = self.E0 * (1 - g) if self._ax else 0.0
        kin = self.fa2 * tp * tp / 2 * F if self._ax else 0.0
        pr = p if inside else 0.0
        dPhi = (m + 4 * math.pi * r**3 * (pr - V + kin)) / (r * r * F)
        dm = 4 * math.pi * r * r * (e + V + kin)
        return e, nb, x, cs2, q1, q2, dPhi, (dm * r - m) / (r * (r - 2 * m)), dm, -0.5 * math.log(F), V, pr, F, inside

    def rhs(self, r, y):
        """[dm, dp, dPhi, dtheta, dtheta'] (the background ODE right-hand side, as Star.rhs)."""
        m, p, Phi, t, tp = y
        e, nb, x, cs2, q1, q2, dPhi, dlam, dm, lam, V, pr, F, inside = self._derived(r, p, m, t, tp)
        dp = (-(e + pr) * dPhi - q1 * tp) if inside else 0.0
        if self._ax:
            src = 4 * self.ma**2 / self.beta * (x - 1) * self.g1(t)
            tpp = (src - 2 * tp / r * (1 - m / r - 2 * math.pi * r * r * (e - pr + 2 * V))) / F
        else:
            tpp = 0.0
        return [dm, dp, dPhi, tp, tpp]

    def bg(self, r):
        """(eps, p, n, c_s^2, q1, q2, Phi', lambda', Phi, lambda, theta, theta') at r (inside the star)."""
        m, p, Phi, t, tp = self.y_at(r)
        d = self._derived(r, p, m, t, tp)
        return d[0], p, d[1], d[3], d[4], d[5], d[6], d[7], Phi + self.Phi_shift, d[9], t, tp
