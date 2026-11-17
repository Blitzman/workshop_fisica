"""N1: the background star as seen by the mode solver, rebuilt from S1's tables (no import from S1).

Inputs (written by S1): eos_bsk26.csv and star_*.csv with their .meta.json. The EOS splines are rebuilt exactly as
bsk26.EOS builds them from the same data; y(r) = (m, p, Phi, theta, theta') is reconstructed with cubic Hermite
interpolation from the stored values and derivatives (S1 guarantees the interpolation error, see its metadata).
`matter` and `rhs` are copies of the two short functions of axion_ns/axstar.py::Star (the TOV + Klein-Gordon
right-hand side), needed for dPhi, dm at arbitrary r; tests/ checks them against S1.
"""
import numpy as np
from scipy.interpolate import CubicHermiteSpline, CubicSpline

from common import artifacts as A


class EOSTable:
    def __init__(self, n, e, p):
        self.n, self.e, self.p = n, e, p
        lp = np.log(p)
        self._le = CubicSpline(lp, np.log(e))
        self._ln = CubicSpline(lp, np.log(n))
        self.p_min = p[0]

    def eps(self, p):
        return np.exp(self._le(np.log(p)))

    def nb(self, p):
        return np.exp(self._ln(np.log(p)))

    def cs2(self, p):
        return p / self.eps(p) / self._le(np.log(p), 1)


class StarTable:
    def __init__(self, eos_csv, star_csv):
        ec, em = A.load_table(eos_csv)
        sc, sm = A.load_table(star_csv)
        self.eos = EOSTable(ec['n_fm3'], ec['eps_MeV_fm3'], ec['p_MeV_fm3'])
        par = sm['parameters']
        self.meta = sm
        self.nu_a, self.p0, self.E0, self.nc = par['nu_a_kHz'], par['p0_MeV_fm3'], par['E0_geo_km2'], par['nc_fm3']
        self.ma = par['ma_per_km']
        self.fa2 = np.inf if par['fa2_geo'] is None else par['fa2_geo']
        self.beta, self.KGEO, self.thc = par['beta'], par['KGEO'], par['theta_c']
        self.MSUN = par['MSUN_km']
        self.R, self.M, self.Phi_shift = par['R_km'], par['M_km'], par['Phi_shift']
        self.rend, self.Mt = par['r_end_km'], par['M_total_km']
        r = sc['r_km']
        Y = np.column_stack([sc['m_km'], sc['p_km2'], sc['Phi'], sc['theta'], sc['thetap']])
        dY = np.column_stack([sc['dm'], sc['dp'], sc['dPhi'], sc['dtheta'], sc['dthetap']])
        self._H = CubicHermiteSpline(r, Y, dY)

    # ---- background functions (same definitions as axion_ns/axstar.py)
    def g(self, t):
        return np.sqrt(1 - self.beta * np.sin(t / 2)**2)

    def g1(self, t):
        return -self.beta * np.sin(t) / (4 * self.g(t))

    def g2(self, t):
        gg = self.g(t)
        return -self.beta * np.cos(t) / (4 * gg) - self.beta**2 * np.sin(t)**2 / (16 * gg**3)

    def y_at(self, r):
        """(m, p, Phi, theta, theta') at r (Phi without the shift), replaces st.sol.sol(r)."""
        return self._H(r)

    def matter(self, p, t):
        """(eps, n_b, x, c_s^2 con el eps axiónico, q1 = n dm*/dtheta, q2 = n d2m*/dtheta2) en km^-2."""
        pm = p / self.KGEO
        if pm <= self.eos.p_min * 1.0001:
            return 0.0, 0.0, 0.0, 1.0, 0.0, 0.0
        e0, nb, cs0 = self.eos.eps(pm) * self.KGEO, self.eos.nb(pm), self.eos.cs2(pm)
        x = nb / self.nc
        e = e0 - self.E0 * x * (1 - self.g(t)) if self.nu_a else e0
        cs2 = (e0 + p) * cs0 / (e + p)                # dp/dn n / (eps + p), con eps(n, theta)
        q1 = self.E0 * x * self.g1(t) if self.nu_a else 0.0
        q2 = self.E0 * x * self.g2(t) if self.nu_a else 0.0
        return e, nb, x, cs2, q1, q2

    def rhs(self, r, y):
        m, p, Phi, t, tp = y
        inside = p / self.KGEO > self.eos.p_min * 1.0001
        e, nb, x, cs2, q1, q2 = self.matter(p, t) if inside else (0, 0, 0, 1, 0, 0)
        pr = p if inside else 0.0
        F = 1 - 2 * m / r
        V = self.E0 * (1 - self.g(t)) if self.nu_a else 0.0
        kin = self.fa2 * tp**2 / 2 * F if self.nu_a else 0.0
        dPhi = (m + 4 * np.pi * r**3 * (pr - V + kin)) / (r**2 * F)
        dm = 4 * np.pi * r**2 * (e + V + kin)
        dp = (-(e + pr) * dPhi - q1 * tp) if inside else 0.0
        if self.nu_a:
            src = 4 * self.ma**2 / self.beta * (x - 1) * self.g1(t)
            tpp = (src - 2 * tp / r * (1 - m / r - 2 * np.pi * r**2 * (e - pr + 2 * V))) / F
        else:
            tpp = 0.0
        return [dm, dp, dPhi, tp, tpp]

    def bg(self, r):
        """(eps, p, n, c_s^2, q1, q2, Phi', lambda', Phi, lambda, theta, theta') at r (inside the star)."""
        m, p, Phi, t, tp = self.y_at(r)
        e, nb, x, cs2, q1, q2 = self.matter(p, t)
        F = 1 - 2 * m / r
        d = self.rhs(r, [m, p, Phi, t, tp])
        dPhi, dm = d[2], d[0]
        lam = -0.5 * np.log(F)
        dlam = (dm * r - m) / (r * (r - 2 * m))
        return e, p, nb, cs2, q1, q2, dPhi, dlam, Phi + self.Phi_shift, lam, t, tp
