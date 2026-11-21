"""N4: coefficient tables of the radial matrices, built on the stored nodes of the star table.

For the radial system A(r, w) = A_0(r) + w^2 A_1(r) exactly (two real coefficient matrices; N3's system also had a 1/w^2 term). They are evaluated EXACTLY at the nodes of S1's
table (vectorised, two values of w and a 2 x 2 solve; a third value checks the structure) and interpolated with a cubic spline in r. What is tabulated is r^2 A (smooth at the regular
singular centre). One right-hand-side call costs one spline lookup of the 2 n^2 columns plus a small complex combination (about 7 us instead of about 200 us for bgvals + A). The
interior entries grow like 1e13 towards the surface (as the polytrope's W in N2); the interpolation of r^2 A keeps the product A y accurate to about 1e-8 relative (measured).
"""
import numpy as np

from .background import bg18
from common.star_table import _FastSpline

W_FIT = (0.05, 0.30)             # km^-1: two real frequencies for the fit; W_CHECK verifies A_0 + w^2 A_1 at a third
W_CHECK = 0.12


def _entries(flat, n, st, r, state, w):
    vals = flat(r, w, *bg18(st, r, *state))
    out = np.empty((len(r), n * n))
    for j, v in enumerate(vals):
        out[:, j] = v
    return out


def _coefficients(flat, n, st, r, state):
    a1_, a2_ = (_entries(flat, n, st, r, state, w) for w in W_FIT)
    A1 = (a2_ - a1_) / (W_FIT[1]**2 - W_FIT[0]**2)
    A0 = a1_ - W_FIT[0]**2 * A1
    chk = _entries(flat, n, st, r, state, W_CHECK)
    err = np.abs(A0 + W_CHECK**2 * A1 - chk).max() / max(np.abs(chk).max(), 1e-300)
    if err > 1e-9:
        raise ValueError(f'the matrix is not A_0 + w^2 A_1 (check at w = {W_CHECK}: {err:.1e})')
    return np.concatenate([A0, A1], axis=1) * (r**2)[:, None], err


class RadialTables:
    def __init__(self, st, flat_i, flat_e):
        self.st, self.flat_e = st, flat_e
        c, i = st.cols, st.i_stop
        r = c['r_km']
        ri, ro = r[:i + 1], r[i:]
        state_i = (c['m_km'][:i + 1], np.exp(c['ln_p'][:i + 1]), c['Phi'][:i + 1], c['theta'][:i + 1], c['thetap'][:i + 1])
        state_o = (c['m_km'][i:], np.zeros(len(ro)), c['Phi'][i:], c['theta'][i:], c['thetap'][i:])
        self._Yi, ei = _coefficients(flat_i, 4, st, ri, state_i)
        self._Yo, eo = _coefficients(flat_e, 2, st, ro, state_o)
        self.structure_error = max(ei, eo)
        self.ri, self.ro = ri, ro
        self._fs = {}

    def _spline(self, region, n, keep):
        key = (region, keep)
        if key not in self._fs:
            Y, r = (self._Yi, self.ri) if region == 'i' else (self._Yo, self.ro)
            kk = range(n) if keep is None else keep
            cols = [k * n * n + a * n + b for k in range(2) for a in kk for b in kk]
            self._fs[key] = _FastSpline(r, Y[:, cols])
        return self._fs[key]

    def _rhs(self, region, n, omega, keep):
        fs = self._spline(region, n, keep)
        m = n if keep is None else len(keep)
        w2, q = omega**2, m * m

        def f(rr, y):
            c = fs(rr)
            return ((c[:q] + c[q:] * w2).reshape(m, m) @ y) / (rr * rr)
        return f

    def rhs_int(self, omega, keep=None):
        return self._rhs('i', 4, omega, keep)

    def rhs_ext(self, omega, keep=None):
        return self._rhs('o', 2, omega, keep)

    def mat_int(self, rr, omega, keep=None):
        fs = self._spline('i', 4, keep)
        m = 4 if keep is None else len(keep)
        c = fs(rr)
        q = m * m
        return (c[:q] + c[q:] * omega**2).reshape(m, m) / (rr * rr)
