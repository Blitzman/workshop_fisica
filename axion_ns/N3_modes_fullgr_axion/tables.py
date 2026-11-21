"""N3: Laurent-coefficient tables of the numerical matrices, built on the stored nodes of the star table.

A(r, w) = A_-1(r)/w^2 + A_0(r) + w^2 A_1(r) with real coefficient matrices that depend only on the star. They are evaluated EXACTLY at the nodes of the S1 table
(vectorised, three values of w^2 and a 3 x 3 solve) and interpolated with a cubic spline in r. What is tabulated is r^2 A (smooth at the regular singular
centre, where A has 1/r and 1/r^2 entries). One right-hand-side call then costs one spline lookup of the 3 n^2 columns plus a small complex combination
(about 10 us instead of about 520 us for the legacy bgvals + A). The nodes of the S1 table are uniform in ln p (dense in r near the surface), uniform in r
and geometric near the centre, so the coefficient tables inherit that resolution.
"""
import numpy as np

from .background import bg18
from common.star_table import _FastSpline

W_FIT = np.array([0.04, 0.12, 0.35])             # km^-1: three real frequencies for the Laurent fit
VINV = np.linalg.inv(np.column_stack([1 / W_FIT**2, np.ones(3), W_FIT**2]))


def _entries(flat, n, st, r, state, w):
    B = bg18(st, r, *state)
    vals = flat(r, w, 2, *B)
    out = np.empty((len(r), n * n))
    for j, v in enumerate(vals):
        out[:, j] = v
    return out


def _laurent(flat, n, st, r, state):
    A3 = np.stack([_entries(flat, n, st, r, state, w) for w in W_FIT])
    c = np.einsum('ij,jnk->ink', VINV, A3)
    return np.concatenate([c[0], c[1], c[2]], axis=1) * (r**2)[:, None]       # r^2 x [A_-1, A_0, A_1]


class LaurentTables:
    def __init__(self, st, flat_i, flat_e):
        self.st = st
        c, i = st.cols, st.i_stop
        r = c['r_km']
        ri, ro = r[:i + 1], r[i:]
        state_i = (c['m_km'][:i + 1], np.exp(c['ln_p'][:i + 1]), c['Phi'][:i + 1], c['theta'][:i + 1], c['thetap'][:i + 1])
        state_o = (c['m_km'][i:], np.zeros(len(ro)), c['Phi'][i:], c['theta'][i:], c['thetap'][i:])
        self._Yi = _laurent(flat_i, 6, st, ri, state_i)
        self._Yo = _laurent(flat_e, 4, st, ro, state_o)
        self.ri, self.ro = ri, ro
        self._fs = {}

    def _spline(self, region, n, keep):
        key = (region, keep)
        if key not in self._fs:
            Y, r = (self._Yi, self.ri) if region == 'i' else (self._Yo, self.ro)
            kk = range(n) if keep is None else keep
            cols = [k * n * n + a * n + b for k in range(3) for a in kk for b in kk]
            self._fs[key] = _FastSpline(r, Y[:, cols])
        return self._fs[key]

    def _rhs(self, region, n, omega, keep):
        fs = self._spline(region, n, keep)
        m = n if keep is None else len(keep)
        w2 = omega**2
        iw2, q = 1 / w2, m * m

        def f(rr, y):
            c = fs(rr)
            return ((c[:q] * iw2 + c[q:2 * q] + c[2 * q:] * w2).reshape(m, m) @ y) / (rr * rr)
        return f

    def rhs_int(self, omega, keep=None):
        return self._rhs('i', 6, omega, keep)

    def rhs_ext(self, omega, keep=None):
        return self._rhs('o', 4, omega, keep)

    def mat_ext(self, rr, omega):
        c = self._spline('o', 4, None)(rr)
        return (c[:16] / omega**2 + c[16:32] + c[32:] * omega**2).reshape(4, 4) / (rr * rr)
