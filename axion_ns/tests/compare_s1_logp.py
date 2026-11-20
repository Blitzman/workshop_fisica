"""SLOW, manual:  python -m tests.compare_s1_logp [--quick] [--stars 6,2.8] [--tight-s]

S1 only. Compares the structure integrated in r (existing Star) with the structure integrated in s = ln p (StarLogP: start at r0 = 10 m from the
series about the centre, stop at p = 1e-10 p_c, then the vacuum axion tail in r). Reference: the r-method at the original precision
(rtol 1e-11, theta_c to 1e-16). Definitions common to both methods: R10 = radius where p = 1e-10 p_c, M10 = m(R10), M40 = m(r = 40 km), profiles
at 40 values of s = ln(p/p_c) between -1e-4 and -23.
"""
import argparse
import sys
import time
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

ROOT = Path(__file__).resolve().parent.parent
sys.path[:0] = [str(ROOT), str(ROOT.parent)]

from S1_star_background.star import KGEO, Star, _eos       # noqa: E402
from S1_star_background.star_logp import StarLogP, series_start   # noqa: E402

PFRAC = 1e-10
SGRID = -np.geomspace(1e-4, 23.0, 40)
T = time.perf_counter


class StarR(Star):
    """The existing r-method with an rhs counter; r0 and the start series selectable (r0 = 1e-5 with the legacy series is the original)."""

    def __init__(self, nu_a=None, r0=1e-5, improved=False, **kw):
        self.nfev = self.nint = 0
        self.r0, self.improved = r0, improved
        super().__init__(nu_a, **kw)

    def rhs(self, r, y):
        self.nfev += 1
        return super().rhs(r, y)

    def _start(self, r0, thc):
        return series_start(self, r0, thc) if self.improved else super()._start(r0, thc)

    def _integrate(self, thc, keep=False, rmax=None):          # Star._integrate with r0 = self.r0
        self.nint += 1
        r0 = self.r0
        rmax = rmax or (60.0 + (40 / self.ma if self.nu_a else 0.0))
        up = lambda r, y: (y[3] - np.pi - 1e-3) if self.nu_a else 1.0
        down = lambda r, y: y[3] + 1e-3 if self.nu_a else 1.0
        turn = lambda r, y: y[4] if (self.nu_a and y[1] / KGEO < _eos.p_min * 1.0001) else -1.0
        for ev in (up, down):
            ev.terminal = True
        turn.terminal, turn.direction = True, 1
        s = solve_ivp(self.rhs, (r0, rmax), self._start(r0, thc), events=[up, down, turn], method='DOP853', rtol=self.rtol, atol=self.atol,
                      dense_output=keep, max_step=self.max_step)
        if s.status < 0 and not s.y[1, -1] / KGEO <= _eos.p_min * 1.0002:
            raise RuntimeError(f'StarR: {s.message}')
        if keep:
            self.sol = s
        if not self.nu_a:
            return 0
        return -1 if s.t_events[1].size else +1


def summarize(st, kind):
    """thc, R10, M10, M40, shift and the profiles at SGRID, with the common definitions."""
    pc = st.p0 * KGEO
    out = {'thc': st.thc, 'shift': st.Phi_shift}
    if kind == 'r':
        f = lambda r, s: np.log(st.sol.sol(r)[1] / pc) - s
        rr = [brentq(f, st.sol.t[0], st.R, args=(s,), xtol=1e-13) for s in SGRID]
        R10 = brentq(f, st.sol.t[0], st.R, args=(np.log(PFRAC),), xtol=1e-13)
        prof = np.array([st.sol.sol(r) for r in rr])                   # (n, 5): m, p, Phi, theta, theta'
        out.update(R10=R10, M10=st.sol.sol(R10)[0], M40=st.sol.sol(40.0)[0] if st.sol.t[-1] >= 40 else np.nan,
                   prof=dict(r=np.array(rr), m=prof[:, 0], Phi=prof[:, 2] + st.Phi_shift, th=prof[:, 3], thp=prof[:, 4]))
    else:
        prof = np.array([st.at_s(s) for s in SGRID])                   # r, m, Phi+shift, theta, theta'
        M40 = (st.at_r_tail(40.0)[0] if (st.nu_a and st.sol_out.t[-1] >= 40) else (st.M_stop if not st.nu_a else np.nan))
        out.update(R10=st.R_stop, M10=st.M_stop, M40=M40, prof=dict(r=prof[:, 0], m=prof[:, 1], Phi=prof[:, 2], th=prof[:, 3], thp=prof[:, 4]))
    return out


def rel(a, b, floor=0.0):
    sc = np.abs(b).max()
    return float(np.max(np.abs(a - b) / np.maximum(np.abs(b), max(floor * sc, 1e-300))))


def profile_err(x, ref):
    return {k: rel(x['prof'][k], ref['prof'][k], 1e-6 if k in ('th', 'thp') else 0.0) for k in ref['prof']}


def build(label, kind, nu, rtol, shoot):
    """Returns (star, summary, timing dict)."""
    kw = dict(rtol=rtol, shoot_tol=shoot)
    t0 = T()
    if kind == 'r':
        st = StarR(nu, r0=1e-5, improved=False, **kw)
    elif kind == 'r10':
        st = StarR(nu, r0=0.01, improved=True, **kw)
    else:
        st = StarLogP(nu, **kw)
    tb = T() - t0
    nint, nf = st.nint, st.nfev
    # one integration at the final theta_c (after the build, so the dense solution exists)
    n0 = st.nfev
    t0 = T(); st._integrate(st.thc, keep=False); t1 = T() - t0
    return st, summarize(st, 'r' if kind in ('r', 'r10') else 's'), dict(build=tb, nint=nint, rhs_per=(nf / max(nint, 1)), one=t1, one_rhs=st.nfev - n0)


def start_check():
    print('A) series start at r0 = 10 m against the integration from 1e-5 km (rtol 1e-12), relative differences of the state at 10 m', flush=True)
    print(f'   {"star":10s} {"m":>9s} {"p":>9s} {"Phi_prime":>10s} {"theta":>9s} {"theta_p":>9s}   (Phi(r0) = 0 is a convention absorbed by the final shift; its derivative is compared instead)', flush=True)
    for label, nu in (('no axion', None), ('6 kHz', 6.0), ('2.8 kHz', 2.8), ('15 kHz', 15.0)):
        st = Star(nu)
        r0 = 0.01
        y0 = Star._start(st, 1e-5, st.thc)
        s = solve_ivp(st.rhs, (1e-5, r0), y0, method='DOP853', rtol=1e-13, atol=1e-40)
        a = np.array(s.y[:, -1])
        b = np.array(series_start(st, r0, st.thc))
        e_c, _, _, _, _, _ = st.matter(st.p0 * KGEO, st.thc)
        V_c = st.E0 * (1 - __import__('S1_star_background.star', fromlist=['g']).g(st.thc)) if st.nu_a else 0.0
        dPhi_series = (4 * np.pi * r0 / 3) * (e_c + 3 * st.p0 * KGEO - 2 * V_c)
        d = [abs(b[0] / a[0] - 1), abs(b[1] / a[1] - 1), abs(dPhi_series / st.rhs(r0, a)[2] - 1)]
        d += [abs(b[3] - a[3]) / max(abs(a[3]), 1e-300), abs(b[4] - a[4]) / max(abs(a[4]), 1e-300) if a[4] != 0 else 0.0]
        print(f'   {label:10s} ' + ' '.join(f'{v:9.1e}' if i != 2 else f'{v:10.1e}' for i, v in enumerate(d)), flush=True)


def report(title, ref, rows):
    print(f'\n{title}', flush=True)
    th = ref['thc']
    print(f'   reference (r-method, rtol 1e-11, original precision): theta_c = {th:.15f}  delta = {np.pi - th:.6e}  R10 = {ref["R10"]:.8f} km  '
          f'M10 = {ref["M10"]:.10f} km  M40 = {ref["M40"]:.10f} km  Phi shift = {ref["shift"]:.8f}', flush=True)
    print(f'   {"method":26s} {"|dtheta_c|/delta":>16s} {"R10":>9s} {"M10":>9s} {"M40":>9s} | profile max rel err: {"r":>8s} {"m":>8s} {"Phi":>8s} {"theta":>8s} {"theta_p":>8s} | '
          f'{"build (s)":>9s} {"integr.":>7s} {"rhs/int":>8s} {"1 int (s)":>9s}', flush=True)
    for label, x, tm in rows:
        pe = profile_err(x, ref)
        dth = abs(x['thc'] - th) / max(np.pi - th, 1e-300) if ref['thc'] else 0.0
        print(f'   {label:26s} {dth:16.1e} {abs(x["R10"] / ref["R10"] - 1):9.1e} {abs(x["M10"] / ref["M10"] - 1):9.1e} '
              f'{abs(x["M40"] / ref["M40"] - 1) if np.isfinite(x["M40"]) and np.isfinite(ref["M40"]) else float("nan"):9.1e} | {"":22s}'
              + ' '.join(f'{pe[k]:8.1e}' for k in ('r', 'm', 'Phi', 'th', 'thp'))
              + f' | {tm["build"]:9.1f} {tm["nint"]:7d} {tm["rhs_per"]:8.0f} {tm["one"]:9.3f}', flush=True)


def case(nu, tight_s):
    name = 'no axion' if nu is None else f'nu_a = {nu:g} kHz'
    print(f'\n######## {name}', flush=True)
    st, ref, tm = build('ref', 'r', nu, 1e-11, None)
    print(f'   reference built in {tm["build"]:.1f} s ({tm["nint"]} integrations)', flush=True)
    rows = [('r  rtol 1e-11 (reference)', ref, tm)]
    plan = [('r   rtol 1e-6 / shoot 1e-6', 'r', 1e-6, 1e-6), ('r10 rtol 1e-6 (10 m start)', 'r10', 1e-6, 1e-6), ('ln p rtol 1e-6', 's', 1e-6, 1e-6),
            ('r   rtol 1e-8 / shoot 1e-8', 'r', 1e-8, 1e-8), ('r10 rtol 1e-8 (10 m start)', 'r10', 1e-8, 1e-8), ('ln p rtol 1e-8', 's', 1e-8, 1e-8)]
    if tight_s:
        plan.append(('ln p rtol 1e-11 (original)', 's', 1e-11, None))
    for label, kind, rt, sh in plan:
        if nu is None:
            sh = None
        st2, x, tm2 = build(label, kind, nu, rt, sh)
        rows.append((label, x, tm2))
    report(f'B) {name}: structure integrated in r vs in ln p', ref, rows)


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--quick', action='store_true')
    ap.add_argument('--stars', default='6,2.8')
    ap.add_argument('--tight-s', action='store_true')
    a = ap.parse_args()
    start_check()
    case(None, True)
    if not a.quick:
        for nu in [float(x) for x in a.stars.split(',') if x]:
            case(nu, a.tight_s)
