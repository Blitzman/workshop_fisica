"""Experimental axion shooting with a finite-radius exterior match.

From the consolidated package root:
  PYTHONPATH=axion_ns python /path/to/axion_brent.py 3 6 12 24 48
  PYTHONPATH=axion_ns python /path/to/axion_brent.py 3 6 12 24 48 --corrector bisect

Frequency arguments are in kHz. The default match is five reduced Compton
wavelengths beyond the surface. Beyond the match, the logarithmic derivative
uses the linear scalar equation on Schwarzschild with the enclosed mass at
the match. The profile ends there. The corrector uses xtol=1e-9 in u_c.
This is an independent driver for the positive condensate branch, using
the consolidated S1 EOS, hydrostatic RHS, and physical conventions.
"""
import math
from time import perf_counter
import numpy as np
from S1_star_background import star_logp as sl
from S1_star_background.star import BETA, MSUN

original_solve = sl.solve_ivp
class LogOffset(sl.StarLogP):
    def _shoot(self):
        lo, hi = np.log(1e-300), np.log(np.pi-0.3)
        def trial(u, keep=False):
            self.u_c = u
            return self._integrate(np.pi-np.exp(u),keep=keep)
        if trial(hi) > 0:
            raise ValueError('sin condensado')
        for _ in range(200):
            mid = (lo+hi)/2
            if trial(mid) > 0:
                lo = mid
            else:
                hi = mid
            if hi-lo < self.shoot_tol:
                break
        self.thc = np.pi-np.exp(lo)
        trial(lo,keep=True)

    def _rhs_logoffset(self, s, z):
        r,m,Phi,u,w = z
        d = math.exp(u)
        t = math.pi-d
        tp = -d*w
        p = self.p0 * sl.KGEO * math.exp(s)
        dm, dp, dPhi, e, x, F, _, _, _ = self._hydro_rhs(r, m, p, t, tp)
        gg = math.sqrt(1-BETA*math.cos(d/2)**2)
        V = self.E0*(1-gg)
        c = 2/r*(1-m/r-2*math.pi*r*r*(e-p+2*V))/F
        sinc = math.sin(d)/d if d else 1.
        wp = -w*w-c*w-self.ma**2*(1-x)*sinc/(gg*F)
        self.nfev += 1
        drds = p/dp
        return [drds,dm*drds,dPhi*drds,w*drds,wp*drds]

def physical(z):
    y = np.array(z)
    d = np.exp(y[3])
    y[3] = np.pi-d
    y[4] *= -d
    return y

def dispatch(fun, span, y0, **kw):
    star = fun.__self__
    if fun.__name__ != '_rhs_s' or not isinstance(star,LogOffset):
        return original_solve(fun,span,y0,**kw)
    r,m,Phi,t,tp = y0
    d = math.exp(star.u_c)
    x = star._background_matter(star.p0 * sl.KGEO, 0.0, 0.0)[3]
    gg = math.sqrt(1-BETA*math.cos(d/2)**2)
    w = -r*star.ma**2*(1-x)*(math.sin(d)/d)/(3*gg)
    u = star.u_c + w*r/2
    def event(e):
        def transformed(s,z):
            return e(s,physical(z))
        transformed.terminal=e.terminal
        return transformed
    kw['events']=[event(e) for e in kw['events']]
    sol=original_solve(star._rhs_logoffset,span,[r,m,Phi,u,w],**kw)
    sol.y=physical(sol.y)
    if sol.sol is not None:
        dense=sol.sol
        sol.sol=lambda s: physical(dense(s))
    return sol


# The interior uses the existing log-offset transformation. The exterior
# matching experiment below has its own shooting and stopping criteria.
from argparse import ArgumentParser
from scipy.optimize import bisect, brentq


class MatchedStar(LogOffset):
    def __init__(self, nu, match, seed=None, corrector=brentq):
        self.match = match
        self.seed = seed
        self.corrector = corrector
        self.trials = 0
        super().__init__(nu, shoot_tol=1e-9)

    def trial(self, u, keep=False):
        if not keep:
            self.trials += 1
        self.u_c = u
        thc = np.pi - np.exp(u)
        m, p, phi, theta, tp = sl.series_start(self, self.r0, thc)
        pc = self.p0 * sl.KGEO

        def interior_cross(s, z):
            return z[3] + 1e-3

        interior_cross.terminal = True
        interior = dispatch(self._rhs_s, (math.log(p/pc), math.log(self.p_frac)),
                            [self.r0,m,phi,theta,tp], events=[interior_cross],
                            method='DOP853', rtol=self.rtol, atol=self.atol,
                            max_step=self.max_step_s, dense_output=keep)
        if not interior.success:
            raise RuntimeError(interior.message)
        if interior.status == 1:
            return -1, None
        R,m,phi,theta,tp = interior.y[:, -1]
        rm = R + self.match / self.ma

        def cross(r, y):
            return y[3] + 1e-3

        def turn(r, y):
            return y[4]

        cross.terminal = turn.terminal = True
        cross.direction, turn.direction = -1, 1
        # A trial already heading away from the decaying branch at R
        # supplies a bracket sign without integrating a long bad tail.
        if tp >= 0 and theta > 0:
            return 1, None
        exterior = original_solve(self._rhs_r, (R,rm), [m,0,phi,theta,tp],
                                  events=[cross,turn], method='DOP853',
                                  rtol=self.rtol, atol=self.atol,
                                  max_step=self.max_step, dense_output=keep)
        if not exterior.success:
            raise RuntimeError(exterior.message)
        if exterior.status == 1:
            return (-1 if exterior.t_events[0].size else 1), None
        mass,_,_,theta,tp = exterior.y[:, -1]
        r_far = rm + 30/self.ma

        def riccati(r, s):
            F = 1 - 2*mass/r
            return self.ma**2/F - 2*(1-mass/r)*s/(r*F) - s*s

        decay = original_solve(riccati, (r_far,rm), [-self.ma-(1+mass*self.ma)/r_far],
                               method='DOP853', rtol=1e-10, atol=1e-12)
        if not decay.success:
            raise RuntimeError(decay.message)
        slope = decay.y[0,-1]
        residual = tp - slope*theta
        if keep:
            self.thc, self.sol_in, self.sol_out = thc,interior,exterior
            self.R_stop, self.M_stop = R,interior.y[1,-1]
            self.theta_match = theta
            self.decay_mismatch = residual/(self.ma*abs(theta))
        return np.sign(residual), residual

    def _shoot(self):
        cache = {}
        def sample(u):
            if u not in cache:
                cache[u] = self.trial(u)
            return cache[u]
        broad = (math.log(1e-300), math.log(np.pi-0.3))
        lo,hi = broad
        if self.seed is not None:
            for width in (0.25, 1, 4, 16):
                left,right = max(lo,self.seed-width), min(hi,self.seed+width)
                if left < right and sample(left)[0] > 0 > sample(right)[0]:
                    lo,hi=left,right
                    break
        if not sample(lo)[0] > 0 > sample(hi)[0]:
            raise ValueError('Could not bracket the condensate')
        # Early events supply signs only. Brent starts once both
        # endpoints actually reach the matching radius.
        for _ in range(64):
            if sample(lo)[1] is not None and sample(hi)[1] is not None:
                break
            mid=(lo+hi)/2
            if sample(mid)[0] > 0:
                lo=mid
            else:
                hi=mid
        else:
            raise RuntimeError('Could not bracket a finite matching residual')
        self.bracket_trials=self.trials

        def residual(u):
            value=sample(u)[1]
            if value is None:
                raise RuntimeError('Brent trial stopped before the matching radius')
            return value

        root = self.corrector(residual, lo, hi, xtol=self.shoot_tol)
        if self.trial(root, keep=True)[1] is None:
            raise RuntimeError("Final trial stopped before the matching radius")


def main():
    parser=ArgumentParser()
    parser.add_argument('frequency', type=float, nargs='+', help='Axion frequency in kHz')
    parser.add_argument('--match', type=float, default=5, help='Matching distance in reduced Compton wavelengths')
    parser.add_argument('--corrector', choices=('brent', 'bisect'), default='brent')
    args=parser.parse_args()
    if args.match <= 0 or any(nu <= 0 for nu in args.frequency):
        parser.error('Matching distance and frequencies must be positive')
    corrector = brentq if args.corrector == 'brent' else bisect
    history=[]
    for nu in args.frequency:
        seed=None
        if history:
            seed=history[-1][1]
        if len(history)>=2:
            (n1,u1),(n2,u2)=history[-2:]
            if n1!=n2:
                seed=u2+(u2-u1)*(nu-n2)/(n2-n1)
        start=perf_counter()
        star=MatchedStar(nu,args.match,seed,corrector)
        history.append((nu,star.u_c))
        print(dict(nu=nu,time=perf_counter()-start,bracket=star.bracket_trials,
                   corrector_trials=star.trials-star.bracket_trials,u=star.u_c,R=star.R,
                   M=star.M/MSUN,rm=star.r_end,theta_match=star.theta_match,
                   relative_mismatch=star.decay_mismatch),flush=True)

if __name__ == '__main__':
    main()
