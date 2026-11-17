"""SLOW, manual:  python -m tests.compare_n2_fullmode      (about 3 minutes on one core)

LMI f-mode (l = 2) of the 1.4 Msun polytrope through full_mode (minimum of |A_out|^2 + curvature), legacy code vs the modular N2 at several
integration tolerances. Reference: Lindblom-Mendell-Ipser 1997, Re w = 1.2016, Im w = 3.95e-4 (units of sqrt(M/R^3)).
"""
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
sys.path[:0] = [str(ROOT), str(ROOT.parent)]

from axion_ns import fmodes as FM                      # noqa: E402
from N2_modes_fullgr_noaxion import modes as M2        # noqa: E402
from N2_modes_fullgr_noaxion import systems as S2      # noqa: E402
from N2_modes_fullgr_noaxion.polytrope import Star     # noqa: E402

G, C, MSUN = 6.6743e-8, 2.99792458e10, 1.98841e33
MASS = 1.4 * G * MSUN / C**2 / np.sqrt(1.3346e5 / G)
LMI = (1.2016, 3.95e-4)
T = time.perf_counter


def main():
    ext, _, _ = S2.load_inputs(ROOT / 'D2_derive_nonradial' / 'outputs', ROOT / 'R_paper_refs' / 'outputs')
    full = S2.full_system(ext, 'energy')
    Lfull = FM.full_system('energy')
    st, sl = Star('energy', mass=MASS), FM.Star('energy', mass=MASS)
    print(f'LMI reference: Re w = {LMI[0]}, Im w = {LMI[1]}', flush=True)
    print(f'{"variant":22s} {"rtol":>6s} {"Re w":>12s} {"Re - LMI":>10s} {"Im w":>12s} {"Im/LMI-1":>9s} {"time (s)":>9s} {"rhs calls":>10s}', flush=True)
    t0 = T(); s0, g0 = FM.full_mode(sl, Lfull, 2, LMI[0], width=0.06); tl = T() - t0
    print(f'{"legacy (no cse)":22s} {"1e-11":>6s} {s0:12.9f} {s0 - LMI[0]:+10.1e} {g0:12.6e} {g0 / LMI[1] - 1:+9.1e} {tl:9.1f} {"-":>10s}', flush=True)
    res = {}
    for rt in (1e-11, 1e-10, 1e-9, 1e-8):
        st.nfev = 0
        t0 = T(); s, g = M2.full_mode(st, full, 2, LMI[0], width=0.06, rtol=rt); tt = T() - t0
        res[rt] = (s, g)
        print(f'{"N2 new (cse)":22s} {rt:6.0e} {s:12.9f} {s - LMI[0]:+10.1e} {g:12.6e} {g / LMI[1] - 1:+9.1e} {tt:9.1f} {st.nfev:10d}', flush=True)
    print('\ndifference from the LEGACY result:')
    for rt, (s, g) in res.items():
        print(f'  rtol {rt:6.0e}:  Re w {s - s0:+.2e} (relative {abs(s / s0 - 1):.1e})   Im w relative {abs(g / g0 - 1):.1e}', flush=True)


if __name__ == '__main__':
    main()
