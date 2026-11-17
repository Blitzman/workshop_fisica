"""SLOW, manual:  python -m tests.compare_n2_atol      (about 1.5 minutes on one core)

Is the tolerance sensitivity of the amplitude |A_out|^2 an absolute-tolerance effect? One amplitude evaluation of the LMI star (l = 2) at
two frequencies for several atol (same atol for the interior integrations and the Zerilli exterior), plus the component sizes and the
cancellation in the X(R) = 0 combination.
"""
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
sys.path[:0] = [str(ROOT), str(ROOT.parent)]

from N2_modes_fullgr_noaxion import modes as M2        # noqa: E402
from N2_modes_fullgr_noaxion import systems as S2      # noqa: E402
from N2_modes_fullgr_noaxion.polytrope import Star     # noqa: E402

G, C, MSUN = 6.6743e-8, 2.99792458e10, 1.98841e33
MASS = 1.4 * G * MSUN / C**2 / np.sqrt(1.3346e5 / G)
T = time.perf_counter


def main():
    ext, _, _ = S2.load_inputs(ROOT / 'D2_derive_nonradial' / 'outputs', ROOT / 'R_paper_refs' / 'outputs')
    full = S2.full_system(ext, 'energy')
    st = Star('energy', mass=MASS)
    u = np.sqrt(st.M / st.R**3)
    print('--- components at the end of the two interior integrations (rtol 1e-11, atol 1e-14)', flush=True)
    for s_ in (1.2016, 1.45):
        om = s_ * u
        y1 = st.integrate(full, st.center(om, 2, 1.0, 0.0), om, 2)[3:]
        y2 = st.integrate(full, st.center(om, 2, 0.0, 1.0), om, 2)[3:]
        c = -y1[3] / y2[3]
        y = y1 + c * y2
        print(f'omega = {s_} u:  y1 = {np.array2string(y1, precision=3)}   y2 = {np.array2string(y2, precision=3)}', flush=True)
        print(f'              c = -X1/X2 = {c:.3e};  y = y1 + c y2 = {np.array2string(y, precision=3)};  cancellation (max(|y1|, |c y2|) / |y|) per component: '
              f'{np.array2string(np.maximum(abs(y1), abs(c * y2)) / np.maximum(abs(y), 1e-300), precision=1)}', flush=True)
    print('\n--- amplitude |A_out|^2: atol scan (value relative to the default rtol 1e-11, atol 1e-14)', flush=True)
    print(f'{"omega":>8s} {"rtol":>6s} {"atol":>7s} {"value":>17s} {"rel diff":>9s} {"time (s)":>9s} {"rhs calls":>10s}', flush=True)
    for s_ in (1.2016, 1.45):
        om = s_ * u
        base = M2.outgoing_amplitude(st, full, om, 2)
        for rt, at in [(1e-11, 1e-12), (1e-11, 1e-14), (1e-11, 1e-16), (1e-11, 1e-18), (1e-11, 1e-20), (1e-11, 1e-24),
                       (1e-10, 1e-14), (1e-10, 1e-18), (1e-10, 1e-24)]:
            st.nfev = 0
            t0 = T(); v = M2.outgoing_amplitude(st, full, om, 2, rtol=rt, atol=at); tt = T() - t0
            print(f'{s_:7.4f}u {rt:6.0e} {at:7.0e} {v:17.9e} {abs(v / base - 1):9.1e} {tt:9.2f} {st.nfev:10d}', flush=True)


if __name__ == '__main__':
    main()
