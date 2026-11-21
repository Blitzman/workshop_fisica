"""SLOW, manual:  python -m tests.verify_vs_legacy_n1     (about 5 min on one core)

Recomputes nu_a = 6 kHz with the ORIGINAL code (axion_ns: in-memory Star, dense output) and compares the modes with
N1's table-based result at full precision. This quantifies what the CSV + Hermite hand-over costs in accuracy.
Needs N1/outputs/modes_cowling.csv (witness) in place. Imports the legacy package on purpose (cross-module test).
"""
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
sys.path[:0] = [str(ROOT), str(ROOT.parent / 'archive')]

from axion_ns import axstar as OA                     # noqa: E402
from axion_ns import step9_axion_cowling as O9        # noqa: E402
from common import artifacts as A                     # noqa: E402


def main():
    cols, _ = A.load_table(ROOT / 'N1_modes_cowling' / 'outputs' / 'modes_cowling.csv')
    sel = cols['nu_a_kHz'] == 6.0
    st = OA.Star(6.0)
    bg = O9._background_modes(st)
    cp = [O9._coupled(st, g) for g in bg]
    print(f'{"mode":4s} {"bg old":>14s} {"bg new":>14s} {"rel diff":>9s} | {"cp old":>26s} {"cp new":>26s} {"rel diff":>9s}')
    worst = 0.0
    for k, nm in enumerate(('f', 'p1', 'p2')):
        new_bg = cols['nu_bg_kHz'][sel][k]
        new_cp = complex(cols['nu_cp_re_kHz'][sel][k], cols['nu_cp_im_kHz'][sel][k])
        d1, d2 = abs(new_bg / bg[k] - 1), abs(new_cp.real / cp[k].real - 1)
        worst = max(worst, d1, d2)
        print(f'{nm:4s} {bg[k]:14.9f} {new_bg:14.9f} {d1:9.1e} | {cp[k]:26.9f} {new_cp:26.9f} {d2:9.1e}')
    print(f'worst relative difference in Re nu: {worst:.2e}')
    print('Im nu (kHz): old', [f'{z.imag:.6e}' for z in cp], 'new', [f'{x:.6e}' for x in cols['nu_cp_im_kHz'][sel]])
    return worst


if __name__ == '__main__':
    main()
