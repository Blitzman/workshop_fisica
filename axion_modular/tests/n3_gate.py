"""N3 against the committed legacy results (results/fig_axgr.json, 6 significant digits): the gate for consolidating N3.

    python -m tests.n3_gate          prints the comparison for the 12 values of nu_a (needs N3/outputs/modes_axion.csv from a full run)

Criterion: every relative difference below 1e-3 (the tolerance set for this comparison): Re nu and damping of the coupled modes, Re nu of the
background-only modes (rebuilt from the stored percent shift and the no-axion reference), and the no-axion references.
"""
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
sys.path[:0] = [str(ROOT), str(ROOT.parent)]
from common import artifacts as A      # noqa: E402

LEGACY = ROOT.parent / 'results' / 'fig_axgr.json'
NEW = ROOT / 'N3_modes_fullgr_axion' / 'outputs' / 'modes_axion.csv'
TOL = 1e-3


def compare():
    old = json.loads(LEGACY.read_text(encoding='utf-8'))
    c, _ = A.load_table(NEW)
    rows = []
    for i, nm in enumerate(('f', 'p1', 'p2')):
        gr_re_old, gr_d_old = old['gr0'][i]
        s0 = (c['mode'] == nm)
        rows.append((f'no axion {nm}', 'Re', gr_re_old, float(c['gr0_re_kHz'][s0][0])))
        rows.append((f'no axion {nm}', 'damping Hz', gr_d_old, -1e3 * float(c['gr0_im_kHz'][s0][0])))
        m = old['modes'][i]
        for j, nu in enumerate(old['nus']):
            s = (c['mode'] == nm) & (c['nu_a_kHz'] == nu)
            if not s.any():
                continue
            rows.append((f'{nm} @ {nu:g} kHz', 'coupled Re', m['re'][j], float(c['cp_re_kHz'][s][0])))
            rows.append((f'{nm} @ {nu:g} kHz', 'coupled damping Hz', m['im'][j], -1e3 * float(c['cp_im_kHz'][s][0])))
            rows.append((f'{nm} @ {nu:g} kHz', 'background Re', gr_re_old * (1 + m['bg'][j] / 100), float(c['bg_re_kHz'][s][0])))
    return rows


def main():
    rows = compare()
    worst = {}
    print(f'{"mode":16s} {"quantity":20s} {"old":>14s} {"new":>14s} {"rel diff":>9s}')
    for nm, q, a, b in rows:
        rel = abs(b / a - 1)
        worst.setdefault(q, (0.0, ''))
        if rel > worst[q][0]:
            worst[q] = (rel, nm)
        flag = '   <-- above 1e-3' if rel >= TOL else ''
        print(f'{nm:16s} {q:20s} {a:14.6f} {b:14.6f} {rel:9.1e}{flag}')
    print('\nworst relative difference by quantity (and where):')
    for q, (v, nm) in worst.items():
        print(f'   {q:20s} {v:9.1e}  ({nm})')
    bad = [r for r in rows if abs(r[3] / r[2] - 1) >= TOL]
    print(f'\n{len(rows)} comparisons, {len(bad)} above {TOL:g}:', 'PASS' if not bad else 'FAIL')
    return not bad


if __name__ == '__main__':
    sys.exit(0 if main() else 1)
