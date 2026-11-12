"""Step 9: full cooling curves of a few representative cases, for the report figure (no new physics)."""
import json
from pathlib import Path

import numpy as np

from ns_cooling import sly4 as S
from ns_cooling import step3_integrals as I3
from ns_cooling import superfluid as SF
from ns_cooling.step7_sf_sweep import Star
from ns_cooling.structure import with_mass
from ns_cooling.tov import Table

RES = Path(__file__).parent / 'results'


def run():
    h8 = json.loads((RES / 'step8_inmedium.json').read_text())
    h7 = json.loads((RES / 'step7_sf_sweep.json').read_text())
    tab, fr, ncc = Table(*S.eos()), I3.sly_fractions(), S.composition()[0][0]
    t = np.logspace(1.5, 5.5, 160)
    stars = {}

    def star(M):
        if M not in stars:
            stars[M] = Star(with_mass(tab, M, hi=15.44, N=600), ncc, fr)
        return stars[M]
    cco = lambda pbf: next(s for s in h8['summary'] if s['object'].startswith('CXOU') and s['pbf'] == pbf)['best']
    cases = [('standard', 1.4, None, None, 'Y01', False), ('inmedium', 1.4, None, None, 'Y01', True)]
    for pbf in ('Y01', 'updated'):
        b = cco(pbf)
        cases.append((f'best_{pbf}', b['M'], b['p'], b['n'], pbf, True))
    curves = {}
    for key, M, gp, gn, pbf, im in cases:
        st = star(M)
        L = st.luminosity(t, *st.tables(SF.GAPS_P.get(gp), SF.GAPS_NT.get(gn), pbf, inmedium=im))
        curves[key] = {'M': M, 'p': gp, 'n': gn, 'pbf': pbf, 'inmedium': im, 'L': L.tolist()}
    # consistency with the sweeps at the source ages
    ia = {a: i for i, a in enumerate(h7['ages'])}
    a = 5000
    ref = {'standard': h7['base']['1.4'][ia[a]], 'inmedium': h8['base']['1.4'][ia[a]]}
    for pbf in ('Y01', 'updated'):
        b = cco(pbf)
        ref[f'best_{pbf}'] = next(r['L'][ia[a]] for r in h8['rows'] if (r['M'], r['p'], r['n'], r['pbf']) == (b['M'], b['p'], b['n'], pbf))
    dev = {k: float(np.exp(np.interp(np.log(a), np.log(t), np.log(curves[k]['L']))) / ref[k] - 1) for k in curves}
    checks = [{'name': 'Curves agree with the sweep values at 5000 yr (within 1 %)', 'ok': bool(max(map(abs, dev.values())) <= 0.01),
               'note': f"max {max(map(abs, dev.values())):.2%}"}]
    (RES / 'step9_figure_curves.json').write_text(json.dumps({'t_yr': t.tolist(), 'curves': curves, 'checks': checks}))
    for c in checks:
        print('PASS' if c['ok'] else 'FAIL', c['name'], f"({c['note']})")


if __name__ == '__main__':
    run()
