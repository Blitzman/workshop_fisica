"""D0 entry point:  python -m D0_derive_action.run [--out DIR] [--no-pdf]

Reads : nothing.
Writes: outputs/field_equations.eq.json (+ .pdf)   general (t, r) field equations from the action
        outputs/background.eq.json      (+ .pdf)   their static limit (TOV + Klein-Gordon form)
        outputs/manifest.json
"""
import argparse
from pathlib import Path

import sympy as sp

from common import artifacts as A
from common.checks import Timer, report
from common.notation import n0, a0, eps_of
from common.pdf import pdf_from_artifact

from .derive import background, derive_field_equations

MODULE = 'D0_derive_action'
LABELS = {'Gtt': 'G_tt - 8 pi G T_tt', 'Gtr': 'G_tr - 8 pi G T_tr', 'Grr': 'G_rr - 8 pi G T_rr',
          'KG': 'Klein-Gordon: box a + V\' + n m*\'', 'FLUID_r': 'radial momentum: div T(total) = 0, with KG eliminated',
          'BAR': 'baryon number conservation', 'ident_ax': 'identity: div T_ax^mu_r - (box a + V\') d_r a (must be 0)'}


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument('--out', default=str(Path(__file__).parent / 'outputs'))
    ap.add_argument('--no-pdf', action='store_true')
    args = ap.parse_args(argv)
    out, tm = Path(args.out), Timer()

    F = derive_field_equations()
    B = {k: background(v) for k, v in F.items()}

    # ---- checks that concern the derivation itself (old step 1, first and last checks)
    res = [('Identidad nabla_mu T^mu_r(axion) = (box a + V\') d_r a', F['ident_ax'] == 0, ''),
           ('Ec. G_tr del fondo estático se anula idénticamente', B['Gtr'] == 0, '')]
    nn = sp.Symbol('n')
    p_exact = (nn * sp.diff(eps_of(nn, a0), nn) - eps_of(nn, a0))
    res.append(('p = n d(eps)/dn - eps no depende de a (la Ec. (2.8) del paper es exacta)',
                sp.diff(p_exact, a0) == 0, 'el paper la presenta como aproximación (≈)'))
    report('D0: acción -> ecuaciones de campo y fondo', res)

    metas, arts = {}, []
    for name, E, title, desc in [
            ('field_equations', F, 'Field equations from the action (general fields of t, r)',
             'Raw (unsimplified) output of the derivation: Einstein, Klein-Gordon, fluid and baryon equations (each = 0) for the '
             'metric diag(e^{2Phi}, -e^{2lambda}, -r^2, -r^2 sin^2 theta) and fields Phi, lambda, a, n, xi of (t, r). '
             'Some terms are very long; the readable forms are the static limit (background.pdf) and the linearizations of D1/D2.'),
            ('background', B, 'Static background equations (TOV + Klein-Gordon form)',
             'Static limit of the field equations: Phi(t,r) -> Phi0(r), lambda -> lambda0, a -> a0, n -> n0, xi -> 0.')]:
        entries = {k: {'expr': v, 'label': LABELS[k]} for k, v in E.items()}
        meta = A.save_equations(out / f'{name}.eq.json', name, entries, title=title, kind='derived',
                                description=desc, module=MODULE, script='run.py',
                                checks=[{'name': n_, 'ok': ok} for n_, ok, _ in res])
        metas[name] = meta
        arts.append(out / f'{name}.eq.json')
        if not args.no_pdf:
            ex = {k: v for k, v in E.items()}
            pr = pdf_from_artifact(out / f'{name}.pdf', meta, ex, per_line=2, max_terms=60)
            print(f'  pdf {name}: ok={pr["ok"]} max_terms={pr["max_terms_used"]}')
            if pr['ok']:
                arts.append(out / f'{name}.pdf')
    A.write_manifest(out / 'manifest.json', MODULE, artifacts=arts, checks=res, seconds=tm.seconds())
    print(f'D0 done in {tm.seconds()} s')
    return res


if __name__ == '__main__':
    main()
