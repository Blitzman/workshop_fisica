"""D1 entry point:  python -m D1_derive_radial.run [--d0-dir DIR] [--r-dir DIR] [--out DIR] [--no-pdf]

Reads : D0 outputs/field_equations.eq.json and background.eq.json   (derived from the action)
        R  outputs/paper1_equations.eq.json                           (Paper I equations, as typed)
Writes: outputs/radial_linear.eq.json (+ .pdf)   radial linearization of the field equations (Gtt, Gtr, Grr, KG, FLUID_r, BAR)
        outputs/toy_model.eq.json     (+ .pdf)   simplified model of Sec. IV derived from (3.3): (4.2), (4.5), (4.12), (4.13)
        outputs/manifest.json                    with the checks of the old steps 1, 2 and 3 (same names)
"""
import argparse
from pathlib import Path

from common import artifacts as A
from common.checks import Timer, report
from common.modulo import Rules
from common.pdf import pdf_from_artifact

from .checks_radial import background_checks, perturbation_checks
from .derive import derive_linear
from .toy import run_toy

MODULE = 'D1_derive_radial'
HERE = Path(__file__).resolve().parent


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument('--d0-dir', default=str(HERE.parent / 'D0_derive_action' / 'outputs'))
    ap.add_argument('--r-dir', default=str(HERE.parent / 'R_paper_refs' / 'outputs'))
    ap.add_argument('--out', default=str(HERE / 'outputs'))
    ap.add_argument('--no-pdf', action='store_true')
    args = ap.parse_args(argv)
    d0, rd, out, tm = Path(args.d0_dir), Path(args.r_dir), Path(args.out), Timer()

    f_path, b_path, p_path = d0 / 'field_equations.eq.json', d0 / 'background.eq.json', rd / 'paper1_equations.eq.json'
    F = A.load_equations(f_path, expect_name='field_equations')
    B = A.load_equations(b_path, expect_name='background')
    P1 = A.load_equations(p_path, expect_name='paper1_equations')
    inputs = [A.input_record(f_path), A.input_record(b_path), A.input_record(p_path)]
    PB = {k[1:]: v for k, v in P1.eq.items() if k.startswith('B')}
    PP = {k[1:]: v for k, v in P1.eq.items() if k.startswith('P')}

    # ---- derive: radial linearization of the field equations from the action
    L = derive_linear(F.eq)
    rules = Rules(B.eq)

    # ---- checks of old step 1 (background) and step 2 (perturbations)
    res = []
    res += report('D1 / Paso 1: fondo (TOV + KG) frente a la Ec. (2.6)', background_checks(B.eq, rules, PB))
    res += report('D1 / Paso 2: perturbaciones radiales', perturbation_checks(L, rules, PP))
    # ---- toy model (old step 3)
    toy, rt = run_toy(PP)
    res += report('D1 / Paso 3: modelo simplificado (Sec. IV)', rt)

    arts = []
    for name, title, kind, desc, eqs, labels in [
        ('radial_linear', 'Radial linearization of the field equations', 'derived',
         'Linear part in (deltaPhi, deltalambda, deltaa, deltan, xi1) with e^{-i omega t} factored out. '
         'Each equation is the linearization of the corresponding D0 field equation.',
         L, {'Gtt': 'tt', 'Gtr': 'tr', 'Grr': 'rr', 'KG': 'Klein-Gordon', 'FLUID_r': 'Euler (radial momentum)',
             'BAR': 'continuity (baryon number)'}),
        ('toy_model', 'Simplified model of Sec. IV (Paper I), derived from (3.3)', 'derived',
         'Planar, gradient-free limit of (3.3) with the approximations (4.1); dimensionless form, matching '
         'determinant f(omega) of Eq. (4.13) in tan form (symbols k_p, k_m, k_a, k_out, x_0). '
         'The roots k_pm^2 / Delta are verified to solve det(M - k^2) = 0.',
         toy, {k: k for k in toy})]:
        entries = {k: {'expr': v, 'label': labels.get(k, k)} for k, v in eqs.items()}
        meta = A.save_equations(out / f'{name}.eq.json', name, entries, title=title, kind=kind, description=desc,
                                module=MODULE, script='run.py', inputs=inputs,
                                checks=[{'name': n_, 'ok': ok} for n_, ok, _ in res])
        arts.append(out / f'{name}.eq.json')
        if not args.no_pdf:
            pr = pdf_from_artifact(out / f'{name}.pdf', meta, dict(eqs), per_line=2, max_terms=60)
            print(f'  pdf {name}: ok={pr["ok"]} max_terms={pr["max_terms_used"]}')
            if pr['ok']:
                arts.append(out / f'{name}.pdf')
    A.write_manifest(out / 'manifest.json', MODULE, artifacts=arts, checks=res, seconds=tm.seconds(), inputs=inputs)
    print(f'D1 done in {tm.seconds()} s')
    return res


if __name__ == '__main__':
    main()
