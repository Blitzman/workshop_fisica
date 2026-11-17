"""D2 entry point:  python -m D2_derive_nonradial.run [--d0-dir DIR] [--r-dir DIR] [--out DIR] [--no-pdf] [--stage S]

Reads : D0 outputs/background.eq.json                              (static background equations, from the action)
        R  outputs/paper2_equations.eq.json, ld_sources.eq.json    (references: Paper II (8a-d); LD as typed + fixes)
Writes: outputs/cowling_linear.eq.json    derived: Eqs. (8a)-(8d) linearized in Cowling + the axion energy flux T^r_t
        outputs/fullgr_linear.eq.json     derived: the 11 linearized equations with the perturbed metric (beyond Cowling)
        outputs/fullgr_extended.eq.json   verified: the extended system LD + axion (written only if its checks pass)
        outputs/manifest.json             with the checks of the old steps 5 and 7 (same names), each artifact has its PDF

Stages (default `all`): `cowling` (old step 5), `fullgr` (old step 7: derivation + checks, the slow one).
The derivation never reads R; R is read only by the checks and by the assembly of the candidate extended system.
"""
import argparse
from pathlib import Path
from types import SimpleNamespace

from common import artifacts as A
from common.checks import Timer, report
from common.modulo import Rules
from common.pdf import pdf_from_artifact

from . import checks_cowling as CC
from . import checks_fullgr as CF
from .derive_cowling import derive as derive_cowling_full, linearize as linearize_cowling
from .derive_fullgr import linear_system
from .ext_system import NEW, build_ext

MODULE = 'D2_derive_nonradial'
HERE = Path(__file__).resolve().parent


def _pdf(out, name, meta, exprs, enabled, arts, max_terms=40):
    if not enabled:
        return
    pr = pdf_from_artifact(out / f'{name}.pdf', meta, exprs, per_line=2, max_terms=max_terms)
    print(f'  pdf {name}: ok={pr["ok"]} max_terms={pr["max_terms_used"]}', flush=True)
    if pr['ok']:
        arts.append(out / f'{name}.pdf')


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument('--d0-dir', default=str(HERE.parent / 'D0_derive_action' / 'outputs'))
    ap.add_argument('--r-dir', default=str(HERE.parent / 'R_paper_refs' / 'outputs'))
    ap.add_argument('--out', default=str(HERE / 'outputs'))
    ap.add_argument('--no-pdf', action='store_true')
    ap.add_argument('--stage', choices=['all', 'cowling', 'fullgr'], default='all')
    args = ap.parse_args(argv)
    d0, rd, out, tm = Path(args.d0_dir), Path(args.r_dir), Path(args.out), Timer()
    pdf = not args.no_pdf

    b_path, p2_path, ld_path = d0 / 'background.eq.json', rd / 'paper2_equations.eq.json', rd / 'ld_sources.eq.json'
    B = A.load_equations(b_path, expect_name='background')
    P2 = A.load_equations(p2_path, expect_name='paper2_equations')
    LDa = A.load_equations(ld_path, expect_name='ld_sources')
    inputs = [A.input_record(b_path), A.input_record(p2_path), A.input_record(ld_path)]
    rules = Rules(B.eq)
    PAPER = dict(P2.eq)
    res, arts = [], []

    # ------------------------------------------------------------ stage cowling (old step 5)
    D = derive_cowling_full()
    L5 = {name: linearize_cowling(D[der], ang) for name, (der, ang, _) in CC.MATCH.items()}
    if args.stage in ('all', 'cowling'):
        res += report('D2 / Paso 5: Paper II, perturbaciones no radiales (simbólico)',
                      CC.equations_checks(L5, rules, PAPER) + CC.origin_checks(PAPER) + CC.appendix_checks()
                      + CC.energy_checks(D['Trt_ax']))
        entries = {k: {'expr': v, 'label': k, 'ref': 'derived from the action in Cowling'} for k, v in L5.items()}
        entries['Trt_ax'] = {'expr': D['Trt_ax'], 'label': 'T^r_t of the axion (energy flux), general fields of (t, r, theta)',
                             'ref': 'used for appendix B (B4)'}
        meta = A.save_equations(out / 'cowling_linear.eq.json', 'cowling_linear', entries,
                                title='Nonradial perturbations in Cowling, derived from the action', kind='derived',
                                description='Linearization (l generic, m = 0, factor (r/R)^l extracted, Legendre eliminated) '
                                            'of KG, radial and angular momentum and baryon conservation with a frozen metric.',
                                module=MODULE, script='run.py', inputs=inputs[:1],
                                checks=[{'name': n_, 'ok': ok} for n_, ok, _ in res])
        arts.append(out / 'cowling_linear.eq.json')
        _pdf(out, 'cowling_linear', meta, {k: v['expr'] for k, v in entries.items()}, pdf, arts)

    # ------------------------------------------------------------ stage fullgr (old step 7)
    if args.stage in ('all', 'fullgr'):
        print('  deriving the linear system with the perturbed metric (slow) ...', flush=True)
        L_raw = linear_system()
        ref = SimpleNamespace(
            LD_ZL={k.split('/', 1)[1]: v for k, v in LDa.eq.items() if k.startswith('ZL/')},
            LD_ZL_FIX={"X' con -2w^2": LDa.eq["ZL_FIX/X'"]},
            LD_BU={k.split('/', 1)[1]: v for k, v in LDa.eq.items() if k.startswith('BU/')},
            LD_BU_FIX={k.split('/', 1)[1]: (LDa.entries[k]['note'], v) for k, v in LDa.eq.items() if k.startswith('BU_FIX/')},
            LD_KU={k.split('/', 1)[1]: v for k, v in LDa.eq.items() if k.startswith('KU/')},
            LD_ZH={k.split('/', 1)[1]: v for k, v in LDa.eq.items() if k.startswith('ZH/')},
            NEW=NEW,
            EXT=build_ext({k.split('/', 1)[1]: v for k, v in LDa.eq.items() if k.startswith('METFIX/')},
                          PAPER['(8d) KG interior']))
        # keep the order of the original dictionaries (the checks compare lists of keys)
        order = ["H1'", "K'", "W'", "X'", 'H0', 'V']
        for nm in ('LD_ZL', 'LD_BU', 'LD_KU', 'LD_ZH'):
            setattr(ref, nm, {k: getattr(ref, nm)[k] for k in order})
        ref.LD_BU_FIX = {k: ref.LD_BU_FIX[k] for k in ["H1'", "K'", "X'", 'H0', 'V']}
        CF.init(L_raw, rules, B.eq, L5, ref)
        r7 = CF.run_all(report)
        res += r7
        entries = {k: {'expr': v, 'label': f'{k} (linearized, perturbed metric)',
                       'ref': 'derived from the action; angular factor: ' + ('dY' if k in ('Etth', 'Erth', 'FLUID_th') else
                                                                             'traceless' if k == 'Etl' else 'Y')}
                   for k, v in L_raw.items()}
        meta = A.save_equations(out / 'fullgr_linear.eq.json', 'fullgr_linear', entries,
                                title='Nonradial perturbations beyond Cowling, derived from the action', kind='derived',
                                description='Linearized Einstein (tt, tr, rr, t-theta, r-theta, angular trace and traceless '
                                            'parts), KG, fluid and baryon equations with the Regge-Wheeler perturbed metric '
                                            '(H0, H1, H2, K) in the Lindblom-Detweiler convention.',
                                module=MODULE, script='run.py', inputs=inputs[:1],
                                checks=[{'name': n_, 'ok': ok} for n_, ok, _ in r7])
        arts.append(out / 'fullgr_linear.eq.json')
        _pdf(out, 'fullgr_linear', meta, {k: v['expr'] for k, v in entries.items()}, pdf, arts, max_terms=40)

        ext_ok = all(ok for n_, ok, _ in r7 if n_.startswith('Sistema extendido'))
        if ext_ok:
            ee = {k: {'expr': v, 'label': f'extended system, equation for {k}',
                      'ref': 'LD (metric Q, corrected sign) + new axion terms; verified against the derived system',
                      'note': 'form: LD_equation - NEW = 0'} for k, v in ref.EXT.items()}
            meta = A.save_equations(out / 'fullgr_extended.eq.json', 'fullgr_extended', ee,
                                    title='Extended system: Lindblom-Detweiler + axion (verified)', kind='candidate',
                                    description='Compact form used by the integrators of the beyond-Cowling modes. Every '
                                                'equation was checked to follow from the action (checks of this module).',
                                    module=MODULE, script='run.py', inputs=inputs,
                                    checks=[{'name': n_, 'ok': ok} for n_, ok, _ in r7])
            arts.append(out / 'fullgr_extended.eq.json')
            _pdf(out, 'fullgr_extended', meta, {k: v['expr'] for k, v in ee.items()}, pdf, arts, max_terms=40)
        else:
            print('  !! extended system checks failed: fullgr_extended.eq.json NOT written', flush=True)

    A.write_manifest(out / 'manifest.json', MODULE, artifacts=arts, checks=res, seconds=tm.seconds(), inputs=inputs,
                     extra={'stage': args.stage})
    print(f'D2 done in {tm.seconds()} s', flush=True)
    return res


if __name__ == '__main__':
    main()
