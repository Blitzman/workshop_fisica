"""N4 entry point:  python -m N4_modes_radial.run [--part fast|all|22,26,15,19,20,16] [--procs N] [--n0 48] [--n-max 384] [--out DIR] ...

Radial (l = 0) modes of the star with the axion condensate in full GR, table-based (old steps 15, 16, 19, 20, 22 and 26; step 17, the Font et al. comparison, is not part of this module).

Reads : D0 outputs/background.eq.json and D1 outputs/radial_linear.eq.json    the equations (the radial system is derived from them, as the legacy radial.py did from pipeline)
        S1 outputs/eos_bsk26.csv, star_noaxion.csv, star_nua<nu>.csv, manifest.json            the 12 stars of Paper II (steps 20, 22, 26)
        S1 outputs_extra/*.csv (python -m S1_star_background.run --extra, about 2 minutes)      nu_a = 3.4, 3.5, the step 19 values, the p_c = 1000 star with condensate and the
                                                                                               no-axion stars for p_c = 700..1100 (steps 15, 16, 19); NOT tracked by git
        R  outputs/paper1_tables.ref.json                                                       Paper I Tables 1-3 (step 26)
Writes: outputs/radial_contour.json (22), radial_paper1.json (26), radial_star.json (15), radial_crossing.json (19), radial_damped.json (20), radial_stability.json (16)
        (same keys as the legacy results/*.json) and outputs/manifest.json with the checks.
Cost  : symbolic compile about 20 s; then (one core, estimates from 0.11 s per det_entire evaluation): 22 and 26 minutes, 15 about 3 minutes, 19 about 3 minutes, 20 and 16
        tens of minutes each (their contour sampling is thousands of det_entire evaluations: the legacy took 6,012 s and 9,561 s on 6 workers). --procs divides it (fork pool).
--part fast = 22, 26, 15; --part analysis = 18, 21, 23, 24 (the old module A, extra parts of N4; they read the outputs of steps 15 and 20 from the output directory).  --n0/--n-max: sampling of the contours of steps 16 and 20 (smoke tests only: results with a smaller n0 are not comparable with the legacy).
"""
import argparse
import json
import os
import time
from pathlib import Path

from common import artifacts as A
from common.checks import Timer, report

from . import analysis as AN
from . import systems as SY
from . import validation as V

MODULE = 'N4_modes_radial'
HERE = Path(__file__).resolve().parent
ORDER = ['16', '22', '26', '15', '19', '20', '18', '21', '23', '24']
FILES = {'22': 'radial_contour.json', '26': 'radial_paper1.json', '15': 'radial_star.json', '19': 'radial_crossing.json', '20': 'radial_damped.json', '16': 'radial_stability.json',
         '18': 'radial_rydberg.json', '21': 'radial_deep.json', '23': 'radial_eos.json', '24': 'radial_scattering.json'}
TITLES = {'22': 'N4 / Paso 22: condición exterior frente al rayo complejo', '26': 'N4 / Paso 26: modos radiales frente a las Tablas 1-3 del Paper I',
          '15': 'N4 / Paso 15: modos radiales de la estrella con condensado (GR completa)', '19': 'N4 / Paso 19: cruce evitado F - átomo gravitatorio',
          '20': 'N4 / Paso 20: recuento en el semiplano amortiguado', '16': 'N4 / Paso 16: estabilidad radial por el principio del argumento',
          '18': 'N4 / Paso 18: ley de Rydberg del átomo gravitatorio', '21': 'N4 / Paso 21: ceros profundos frente al radio de empalme',
          '23': 'N4 / Paso 23: sensibilidad a la consistencia termodinámica de la EOS', '24': 'N4 / Paso 24: anchura de los modos desde la dispersión a frecuencia real'}


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument('--part', default='fast', help="'fast' (22, 26, 15), 'analysis' (18, 21, 23, 24: need the outputs of steps 15 and 20), 'all', or a comma list such as 22,26")
    ap.add_argument('--procs', type=int, default=min(6, os.cpu_count() or 1))
    ap.add_argument('--n0', type=int, default=48)
    ap.add_argument('--n-max', type=int, default=384)
    ap.add_argument('--d0-dir', default=str(HERE.parent / 'D0_derive_action' / 'outputs'))
    ap.add_argument('--d1-dir', default=str(HERE.parent / 'D1_derive_radial' / 'outputs'))
    ap.add_argument('--s1-dir', default=str(HERE.parent / 'S1_star_background' / 'outputs'))
    ap.add_argument('--extra-dir', default=str(HERE.parent / 'S1_star_background' / 'outputs_extra'))
    ap.add_argument('--r-dir', default=str(HERE.parent / 'R_paper_refs' / 'outputs'))
    ap.add_argument('--out', default=str(HERE / 'outputs'))
    args = ap.parse_args(argv)
    parts = ORDER if args.part == 'all' else (['22', '26', '15'] if args.part == 'fast' else (['18', '21', '23', '24'] if args.part == 'analysis' else [p for p in args.part.split(',') if p]))
    parts = [p for p in ORDER if p in parts]
    out, tm = Path(args.out), Timer()
    out.mkdir(parents=True, exist_ok=True)

    L, B, inputs = SY.load_inputs(args.d0_dir, args.d1_dir)
    tpath = Path(args.r_dir) / 'paper1_tables.ref.json'
    inputs.append(A.input_record(tpath))
    s1 = Path(args.s1_dir)
    inputs += [A.input_record(s1 / 'eos_bsk26.csv'), A.input_record(s1 / 'manifest.json')]
    inputs += [A.input_record(p) for p in sorted(s1.glob('star_*.csv'))]
    if (Path(args.extra_dir) / 'manifest.json').exists():
        inputs.append(A.input_record(Path(args.extra_dir) / 'manifest.json'))
    t0 = time.time()
    flat_i, flat_e, _, _ = SY.compile_flat(L, B)
    V.setup(flat_i, flat_e, str(s1 / 'eos_bsk26.csv'), s1, args.extra_dir, args.procs)
    print(f'  symbolic compile from the artifacts: {time.time() - t0:.1f} s', flush=True)

    res, arts, extra = [], [], {'parts': parts, 'procs': args.procs, 'n0': args.n0, 'n_max': args.n_max, 'seconds': {}}
    for p in parts:
        t0 = time.time()
        if p == '22':
            r_, o_ = V.step22()
        elif p == '26':
            r_, o_ = V.step26(V.load_paper1_tables(tpath))
        elif p == '15':
            r_, o_ = V.step15()
        elif p == '19':
            r_, o_ = V.step19()
        elif p == '20':
            r_, o_ = V.step20(args.n0, args.n_max)
        elif p == '16':
            r_, o_ = V.step16(args.n0, args.n_max)
        elif p == '18':
            r_, o_ = AN.step18(out)
        elif p == '21':
            r_, o_ = AN.step21(out)
        elif p == '23':
            r_, o_ = AN.step23(out)
        else:
            r_, o_ = AN.step24(out)
        extra['seconds'][p] = round(time.time() - t0, 1)
        print(f'  step {p}: {extra["seconds"][p]} s', flush=True)
        res += report(TITLES[p], r_)
        V.dump(out / FILES[p], o_)
        arts.append(out / FILES[p])
    # the manifest lists the checks of EVERY step whose JSON is in `out` (a partial run adds to what earlier runs wrote; each JSON stores its checks)
    res_all, arts_all = [], []
    for p in ORDER:
        f = out / FILES[p]
        if f.exists():
            res_all += [(c['name'], c['ok'], c['note']) for c in json.loads(f.read_text(encoding='utf-8'))['checks']]
            arts_all.append(f)
    extra['steps_in_manifest'] = [p for p in ORDER if (out / FILES[p]).exists()]
    A.write_manifest(out / 'manifest.json', MODULE, artifacts=arts_all, checks=res_all, seconds=tm.seconds(), inputs=inputs, extra=extra)
    print(f'N4 done in {tm.seconds()} s', flush=True)
    return res


if __name__ == '__main__':
    main()
