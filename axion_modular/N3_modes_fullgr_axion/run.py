"""N3 entry point:  python -m N3_modes_fullgr_axion.run [--part all|10|25] [--quick] [--procs N] [--d2-dir D] [--s1-dir D] [--r-dir D] [--n1-dir D] [--out D]

Modes with the axion and the dynamic metric, beyond Cowling (old steps 10 and 25), table-based.

Reads : D2 outputs/fullgr_extended.eq.json                the extended system (verified against the action)
        S1 outputs/eos_bsk26.csv, star_noaxion.csv, star_nua<nu>.csv, manifest.json
        R  outputs/paper2_tables.ref.json                 Paper II Tables I and II (starting points)
        N1 outputs/modes_cowling.csv                      (optional) no-axion Cowling modes and Cowling shifts, for two of the checks
Writes: outputs/modes_axion.csv(+meta)        f, p1, p2 for every nu_a: background only (delta theta off) and coupled, Re and Im in kHz, plus the no-axion references
        outputs/exterior_convergence.json     step 25 details (part 25)
        outputs/manifest.json                 the checks of steps 10 and 25 (same texts where they exist)
Cost  : symbolic compile about 25 s; about 4 s per mode (about 70 modes in step 10, 26 in step 25); --procs divides it (fork pool).
--quick: only the 6 kHz star (6 modes), no convergence study; about 1 minute.
Not ported: the LMI polytrope end-to-end check (needs a polytrope EOS table) and the comparison with N2's matrix (tests/compare_n3_n2_matrix.py).
"""
import argparse
import json
import os
from pathlib import Path

import numpy as np

from common import artifacts as A
from common.checks import Timer, report

from . import systems as SY
from . import validation as V

MODULE = 'N3_modes_fullgr_axion'
HERE = Path(__file__).resolve().parent


def _n1_info(n1_dir):
    p = Path(n1_dir) / 'modes_cowling.csv'
    if not p.exists():
        return None, None
    cols, meta = A.load_table(p)
    nu0 = cols['nu_a_kHz'][0]
    s0 = cols['nu_a_kHz'] == nu0
    cow0 = [float(cols['nu_noaxion_kHz'][s0][k]) for k in range(3)]
    f = cols['mode'] == 'f'
    sh = 100 * (cols['nu_cp_re_kHz'][f] / cols['nu_noaxion_kHz'][f] - 1)
    p2 = cols['mode'] == 'p2'
    p2im = {float(n): abs(float(z)) * 1e3 for n, z in zip(cols['nu_a_kHz'][p2], cols['nu_cp_im_kHz'][p2])}
    return {'cow0': cow0, 'f_cow_max': float(np.nanmax(sh)), 'p2_cow_im': p2im}, A.input_record(p)


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument('--part', choices=['all', '10', '25'], default='all')
    ap.add_argument('--quick', action='store_true')
    ap.add_argument('--procs', type=int, default=min(6, os.cpu_count() or 1))
    ap.add_argument('--d2-dir', default=str(HERE.parent / 'D2_derive_nonradial' / 'outputs'))
    ap.add_argument('--s1-dir', default=str(HERE.parent / 'S1_star_background' / 'outputs'))
    ap.add_argument('--r-dir', default=str(HERE.parent / 'R_paper_refs' / 'outputs'))
    ap.add_argument('--n1-dir', default=str(HERE.parent / 'N1_modes_cowling' / 'outputs'))
    ap.add_argument('--out', default=str(HERE / 'outputs'))
    args = ap.parse_args(argv)
    out, tm = Path(args.out), Timer()
    out.mkdir(parents=True, exist_ok=True)

    ext, rec = SY.load_ext(args.d2_dir)
    inputs = [rec]
    tpath = Path(args.r_dir) / 'paper2_tables.ref.json'
    inputs.append(A.input_record(tpath))
    T1, T2 = _load_tables(tpath)
    man = json.loads((Path(args.s1_dir) / 'manifest.json').read_text(encoding='utf-8'))
    nus = [6.0] if args.quick else [float(x) for x in man['extra']['nu_a_list']]
    inputs += [A.input_record(Path(args.s1_dir) / 'eos_bsk26.csv'), A.input_record(Path(args.s1_dir) / 'star_noaxion.csv'),
               A.input_record(Path(args.s1_dir) / 'manifest.json')] + [A.input_record(Path(args.s1_dir) / f'star_nua{nu:g}.csv') for nu in nus]
    cow_n1, rec1 = _n1_info(args.n1_dir)
    if rec1:
        inputs.append(rec1)

    import time
    t0 = time.time()
    flats = SY.compile_flat(ext)
    subs = SY.substitution_checks()
    V.setup(flats, T1, T2, str(Path(args.s1_dir) / 'eos_bsk26.csv'), args.s1_dir)
    print(f'  symbolic compile + substitution checks: {time.time() - t0:.1f} s', flush=True)

    res, arts, extra = [], [], {'quick': args.quick, 'procs': args.procs, 'matching_radius_rule': 'R + 0.7 (r_end - R)'}
    rows = None
    if args.part in ('all', '10'):
        t0 = time.time()
        ref_gr, rows, conv = V.compute(nus, args.procs, args.quick)
        print(f'  step 10: {time.time() - t0:.1f} s (procs = {args.procs})', flush=True)
        res += report('N3 / Paso 10: con axión y métrica dinámica (etapa B)', V.checks(ref_gr, rows, conv, subs, cow_n1))
        cols = {k: [] for k in ('nu_a_kHz', 'mode', 'bg_re_kHz', 'bg_im_kHz', 'cp_re_kHz', 'cp_im_kHz', 'gr0_re_kHz', 'gr0_im_kHz')}
        for nu in sorted(rows):
            for i, nm in enumerate(V.MODES):
                for k, v in zip(cols, (nu, nm, rows[nu]['bg'][i].real, rows[nu]['bg'][i].imag, rows[nu]['cp'][i].real, rows[nu]['cp'][i].imag,
                                       ref_gr[i].real, ref_gr[i].imag)):
                    cols[k].append(v)
        cols = {k: (np.array(v, dtype=str) if k == 'mode' else np.array(v, dtype=float)) for k, v in cols.items()}
        A.save_table(out / 'modes_axion.csv', cols, name='modes_axion',
                     units={'nu_a_kHz': 'kHz', **{k: 'kHz' for k in cols if k.endswith('kHz') and k != 'nu_a_kHz'}},
                     descriptions={'bg_re_kHz': 'background only (delta theta off), Re nu', 'cp_re_kHz': 'coupled mode, Re nu',
                                   'cp_im_kHz': 'Im nu of the coupled mode (negative = damped); damping in Hz = -1e3 Im',
                                   'gr0_re_kHz': 'no-axion mode with the dynamic metric (reference)'},
                     parameters={'l': V.L, 'rtol': 1e-10, 'matching_radius': 'R + 0.7 (r_end - R)', 'quick': args.quick},
                     module=MODULE, script='run.py', inputs=inputs, title='Modes with the axion and the dynamic metric')
        arts.append(out / 'modes_axion.csv')
        extra['step10_seconds'] = {f'{nu:g}': round(rows[nu]['secs'], 1) for nu in rows}
    if args.part in ('all', '25') and not args.quick:
        if rows is None:
            c, _ = A.load_table(out / 'modes_axion.csv')
            rows = {float(nu): {'cp': [complex(float(c['cp_re_kHz'][(c['nu_a_kHz'] == nu) & (c['mode'] == m)][0]), float(c['cp_im_kHz'][(c['nu_a_kHz'] == nu) & (c['mode'] == m)][0]))
                                       for m in V.MODES]} for nu in sorted(set(c['nu_a_kHz']))}
        t0 = time.time()
        r25, o25 = V.step25(rows, args.procs)
        print(f'  step 25: {time.time() - t0:.1f} s', flush=True)
        res += report('N3 / Paso 25: rama del exterior y convergencia de los modos delicados', r25)
        (out / 'exterior_convergence.json').write_text(json.dumps(o25, indent=1, default=str))
        arts.append(out / 'exterior_convergence.json')
    A.write_manifest(out / 'manifest.json', MODULE, artifacts=arts, checks=res, seconds=tm.seconds(), inputs=inputs, extra=extra)
    print(f'N3 done in {tm.seconds()} s', flush=True)
    return res


def _load_tables(path):
    d = json.loads(Path(path).read_text(encoding='utf-8'))
    key = lambda k: None if k == 'none' else float(k)
    t1 = {key(k): tuple(v) for k, v in d['tables']['T1']['nu_a_kHz_rows'].items() if k != 'none'}
    t2 = {}
    for k, v in d['tables']['T2']['nu_a_kHz_rows'].items():
        t2[key(k)] = ((v[0], v[1]), (v[2], v[3]), (v[4], v[5]))
    return t1, t2


if __name__ == '__main__':
    main()
