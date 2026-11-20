"""N1 entry point:  python -m N1_modes_cowling.run [--part all|tables|atom] [--nu-a 6,8,...] [--s1-dir DIR] [--r-dir DIR] [--out DIR]

Reads : S1 outputs/eos_bsk26.csv, star_noaxion.csv, star_nua<nu>.csv, manifest.json   (part `tables`)
        R  outputs/paper2_tables.ref.json                                             (both parts)
Writes: outputs/modes_cowling.csv(+meta)       full-precision f, p1, p2: no axion, background only, coupled (Re, Im)
        outputs/atom_levels.csv(+meta), paper2_atom_levels.png, levels.json           (part `atom`)
        outputs/manifest.json                  checks of the old steps 9 and 6
Default --nu-a: every star listed in S1's manifest. With the S1 witness only nu_a = 6 kHz is available.
The per-star cost is dominated by the coupled-mode searches (minutes per star on one core).
"""
import argparse
import json
from pathlib import Path

import numpy as np

from common import artifacts as A
from common.checks import Timer, report

from . import atom_levels as AL
from . import run_tables as RT
from .star_table import StarTable
from .tables import load_paper2_tables, step9_shapes

MODULE = 'N1_modes_cowling'
HERE = Path(__file__).resolve().parent


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument('--part', choices=['all', 'tables', 'atom'], default='all')
    ap.add_argument('--nu-a', default=None, help='comma list of nu_a (kHz); default: all stars in the S1 manifest')
    ap.add_argument('--s1-dir', default=str(HERE.parent / 'S1_star_background' / 'outputs'))
    ap.add_argument('--r-dir', default=str(HERE.parent / 'R_paper_refs' / 'outputs'))
    ap.add_argument('--out', default=str(HERE / 'outputs'))
    args = ap.parse_args(argv)
    s1, out, tm = Path(args.s1_dir), Path(args.out), Timer()
    out.mkdir(parents=True, exist_ok=True)
    tpath = Path(args.r_dir) / 'paper2_tables.ref.json'
    T1, T2, T3, T4 = load_paper2_tables(tpath)
    inputs = [A.input_record(tpath)]
    res, arts, extra = [], [], {}

    if args.part in ('all', 'tables'):
        man = json.loads((s1 / 'manifest.json').read_text(encoding='utf-8'))
        nus = [float(x) for x in args.nu_a.split(',')] if args.nu_a else list(man['extra']['nu_a_list'])
        eos = s1 / 'eos_bsk26.csv'
        inputs += [A.input_record(eos), A.input_record(s1 / 'star_noaxion.csv'), A.input_record(s1 / 'manifest.json')]
        for nu in nus:
            inputs.append(A.input_record(s1 / f'star_nua{nu:g}.csv'))
        ref_star = StarTable(eos, s1 / 'star_noaxion.csv')
        star_for = lambda nu: StarTable(eos, s1 / f'star_nua{nu:g}.csv')
        T1p, T2p, T4M = step9_shapes(T1, T2, T4)
        ref, rows = RT.compute(star_for, nus, set(T1p), set(T2p), ref_star)
        cond = man['extra'].get('condensate_exists', {})
        r9 = RT.checks(ref, rows, cond.get('2.5') is False if '2.5' in cond else None, T1p, T2p, T4M)
        res += report('N1 / Paso 9: notas con axión en Cowling (Tablas I, II y IV del Paper II)', r9)
        cols = RT.table_columns(ref, rows)
        A.save_table(out / 'modes_cowling.csv', cols, name='modes_cowling',
                     units={'nu_a_kHz': 'kHz', 'nu_noaxion_kHz': 'kHz', 'nu_bg_kHz': 'kHz', 'nu_cp_re_kHz': 'kHz',
                            'nu_cp_im_kHz': 'kHz', 'M_msun': 'Msun', 'R_km': 'km'},
                     descriptions={'nu_bg_kHz': 'mode with the axion background only (delta theta off)',
                                   'nu_cp_re_kHz': 'Re nu of the coupled mode (nan if not computed: Table II stars only)',
                                   'nu_cp_im_kHz': 'Im nu of the coupled mode (negative = damped)'},
                     parameters={'l': RT.L, 'reference_noaxion': ref, 'stars': {f'{k:g}': {'M': v['M'], 'R': v['R'],
                                 'thc': v['thc'], 'seconds': v['secs']} for k, v in rows.items()}},
                     module=MODULE, script='run.py', inputs=inputs, title='Cowling modes of the stars with condensate')
        arts.append(out / 'modes_cowling.csv')
        extra['nu_a_computed'] = nus

    if args.part in ('all', 'atom'):
        AL.init((T1, T2, T3, T4), out)
        r6 = AL.numbers_checks()
        lv, rows6 = AL.levels_checks()
        r6 += lv
        path = AL.figure(rows6)
        r6.append((f'Figura guardada en results/{path.name}', path.exists(), ''))
        res += report('N1 / Paso 6: Paper II, números y niveles ligados', r6)
        A.save_table(out / 'atom_levels.csv', AL.levels_columns(rows6), name='atom_levels',
                     units={'nu_a_kHz': 'kHz', 'R_km': 'km', 'nu_n_kHz': 'kHz'},
                     parameters={'variants': AL.VARIANTS, 'l': AL.L_}, module=MODULE, script='run.py',
                     inputs=inputs[:1], title='Gravitational-atom levels from the Schwarzschild exterior (Table IV)')
        arts += [out / 'atom_levels.csv', path, out / 'levels.json']

    A.write_manifest(out / 'manifest.json', MODULE, artifacts=arts, checks=res, seconds=tm.seconds(), inputs=inputs,
                     extra=extra)
    print(f'N1 done in {tm.seconds()} s', flush=True)
    return res


if __name__ == '__main__':
    main()
