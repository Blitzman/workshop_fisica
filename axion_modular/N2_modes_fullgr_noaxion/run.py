"""N2 entry point:  python -m N2_modes_fullgr_noaxion.run [--quick] [--procs N] [--cowling-rtol 1e-11] [--d2-dir DIR] [--r-dir DIR] [--out DIR]

Validation of the extended Lindblom-Detweiler system with the axion switched off against the literature (old step 8, without the axion
background and the figures).

Reads : D2 outputs/fullgr_extended.eq.json      the extended system, verified against the action
        R  outputs/paper2_equations.eq.json     (8a)-(8c) of Paper II (Cowling)
        R  outputs/literature_modes.ref.json    Lindblom-Mendell-Ipser 1997 and Yoshida-Kojima 1997 mode frequencies
Writes: outputs/modes_noaxion.csv(+meta)        every mode computed, with the literature value
        outputs/scan_amplitude.csv(+meta)       |A_out|^2 against omega (l = 2, LMI star), extended LD and ZL as printed (full run only)
        outputs/manifest.json                   the 7 checks of the old step 8 (full run) or 4 smoke checks (--quick)

Cost  : the symbolic solve (about 40 s) plus about 5,800 ODE solves; with CSE about 1,000 s of CPU, divided by --procs (default
        min(6, cores); the pool needs the 'fork' start method, otherwise it runs serially). --quick: LMI star only, about 2 minutes.
Tolerances: rtol 1e-11 and atol 1e-14 everywhere (the amplitude must not be loosened, see validation.py); --cowling-rtol relaxes
        only the Cowling root search (try 1e-8: 2.2x faster per evaluation, and compare the 24 Yoshida-Kojima roots).
"""
import argparse
import os
import time
from pathlib import Path

import numpy as np

from common import artifacts as A
from common.checks import Timer, report

from . import systems as S2
from . import validation as V
from .refs import load_literature

MODULE = 'N2_modes_fullgr_noaxion'
HERE = Path(__file__).resolve().parent


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument('--quick', action='store_true', help='smoke run: LMI star only, l = 2')
    ap.add_argument('--procs', type=int, default=min(6, os.cpu_count() or 1))
    ap.add_argument('--cowling-rtol', type=float, default=1e-11, help='rtol of the Cowling root search only (default 1e-11)')
    ap.add_argument('--d2-dir', default=str(HERE.parent / 'D2_derive_nonradial' / 'outputs'))
    ap.add_argument('--r-dir', default=str(HERE.parent / 'R_paper_refs' / 'outputs'))
    ap.add_argument('--out', default=str(HERE / 'outputs'))
    args = ap.parse_args(argv)
    out, tm = Path(args.out), Timer()
    out.mkdir(parents=True, exist_ok=True)

    ext, p2, inputs = S2.load_inputs(args.d2_dir, args.r_dir)
    litpath = Path(args.r_dir) / 'literature_modes.ref.json'
    lit = load_literature(litpath)
    inputs.append(A.input_record(litpath))

    t0 = time.time()
    V.setup(ext, p2, lit, cow_rtol=args.cowling_rtol, quick=args.quick)
    t_setup = time.time() - t0
    print(f'  setup (symbolic solve + compile + stars): {t_setup:.1f} s', flush=True)
    t0 = time.time()
    cow, full, conv, scan, grid = V.compute(args.procs, args.quick)
    t_comp = time.time() - t0
    print(f'  computation: {t_comp:.1f} s (procs = {args.procs})', flush=True)

    res = V.quick_checks(cow, full, lit) if args.quick else V.checks(cow, full, conv, lit)
    report('N2: sistema extendido sin axión frente a la literatura (LMI 1997, Yoshida-Kojima 1997)', res)

    arts = []
    A.save_table(out / 'modes_noaxion.csv', V.mode_table(cow, full, lit), name='modes_noaxion',
                 units={'re_w': 'sqrt(M/R^3)', 'im_w': 'sqrt(M/R^3)', 'ref_re_w': 'sqrt(M/R^3)', 'ref_im_w': 'sqrt(M/R^3)'},
                 descriptions={'star': "'LMI' (1.4 Msun, kappa = 1.3346e5 cgs) or the compactness M/R",
                               'system': 'Cowling (8a-8c), extended LD (no axion) or ZL as printed (sign erratum)',
                               'ref_re_w': 'published value (nan if none; Cowling: Yoshida-Kojima Cowling column; extended: LMI or YK exact)'},
                 parameters={'quick': args.quick, 'rtol': 1e-11, 'atol': 1e-14, 'cowling_rtol': args.cowling_rtol},
                 module=MODULE, script='run.py', inputs=inputs, title='Mode frequencies of polytropes without axion')
    arts.append(out / 'modes_noaxion.csv')
    cols = {'system': np.array([w for w in scan for _ in grid]), 'omega_u': np.array([s for _ in scan for s in grid], dtype=float),
            'amp': np.concatenate([scan[w] for w in scan])}
    A.save_table(out / 'scan_amplitude.csv', cols, name='scan_amplitude', units={'omega_u': 'sqrt(M/R^3)'},
                 descriptions={'amp': '|A_out|^2 of the exterior Zerilli wave (l = 2, LMI star); noisy at the 1e-3..1e-2 level (see validation.py)'},
                 parameters={'quick': args.quick}, module=MODULE, script='run.py', inputs=inputs, title='Outgoing amplitude against frequency')
    arts.append(out / 'scan_amplitude.csv')
    A.write_manifest(out / 'manifest.json', MODULE, artifacts=arts, checks=res, seconds=tm.seconds(), inputs=inputs,
                     extra={'quick': args.quick, 'procs': args.procs, 'cowling_rtol': args.cowling_rtol, 'rtol': 1e-11, 'atol': 1e-14,
                            'setup_seconds': round(t_setup, 1), 'compute_seconds': round(t_comp, 1)})
    print(f'N2 done in {tm.seconds()} s', flush=True)
    return res


if __name__ == '__main__':
    main()
