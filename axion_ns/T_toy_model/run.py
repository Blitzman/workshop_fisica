"""T entry point:  python -m T_toy_model.run [--part all|toy|branches|bridge] [--d1-dir DIR] [--out DIR]

Reads : D1 outputs/toy_model.eq.json   (only for the bridge check; toy and branches need no input)
Writes: outputs/toy_qnm_landscape.png, fig_toy.json       Fig. 3 of Paper I (old step 4) and its data
        outputs/radial_branches.json                      threshold, branches, avoided crossing (old step 14)
        outputs/manifest.json                             checks of old steps 4 and 14 + the bridge checks
Cost  : toy ~ minutes (mpmath landscape 600x250 twice), branches ~ 2 minutes, bridge seconds.
The 3 known FAILs of old step 4 (m_a R, m_a^2 f_a^2, lambda^2 quoted in Paper I) are disagreements with the paper, not bugs.
"""
import argparse
from pathlib import Path

from common import artifacts as A
from common.checks import Timer, report

from . import branches, toy
from .bridge import bridge_checks

MODULE = 'T_toy_model'
HERE = Path(__file__).resolve().parent


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument('--part', choices=['all', 'toy', 'branches', 'bridge'], default='all')
    ap.add_argument('--d1-dir', default=str(HERE.parent / 'D1_derive_radial' / 'outputs'))
    ap.add_argument('--out', default=str(HERE / 'outputs'))
    args = ap.parse_args(argv)
    out, tm = Path(args.out), Timer()
    out.mkdir(parents=True, exist_ok=True)
    toy.OUT = branches.OUT = out
    res, arts, inputs = [], [], []
    if args.part in ('all', 'bridge'):
        p = Path(args.d1_dir) / 'toy_model.eq.json'
        inputs.append(A.input_record(p))
        res += report('T / Puente ecuación -> numérico (D1 -> T)', bridge_checks(p))
    if args.part in ('all', 'toy'):
        res += report('T / Paso 4: números y figura', toy.run())
        arts += [out / 'toy_qnm_landscape.png', out / 'fig_toy.json']
    if args.part in ('all', 'branches'):
        res += report('T / Paso 14: modos radiales, umbral cuasi-normal/evanescente y ramas (modelo simplificado)',
                      branches.run())
        arts.append(out / 'radial_branches.json')
    A.write_manifest(out / 'manifest.json', MODULE, artifacts=arts, checks=res, seconds=tm.seconds(), inputs=inputs,
                     extra={'part': args.part})
    print(f'T done in {tm.seconds()} s', flush=True)
    return res


if __name__ == '__main__':
    main()
