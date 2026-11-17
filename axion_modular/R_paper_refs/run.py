"""R entry point:  python -m R_paper_refs.run [--out DIR] [--no-pdf]

Reads : nothing (the equations and tables are typed in paper1.py, paper2.py, ld_sources.py, tables.py).
Writes: outputs/paper1_equations.eq.json    Paper I (2510.16245): Eqs. (2.6a-d), (3.3a-e)      (+ .pdf)
        outputs/paper2_equations.eq.json    Paper II (2610.37308): Eqs. (8a-d)                  (+ .pdf)
        outputs/ld_sources.eq.json          Lindblom-Detweiler as typed in 4 sources + fixes    (+ .pdf)
        outputs/paper2_tables.ref.json      Paper II Tables I-IV as data
        outputs/manifest.json
These are REFERENCES (kind = 'transcribed'); they are compared against derived equations in D1/D2, never used to derive.
"""
import argparse
import json
from pathlib import Path

from common import artifacts as A
from common.checks import Timer, report
from common.pdf import pdf_from_artifact

from . import ld_sources as LD
from . import paper1, paper2, tables

MODULE = 'R_paper_refs'


def _ld_entries():
    src = {'ZL': ('Zhao & Lattimer, arXiv:2204.03037 Sec. II', LD.LD_ZL),
           'BU': ('Burgio et al., arXiv:1106.2736 appendix A', LD.LD_BU),
           'KU': ('Kunjipurayil et al., arXiv:2205.02081 Eqs. 9-14', LD.LD_KU),
           'ZH': ('Zhang et al., arXiv:2601.21911 appendix A', LD.LD_ZH),
           'METFIX': ('ZL with Q from the metric, sign of 2 w^2 W corrected (base of the extended system)', LD.LD_METFIX)}
    out = {}
    for tag, (ref, dct) in src.items():
        for k, e in dct.items():
            out[f'{tag}/{k}'] = {'expr': e, 'label': f'{tag}: equation for {k}', 'ref': ref,
                                 'note': 'form: typed equation lhs - rhs = 0'}
    out["ZL_FIX/X'"] = {'expr': LD.LD_ZL_FIX["X' con -2w^2"], 'label': "ZL X' with -2 w^2 in the W term",
                        'ref': 'Zhao & Lattimer with the sign erratum corrected', 'note': 'erratum fix'}
    for k, (msg, e) in LD.LD_BU_FIX.items():
        out[f'BU_FIX/{k}'] = {'expr': e, 'label': f'BU minimal correction for {k}', 'ref': 'Burgio et al. + fix',
                              'note': msg}
    return out


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument('--out', default=str(Path(__file__).parent / 'outputs'))
    ap.add_argument('--no-pdf', action='store_true')
    args = ap.parse_args(argv)
    out, tm = Path(args.out), Timer()

    sets = [
        ('paper1_equations', 'Paper I (arXiv:2510.16245): background (2.6) and radial perturbations (3.3)',
         {**{f'B{k}': {'expr': v, 'label': k, 'ref': 'Paper I'} for k, v in paper1.BACKGROUND.items()},
          **{f'P{k}': {'expr': v, 'label': k, 'ref': 'Paper I'} for k, v in paper1.PERTURBATIONS.items()}}),
        ('paper2_equations', 'Paper II (arXiv:2610.37308): nonradial perturbations in Cowling, Eqs. (8a-d)',
         {k: {'expr': v, 'label': k, 'ref': 'Paper II'} for k, v in paper2.PAPER.items()}),
        ('ld_sources', 'Lindblom-Detweiler equations as typed in four sources (with corrections)', _ld_entries()),
    ]
    res, arts = [], []
    expect = {'paper1_equations': 9, 'paper2_equations': 4, 'ld_sources': 5 * 6 + 1 + 5}
    for name, title, entries in sets:
        meta = A.save_equations(out / f'{name}.eq.json', name, entries, title=title, kind='transcribed',
                                description='Transcribed as typed in the source; not derived here.', module=MODULE,
                                script='run.py')
        back = A.load_equations(out / f'{name}.eq.json')
        ok = len(back.eq) == expect[name] and all(back.eq[k] == entries[k]['expr'] for k in entries)
        res.append((f'{name}: {len(entries)} equations written and read back identically', ok, ''))
        arts.append(out / f'{name}.eq.json')
        if not args.no_pdf:
            pr = pdf_from_artifact(out / f'{name}.pdf', meta, {k: v['expr'] for k, v in entries.items()},
                                   per_line=2, max_terms=80)
            print(f'  pdf {name}: ok={pr["ok"]}')
            if pr['ok']:
                arts.append(out / f'{name}.pdf')
    tj = out / 'paper2_tables.ref.json'
    tj.write_text(json.dumps(tables.as_json(), indent=1), encoding='utf-8')
    arts.append(tj)
    res.append(('paper2_tables.ref.json: Tables I-IV written (rows: %s)' % {k: len(v['nu_a_kHz_rows'])
                for k, v in tables.as_json()['tables'].items()}, True, ''))
    report('R: referencias del paper', res)
    A.write_manifest(out / 'manifest.json', MODULE, artifacts=arts, checks=res, seconds=tm.seconds())
    print(f'R done in {tm.seconds()} s')
    return res


if __name__ == '__main__':
    main()
