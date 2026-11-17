"""Retypeset equation artifacts to PDF from the JSON alone (no derivation is rerun).

    python -m common.make_pdf D2_derive_nonradial/outputs/fullgr_extended.eq.json [more.eq.json ...] [--max-terms N]
    python -m common.make_pdf --all            # every *.eq.json under the module outputs folders

The JSON is the single source: this is also a check that the PDF can always be rebuilt from it.
"""
import argparse
from pathlib import Path

from .artifacts import ROOT, load_equations
from .pdf import pdf_from_artifact


def retypeset(path, max_terms=60, per_line=2):
    path = Path(path)
    art = load_equations(path)
    pdf = path.with_name(path.name.replace('.eq.json', '.pdf'))
    r = pdf_from_artifact(pdf, art.meta, dict(art.eq), per_line=per_line, max_terms=max_terms)
    print(f'{path.relative_to(ROOT) if path.is_absolute() and ROOT in path.parents else path}: pdf ok={r["ok"]} '
          f'max_terms={r["max_terms_used"]}')
    return r


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('files', nargs='*')
    ap.add_argument('--all', action='store_true')
    ap.add_argument('--max-terms', type=int, default=60)
    a = ap.parse_args()
    files = sorted(ROOT.glob('*/outputs/*.eq.json')) if a.all else [Path(f) for f in a.files]
    for f in files:
        retypeset(f, a.max_terms)
