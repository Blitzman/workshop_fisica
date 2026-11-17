"""Modules must not import each other: only `common`, the standard library, third-party libraries and their own files."""
import ast
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
MODULES = ['R_paper_refs', 'D0_derive_action', 'D1_derive_radial', 'D2_derive_nonradial', 'S1_star_background',
           'N1_modes_cowling', 'T_toy_model', 'report']
THIRD = {'numpy', 'scipy', 'sympy', 'mpmath', 'matplotlib'}


def _imports(path):
    tree = ast.parse(path.read_text(encoding='utf-8'))
    for n in ast.walk(tree):
        if isinstance(n, ast.Import):
            for a in n.names:
                yield a.name.split('.')[0], 0
        elif isinstance(n, ast.ImportFrom):
            yield (n.module or '').split('.')[0], n.level


def test_no_cross_module_imports():
    stdlib = set(sys.stdlib_module_names)
    bad = []
    for m in MODULES:
        for f in (ROOT / m).glob('*.py'):
            for top, level in _imports(f):
                if level > 0 or top in stdlib | THIRD | {'common', ''}:
                    continue
                bad.append(f'{m}/{f.name} imports {top}')
    assert not bad, bad


def test_common_does_not_import_modules():
    bad = []
    for f in (ROOT / 'common').glob('*.py'):
        for top, level in _imports(f):
            if top in MODULES:
                bad.append(f'common/{f.name} imports {top}')
    assert not bad, bad


def test_derivation_files_do_not_import_references():
    """Rule of axion_ns/AGENTS.md section 1: derivation code never knows the equations typed from the papers."""
    for m, files in {'D0_derive_action': ['derive.py'], 'D1_derive_radial': ['derive.py'],
                     'D2_derive_nonradial': ['derive_cowling.py', 'derive_fullgr.py']}.items():
        for fn in files:
            src = (ROOT / m / fn).read_text(encoding='utf-8')
            assert 'R_paper_refs' not in src and 'paper1' not in src and 'paper2' not in src and 'ld_sources' not in src, (m, fn)
