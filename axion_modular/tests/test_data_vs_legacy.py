"""Data files of T and N1 (atom levels) must equal the legacy committed ones."""
import json

import pytest

from tests.conftest import REPO, ROOT


def _diff(a, b, tol=1e-9, path=''):
    out = []
    if isinstance(a, dict) and isinstance(b, dict):
        for k in set(a) | set(b):
            out += [f'{path}/{k}: missing'] if (k not in a or k not in b) else _diff(a[k], b[k], tol, f'{path}/{k}')
    elif isinstance(a, list) and isinstance(b, list):
        out += [f'{path}: length'] if len(a) != len(b) else sum((_diff(x, y, tol, f'{path}[{i}]') for i, (x, y) in
                                                                  enumerate(zip(a, b))), [])
    elif isinstance(a, (int, float)) and isinstance(b, (int, float)) and not isinstance(a, bool):
        if abs(a - b) > tol * max(1, abs(a), abs(b)):
            out.append(f'{path}: {a} vs {b}')
    elif a != b:
        out.append(f'{path}: {a!r} vs {b!r}')
    return out


@pytest.mark.parametrize('new,old', [('T_toy_model/outputs/fig_toy.json', 'results/fig_toy.json'),
                                     ('T_toy_model/outputs/radial_branches.json', 'results/radial_branches.json'),
                                     ('N1_modes_cowling/outputs/levels.json', 'results/fig_levels.json')])
def test_equals_legacy(new, old):
    n, o = ROOT / new, REPO / old
    if not (n.exists() and o.exists()):
        pytest.skip('missing file')
    d = _diff(json.loads(n.read_text()), json.loads(o.read_text()))
    assert not d, d[:5]
