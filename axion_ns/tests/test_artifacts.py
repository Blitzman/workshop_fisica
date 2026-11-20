"""Round trip of equation and table artifacts, including awkward SymPy objects."""
import numpy as np
import sympy as sp

from common import artifacts as A
from common.notation import r, a0, n0, Phi0, eps0, V, l, fa, p, eps, dpdn, mstar


def test_expression_roundtrip(tmp_path):
    x = sp.Dummy('x')
    exprs = {
        'deriv': sp.Derivative(eps0(n0), (n0, 2)) * sp.diff(n0, r) + sp.Rational(3, 7) * sp.sqrt(2) * sp.pi,
        'p_and_dpdn': p + dpdn,
        'func': sp.exp(-2 * Phi0) * V(a0) / (l * fa**2) + sp.I * sp.sin(r) * sp.log(r) - sp.Float('0.25', 30),
        'subs': sp.Subs(sp.Derivative(eps0(x), x), x, n0) + mstar(a0),
        'zero': sp.Integer(0),
    }
    meta = A.save_equations(tmp_path / 't.eq.json', 't', exprs, module='test')
    back = A.load_equations(tmp_path / 't.eq.json')
    for k, v in exprs.items():
        assert sp.expand(back[k] - v) == 0, k
    assert back.sha256 == meta['content_sha256']


def test_tamper_detected(tmp_path):
    A.save_equations(tmp_path / 't.eq.json', 't', {'e': r + 1}, module='test')
    txt = (tmp_path / 't.eq.json').read_text()
    (tmp_path / 't.eq.json').write_text(txt.replace('"n",1', '"n",2', 1) if '"n",1' in txt else txt.replace('1', '2', 1))
    try:
        A.load_equations(tmp_path / 't.eq.json')
    except ValueError:
        return
    raise AssertionError('tampered artifact was accepted')


def test_table_roundtrip(tmp_path):
    cols = {'a': np.array([1.0, 2.5e-300, np.pi]), 'name': np.array(['f', 'p1', 'p2']), 'z': np.array([np.nan, 1.0, 2.0])}
    A.save_table(tmp_path / 't.csv', cols, name='t', module='test')
    back, meta = A.load_table(tmp_path / 't.csv')
    assert np.array_equal(back['a'], cols['a']) and list(back['name']) == ['f', 'p1', 'p2']
    assert np.isnan(back['z'][0]) and meta['n_rows'] == 3
