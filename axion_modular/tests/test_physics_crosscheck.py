"""N1 re-implements the two short background functions (matter, rhs). They must agree with S1's Star (cross-module on purpose)."""
import numpy as np
import pytest

from tests.conftest import ROOT

S1 = ROOT / 'S1_star_background' / 'outputs'


def _pairs():
    out = []
    for name in ('star_noaxion', 'star_nua6'):
        if (S1 / f'{name}.csv').exists():
            out.append(name)
    return out


@pytest.mark.parametrize('name', _pairs() or [pytest.param('none', marks=pytest.mark.skip(reason='S1 not run'))])
def test_rhs_and_matter_agree_with_S1(name):
    from N1_modes_cowling.star_table import StarTable
    from S1_star_background import star as S
    st = StarTable(S1 / 'eos_bsk26.csv', S1 / f'{name}.csv')
    ref = object.__new__(S.Star)
    ref.nu_a, ref.p0, ref.E0, ref.nc, ref.ma, ref.fa2 = st.nu_a, st.p0, st.E0, st.nc, st.ma, st.fa2
    rng = np.random.default_rng(3)
    for r in rng.uniform(0.05, 0.98 * st.R, 60):
        y = st.y_at(r)
        a, b = np.array(st.rhs(r, y)), np.array(ref.rhs(r, list(y)))
        assert np.allclose(a, b, rtol=1e-12, atol=0), (r, a, b)
        assert np.allclose(st.matter(y[1], y[3]), ref.matter(y[1], y[3]), rtol=1e-10, atol=0)
