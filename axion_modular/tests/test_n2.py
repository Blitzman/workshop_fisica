"""N2: fast regression of the Cowling path, literature data and the polytrope star (the 29-mode validation is `make N2`)."""
import numpy as np
import pytest

from common import artifacts as A
from tests.conftest import ROOT

P2 = ROOT / 'R_paper_refs' / 'outputs' / 'paper2_equations.eq.json'
LIT = ROOT / 'R_paper_refs' / 'outputs' / 'literature_modes.ref.json'
G, C, MSUN = 6.6743e-8, 2.99792458e10, 1.98841e33
MASS = 1.4 * G * MSUN / C**2 / np.sqrt(1.3346e5 / G)


def test_polytrope_star_by_compactness():
    from N2_modes_fullgr_noaxion.polytrope import Star
    st = Star('energy', 0.1)
    assert abs(st.M / st.R - 0.1) < 1e-10 and st.R > 0


@pytest.mark.skipif(not LIT.exists(), reason='R literature artifact not produced')
def test_literature_artifact_matches_the_typed_tables():
    from N2_modes_fullgr_noaxion.refs import load_literature
    from R_paper_refs import literature as L
    lmi, yk, rkm = load_literature(LIT)
    assert lmi == L.LMI97 and yk == L.YK97 and rkm == L.LMI_R_KM


@pytest.mark.skipif(not P2.exists(), reason='R paper2 artifact not produced')
def test_cowling_evaluation_regression():
    """X(R) at omega = 1.3 sqrt(M/R^3), l = 2, LMI star: value measured with the legacy code and with this module (agree to 3e-14)."""
    from N2_modes_fullgr_noaxion import modes as M2
    from N2_modes_fullgr_noaxion import systems as S2
    from N2_modes_fullgr_noaxion.polytrope import Star
    p2 = dict(A.load_equations(P2, expect_name='paper2_equations').eq)
    cow = S2.cowling_system(p2, 'energy')
    st = Star('energy', mass=MASS)
    assert abs(st.M - 0.146192668557) < 1e-11 and abs(st.R - 1.021808063106) < 1e-11
    u = np.sqrt(st.M / st.R**3)
    assert abs(M2.cowling_X(st, cow, 1.3 * u, 2) / 1.179858973801e-03 - 1) < 1e-9
    # the compiled A*y and the compiled matrix agree
    r0 = 0.5 * st.R
    b0 = st._start(st.nc, 1e-4 * st.R)
    b0[2] = st.Phic
    B = st.bg(r0, *b0)[0]
    y = np.array([1.0, 0.4])
    assert np.allclose(cow.ay(r0, 1.3 * u, 2, *B, *y), cow.mat(r0, 1.3 * u, 2, *B) @ y, rtol=1e-12, atol=0)
