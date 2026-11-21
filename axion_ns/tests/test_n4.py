"""N4: fast tests (no symbolic compile): the Paper I tables artifact, the legacy winding function, the exterior wave on an analytic exterior matrix, and the coefficient tables with a
synthetic 1/r matrix."""
import numpy as np
import pytest

from tests.conftest import ROOT

R_TABLES = ROOT / 'R_paper_refs' / 'outputs' / 'paper1_tables.ref.json'


@pytest.mark.skipif(not R_TABLES.exists(), reason='R paper1 artifact not produced')
def test_paper1_tables_artifact_matches_the_typed_tables():
    from N4_modes_radial.validation import load_paper1_tables
    from R_paper_refs import paper1_radial_tables as P
    T = load_paper1_tables(R_TABLES)
    assert T['NUA'] == P.NUA and T['T1'] == P.T1 and T['T2'] == P.T2 and T['T3'] == P.T3


def test_winding_counts_the_zeros_inside_a_rectangle():
    from N4_modes_radial import modes as M
    f = lambda z: (z - (1 + 0.5j)) * (z - (2 - 1j)) * (z + 3 + 0j)
    rect = lambda x0, x1, y0, y1: [complex(x0, y0), complex(x1, y0), complex(x1, y1), complex(x0, y1)]
    assert M.winding_converged(f, rect(0, 3, -2, 1))[0] == 2
    assert M.winding_converged(f, rect(0.5, 1.5, 0, 1))[0] == 1
    assert M.winding_converged(f, rect(-5, 5, -3, 3))[0] == 3
    assert M.winding_converged(f, rect(6, 8, 1, 2))[0] == 0


def test_exterior_wave_does_not_depend_on_the_ray_and_is_analytic():
    """With the analytic exterior matrix of Eq. (3.7) (Schwarzschild + massive Klein-Gordon) the outgoing condition gives the same L on rotated rays and satisfies Cauchy-Riemann."""
    from N4_modes_radial import modes as M
    Mass, ma = 2.17, 0.12566
    f_ = lambda r: 1 - 2 * Mass / r

    def flat_e(r, w, *B):
        f = f_(r)
        return [0, 1, (ma**2 - w**2 / f) / f, -(2 / r) * (1 - Mass / r) / f]
    om = lambda nu: M.omega_of(nu)
    L = lambda nu, **kw: (lambda p: p[1] / p[0])(M.exterior_psi(flat_e, om(nu), Mass, 60.0, ma, 1.0, 1e-6, 0.88, rtol=1e-12, **kw))
    for nu in (8.1192 - 0.0036j, 6.585 - 0.563j, 5.9 + 0j):
        L0 = L(nu)
        assert abs(L(nu, dphi=0.4) / L0 - 1) < 1e-9 and abs(L(nu, dphi=-0.4) / L0 - 1) < 1e-9 and abs(L(nu, n_decay=60.0) / L0 - 1) < 1e-9
    z0, h = om(6.585 - 0.563j), om(1e-3)
    Lf = lambda o: (lambda p: p[1] / p[0])(M.exterior_psi(flat_e, o, Mass, 60.0, ma, 1.0, 1e-6, 0.88, rtol=1e-12))
    d1, d2 = (Lf(z0 + h) - Lf(z0)) / h, (Lf(z0 + 1j * h) - Lf(z0)) / (1j * h)
    assert abs(d1 - d2) / abs(d1) < 1e-2


def test_coefficient_tables_interpolate_a_matrix_with_1_over_r_entries(table_noaxion):
    """A = A_0(r) + w^2 A_1(r) with 1/r and 1/r^2 entries (like the regular singular centre): the r^2 A tabulation must keep A y accurate."""
    from N4_modes_radial.tables import RadialTables
    tab_star, _ = table_noaxion

    def entry(n, r, w):
        out = []
        for i in range(n):
            for j in range(n):
                out.append((1.0 / r**2 + (i + 1) * np.sin(r / 3) / r + 0.1 * j) * (1 if (i + j) % 2 else 0.5) + w**2 * (np.exp(-r / 5) / r + (j + 1) * 0.01 * i))
        return out
    fi = lambda r, w, *B: entry(4, r, w)
    fe = lambda r, w, *B: entry(2, r, w)
    tabs = RadialTables(tab_star, fi, fe)
    assert tabs.structure_error < 1e-12
    # matching at the surface (r_m <= R): no tail, the state is returned unchanged (control of step 21)
    from N4_modes_radial import modes as M
    y_R = np.array([1 + 0j, 2j])
    assert np.array_equal(M.outward(tab_star, tabs, 0.2, y_R, tab_star.R).y[:, -1], y_R)
    om = 0.2 - 0.03j
    f = tabs.rhs_int(om)
    rng = np.random.default_rng(1)
    for r in rng.uniform(1.2e-3, tab_star.R, 60):
        y = rng.normal(size=4) + 1j * rng.normal(size=4)
        exact = np.array(entry(4, r, om), dtype=complex).reshape(4, 4) @ y
        assert np.abs(f(r, y) - exact).max() <= 1e-7 * np.abs(exact).max(), r


def test_eos_cs2_definitions_and_the_thermodynamic_mismatch(table_noaxion):
    """The three c_s^2 definitions of step 23 on the S1 EOS table: 'mid' is the geometric mean, 'n' differs from 'eps' by ~1e-3 and the first law is violated by ~2e-4 at 100 MeV/fm^3."""
    from common import artifacts as A
    from common.star_table import _EOSVec
    cols, _ = A.load_table(table_noaxion[1] / 'eos_bsk26.csv')
    E = _EOSVec(cols['n_fm3'], cols['eps_MeV_fm3'], cols['p_MeV_fm3'])
    p = np.array([1.0, 10.0, 100.0, 300.0])
    c_eps, c_n = E._cs2_eps(p), E._cs2_n(p)
    E.cs2_mode = 'mid'
    assert np.allclose(E.cs2(p), np.sqrt(c_eps * c_n), rtol=1e-14)
    E.cs2_mode = 'n'
    assert np.array_equal(E.cs2(p), c_n)
    E.cs2_mode = 'eps'
    assert np.array_equal(E.cs2(p), c_eps)
    assert np.all(np.abs(c_n / c_eps - 1) < 5e-3)
    assert abs(E.mismatch(100.0)) < 1e-3 and abs(E.mismatch(100.0)) > 1e-5       # the published fits are consistent to ~1e-4, not exactly


def test_scattering_phase_ladder_recovers_width_and_position():
    """fit_ladder on a synthetic arg S = polynomial background + 2 arctan((x - nu_r)/g): the width and position come back."""
    from N4_modes_radial import analysis as AN
    nu0, g0 = 8.1148, 3.65e-6                      # kHz: resonance at nu0 with half-width g0
    nus = nu0 + g0 * np.tan(np.linspace(-1.47, 1.47, 61))
    x = (nus - nu0) * 1e3                          # Hz
    ph = 0.3 + 0.002 * x - 1e-5 * x**2 + 5e-8 * x**3 + 2 * np.arctan(x / (g0 * 1e3))
    tab, best, conv = AN.fit_ladder(nus, ph, nu0, g0)
    assert abs(best['half_width'] / g0 - 1) < 1e-6 and abs(best['nu_r'] - nu0) < 1e-9 * 1e3 and conv < 1e-3


def test_rydberg_step_on_an_exact_ladder(tmp_path):
    """step18 on a ladder built from kappa_n = n - delta with the relativistic Coulomb parameter: the Whittaker defect is constant and the two negative controls vary."""
    import json
    from N4_modes_radial import analysis as AN
    a, delta, nua = 0.2748, 0.0566, 6.0
    ns = np.arange(2, 7)
    kap = ns - delta
    q = (-kap + np.sqrt(kap**2 + 8 * a**2)) / (4 * a)
    nus = nua * np.sqrt(1 - q**2)
    (tmp_path / 'radial_star.json').write_text(json.dumps({'nu6': {'alpha_g': a, 'ladder': [{'n': int(n), 'nu': float(v)} for n, v in zip(ns, nus)]}}))
    res, out = AN.step18(tmp_path)
    assert [r[1] for r in res[:2]] == [True, True] and np.allclose(out['n_minus_kappa'], delta, atol=1e-12)
