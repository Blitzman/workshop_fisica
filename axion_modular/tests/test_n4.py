"""N4: fast tests (no symbolic compile): the Paper I tables artifact, the legacy winding function, the exterior wave on an analytic exterior matrix, and the coefficient tables with a
synthetic 1/r matrix. The comparison with the legacy results is `python -m tests.n4_gate`."""
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


def test_coefficient_tables_interpolate_a_matrix_with_1_over_r_entries(tmp_path):
    """A = A_0(r) + w^2 A_1(r) with 1/r and 1/r^2 entries (like the regular singular centre): the r^2 A tabulation must keep A y accurate."""
    from S1_star_background.export import write_eos, write_star
    from S1_star_background.star_logp import StarLogP
    from N4_modes_radial.star_table import StarTable
    from N4_modes_radial.tables import RadialTables
    st = StarLogP(None)
    write_eos(tmp_path, 'test')
    write_star(tmp_path, st, 'star_noaxion', 'test')
    tab_star = StarTable(tmp_path / 'eos_bsk26.csv', tmp_path / 'star_noaxion.csv')

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
    om = 0.2 - 0.03j
    f = tabs.rhs_int(om)
    rng = np.random.default_rng(1)
    for r in rng.uniform(1.2e-3, tab_star.R, 60):
        y = rng.normal(size=4) + 1j * rng.normal(size=4)
        exact = np.array(entry(4, r, om), dtype=complex).reshape(4, 4) @ y
        assert np.abs(f(r, y) - exact).max() <= 1e-7 * np.abs(exact).max(), r
