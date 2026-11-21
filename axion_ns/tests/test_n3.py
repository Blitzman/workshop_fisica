"""N3: fast tests of the table-based background pieces (no symbolic compile), on S1's 6 kHz star read through the shared reader."""


def test_analytic_phi_second_derivative_matches_a_finite_difference(table_6):
    from N3_modes_fullgr_axion.background import bg18
    tab, _ = table_6
    h = 1e-4
    dphi = lambda x: tab.rhs(x, list(tab.y_at(x)))[2]
    for r in (0.3, 2.0, 5.0, 9.0, 11.0, 11.3):
        m, p, Phi, t, tp = tab.y_at(r)
        an = bg18(tab, r, m, p, Phi, t, tp)[2][0]
        fd = (dphi(r + h) - dphi(r - h)) / (2 * h)
        assert abs(an / fd - 1) < 2e-5, (r, an, fd)


def test_bg18_is_consistent_with_the_star_reader(table_6):
    from N3_modes_fullgr_axion.background import bg18
    tab, _ = table_6
    for r in (0.5, 4.0, 8.0, 11.0):
        m, p, Phi, t, tp = tab.y_at(r)
        b = bg18(tab, r, m, p, Phi, t, tp)
        dm, dp, dPhi, dt, dtp = tab.rhs(r, [m, p, Phi, t, tp])
        assert abs(b[1][0] / dPhi - 1) < 1e-12 and abs(b[8][0] / tp - 1) < 1e-12
        bg = tab.bg(r)                                   # (e, p, n, cs2, q1, q2, dPhi, dlam, Phi + shift, lam, t, tp)
        assert abs(b[4][0] / bg[7] - 1) < 1e-10 and abs(b[3][0] / bg[9] - 1) < 1e-10
