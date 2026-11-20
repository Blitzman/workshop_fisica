"""Numbers published in Paper II (arXiv:2610.37308), Tables I-IV, as typed in the paper.

Moved verbatim from axion_ns/step6_nonradial_numerics.py (superset of the tables in step 9). Stored once here
and shipped as paper2_tables.ref.json for the numerical modules (N1 ...) that compare with them.
"""
# ---------------------------------------------------------------- tablas del paper
T1 = {None: (2.43, 7.15, 10.18), 2.8: (2.44, 7.19, 10.22), 3: (2.49, 7.35, 10.42), 4: (2.64, 7.80, 11.17),
      5: (2.70, 7.95, 11.49), 6: (2.73, 8.01, 11.64), 8: (2.78, 8.08, 11.82), 12: (2.84, 8.15, 12.01),
      15: (2.87, 8.19, 12.11)}
# nu_a: (f Re, f Im, p1 Re, p1 Im, p2 Re, p2 Im)  [Re en kHz, Im en Hz]
T2 = {2.8: (2.44, 0, 7.18, -0.28, 10.22, -0.74), 3: (2.47, 0, 7.32, -1.4, 10.40, -2.7),
      4: (2.55, 0, 7.72, -1.1, 11.14, -0.75), 5: (2.56, 0, 7.88, -0.074, 11.44, -12),
      6: (2.55, 0, 7.95, -1.4, 11.58, -13), 7: (2.54, 0, 7.98, -1.6, 11.68, -7.1),
      8: (2.54, 0, 7.99, 0, 11.75, -1.7), 10: (2.52, 0, 7.99, 0, 11.86, -1.5),
      11: (2.51, 0, 7.98, 0, 11.90, -3.8), 12: (2.51, 0, 7.98, 0, 11.91, 0),
      13: (2.50, 0, 7.97, 0, 11.95, 0), 15: (2.50, 0, 7.96, 0, 12.00, 0)}
# nu_a: (a0 Re, a0 Im, a1 Re, a1 Im, asterisco)
T3 = {2.8: (9.920, -2.623, 18.87, -5.49, False), 3: (10.15, -2.86, 19.56, -5.61, False),
      4: (10.2, -3.8, 22.42, -5.74, False), 5: (9.1, -4.2, 25.29, -4.72, False),
      6: (9.0, -3.8, 27.95, -3.61, False), 8: (9.3, -3.3, 33.18, -2.39, False),
      12: (10.9, -2.0, 44.68, -1.87, True), 15: (12.5, -0.5, 53.67, -2.13, True)}
# nu_a: (M/M_sun, alpha, nu_3, ratio_3, nu_4, ratio_4)
T4 = {8: (1.494, .370, 7.94, 1.06, 7.96, 1.05), 9: (1.501, .418, 8.91, 1.07, 8.95, 1.07),
      10: (1.507, .466, 9.87, 1.10, 9.93, 1.09), 11: (1.512, .515, 10.82, 1.12, 10.90, 1.12),
      12: (1.515, .563, 11.75, 1.15, 11.86, 1.15), 13: (1.519, .611, 12.68, 1.19, 12.82, 1.19),
      14: (1.521, .659, 13.58, 1.22, 13.77, 1.22), 15: (1.524, .707, 14.47, 1.25, 14.71, 1.24)}


COLUMNS = {
    'T1': ['f_kHz', 'p1_kHz', 'p2_kHz'],
    'T2': ['f_re_kHz', 'f_im_Hz', 'p1_re_kHz', 'p1_im_Hz', 'p2_re_kHz', 'p2_im_Hz'],
    'T3': ['a0_re', 'a0_im', 'a1_re', 'a1_im', 'asterisk'],
    'T4': ['M_msun', 'alpha', 'nu3_kHz', 'ratio3', 'nu4_kHz', 'ratio4'],
}
TITLES = {
    'T1': 'Table I: f, p1, p2 frequencies with the axion background only (delta theta off); key None = no axion',
    'T2': 'Table II: coupled modes, Re nu (kHz) and Im nu (Hz), l = 2',
    'T3': 'Table III: gravitational-atom quasi-bound state frequencies (a0, a1), kHz; asterisk if Re nu < nu_a',
    'T4': 'Table IV: star mass, alpha = G M m_a and bound levels n = 3, 4 (nu in kHz, ratio = -k^2 n^2/alpha^2)',
}


def as_json():
    out = {'schema': 'axion-reftable/1', 'source': 'arXiv:2610.37308', 'tables': {}}
    for key, tab in (('T1', T1), ('T2', T2), ('T3', T3), ('T4', T4)):
        rows = {('none' if k is None else str(k)): [(bool(x) if isinstance(x, bool) else float(x)) for x in v]
                for k, v in tab.items()}
        out['tables'][key] = {'title': TITLES[key], 'columns': COLUMNS[key], 'nu_a_kHz_rows': rows}
    return out
