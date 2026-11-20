"""Tables 1-3 of Paper I (Gomez-Banon, Pnigouras and Pons, PRD 113, 023053, arXiv:2510.16245v2, Sec. V): radial modes of the BSk26 star with the condensate
(eps = 0.1, p_0 = 100 MeV/fm^3, full GR, e^{-i omega t}), as typed in the paper (moved from axion_ns/step26_paper1_tables.py).

  T1: background with delta theta switched off, f0..f3 in kHz (None = no axion).
  T2: full problem, (Re in kHz, Im in Hz) of f0..f3.
  T3: axion modes, very damped, (Re in kHz, Im in kHz) of f1a and f2a (f0a is given with one digit and declared unreliable by the authors).
These are REFERENCES: nothing in a derivation module may import this file.
"""
NUA = (2.8, 3.0, 4.0, 5.0, 6.0, 8.0, 12.0)
T1 = {None: (3.10, 7.27, 10.22, 11.63), 2.8: (3.14, 7.32, 10.26, 11.67), 3.0: (3.32, 7.54, 10.46, 11.86),
      4.0: (3.63, 8.03, 11.25, 12.84), 5.0: (3.67, 8.14, 11.56, 13.45), 6.0: (3.68, 8.18, 11.71, 13.87),
      8.0: (3.69, 8.22, 11.88, 14.45), 12.0: (3.72, 8.27, 12.07, 15.16)}
T2 = {2.8: ((3.13, -8.80), (7.31, -2.60), (10.25, -1.04), (11.66, -0.394)),
      3.0: ((3.26, -43.7), (7.49, -6.10), (10.46, -2.15), (11.84, -0.854)),
      4.0: ((3.49, 0.0), (7.96, -38.9), (11.21, -4.55), (12.82, -1.06)),
      5.0: ((3.50, 0.0), (8.05, -8.38), (11.52, -7.74), (13.44, -4.02)),
      6.0: ((3.47, 0.0), (8.11, -2.05), (11.66, -22.5), (13.85, -2.01)),
      8.0: ((3.42, 0.0), (8.12, -36.2), (11.80, -6.71), (14.42, -5.22)),
      12.0: ((3.36, 0.0), (8.09, 0.0), (11.98, -12.0), (15.08, -11.4))}
T3 = {2.8: ((10.50, -4.30), (19.29, -5.82)), 3.0: ((11.23, -4.31), (20.06, -5.94)),
      4.0: ((16.36, -3.45), (22.99, -5.37)), 5.0: ((20.30, -2.22), (25.54, -4.70)),
      6.0: ((23.51, -1.50), (28.22, -3.81)), 8.0: ((29.86, -0.98), (33.55, -2.50)),
      12.0: ((42.97, -1.08), (45.19, -1.50))}


def as_json():
    key = lambda k: 'none' if k is None else f'{k:g}'
    return {'schema': 'axion-reftable/1', 'source': 'Paper I (arXiv:2510.16245v2), Tables 1-3', 'NUA': list(NUA),
            'units': {'T1': 'kHz (f0..f3)', 'T2': '[Re kHz, Im Hz] (f0..f3)', 'T3': '[Re kHz, Im kHz] (f1a, f2a)'},
            'T1': {key(k): list(v) for k, v in T1.items()}, 'T2': {key(k): [list(x) for x in v] for k, v in T2.items()},
            'T3': {key(k): [list(x) for x in v] for k, v in T3.items()}}
