"""Published mode frequencies of polytropes, as typed in the papers (moved from axion_ns/step8_numerics.py and fmodes.py).

  LMI97: Lindblom, Mendell and Ipser (1997, PRD 56, 2118; gr-qc/9704046), Table I, column "Exact": p = kappa rho^2 (rho: energy
         density), kappa = 1.3346e5 cm^5 g^-1 s^-2, 1.4 Msun, R = 14.45 km. {l: (Re w, Im w)} in units of sqrt(M/R^3).
  YK97:  Yoshida and Kojima (1997, MNRAS 289, 117; gr-qc/9705081), Tables 1-2: n = 1 polytrope, M/R = 0.1 and 0.2, modes f, p1, p2 with
         l = 2..5, {M/R: {l: {mode: (exact, Cowling)}}} in units of sqrt(M/R^3).
These are REFERENCES: nothing in a derivation module may import this file.
"""
LMI97 = {2: (1.2016, 3.95e-4), 3: (1.5862, 3.96e-6), 4: (1.8738, 4.01e-8), 5: (2.1129, 3.63e-10)}
LMI_R_KM = 14.45
LMI_KAPPA_CGS = 1.3346e5
LMI_MASS_MSUN = 1.4

YK97 = {0.1: {2: {'f': (1.21, 1.52), 'p1': (3.11, 3.43), 'p2': (4.83, 5.09)},
              3: {'f': (1.62, 1.82), 'p1': (3.61, 3.85), 'p2': (5.36, 5.58)},
              4: {'f': (1.93, 2.07), 'p1': (4.03, 4.22), 'p2': (5.84, 6.01)},
              5: {'f': (2.18, 2.29), 'p1': (4.40, 4.55), 'p2': (6.27, 6.42)}},
        0.2: {2: {'f': (1.17, 1.35), 'p1': (2.70, 2.99), 'p2': (4.12, 4.38)},
              3: {'f': (1.52, 1.64), 'p1': (3.19, 3.39), 'p2': (4.64, 4.84)},
              4: {'f': (1.78, 1.87), 'p1': (3.59, 3.74), 'p2': (5.10, 5.26)},
              5: {'f': (2.01, 2.08), 'p1': (3.94, 4.06), 'p2': (5.52, 5.65)}}}


def as_json():
    return {'schema': 'axion-reftable/1', 'source': 'LMI 1997 (gr-qc/9704046) Table I; YK 1997 (gr-qc/9705081) Tables 1-2',
            'units': 'sqrt(M/R^3)',
            'LMI97': {'star': {'kappa_cgs': LMI_KAPPA_CGS, 'mass_msun': LMI_MASS_MSUN, 'R_km': LMI_R_KM, 'eos': 'p = kappa rho^2, rho = energy density'},
                      'l': {str(k): {'re_w': v[0], 'im_w': v[1]} for k, v in LMI97.items()}},
            'YK97': {'eos': 'n = 1 polytrope', 'M_over_R': {str(mr): {str(ll): {m: {'exact': v[0], 'cowling': v[1]} for m, v in d.items()}
                                                                     for ll, d in t.items()} for mr, t in YK97.items()}}}
