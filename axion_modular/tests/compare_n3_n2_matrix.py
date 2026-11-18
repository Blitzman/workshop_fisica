"""SLOW, manual:  python -m tests.compare_n3_n2_matrix      (about 1 minute)

With the axion switched off, N3's numerical interior matrix must equal N2's (old step 10, matrix_checks): both are compiled from D2's artifact, N3's with the axion
background symbols and N2's with a concrete polytrope EOS. Cross-module on purpose (both modules are imported here and nowhere else).
"""
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
sys.path[:0] = [str(ROOT), str(ROOT.parent)]

from N2_modes_fullgr_noaxion import systems as S2      # noqa: E402
from N3_modes_fullgr_axion import systems as S3        # noqa: E402


def main():
    d2 = ROOT / 'D2_derive_nonradial' / 'outputs'
    ext, _ = S3.load_ext(d2)
    flat_i, _ = S3.compile_flat(ext)
    n2 = S2.full_system(ext, 'rest')                  # p = K rho^2 with rho the rest-mass density: eps = n + n^2
    rng, worst = np.random.default_rng(1), 0.0
    for _ in range(5):
        rr, om = rng.uniform(0.1, 1), rng.uniform(0.5, 2)
        n, n1, n2_ = rng.uniform(0.05, 0.3), -rng.uniform(0.01, 0.5), rng.uniform(-1, 1)
        P, P1, P2, Lm, L1 = rng.uniform(-1, 0), rng.uniform(0.01, 1), rng.uniform(-1, 1), rng.uniform(0, 0.5), rng.uniform(0.01, 1)
        ref = np.array(n2.mat(rr, om, 2, n, n1, n2_, P, P1, P2, Lm, L1), dtype=float)
        vals = flat_i(rr, om, 2, P, P1, P2, Lm, L1, n, n1, 0, 0, n + n * n, n * n, 2 * n / (1 + 2 * n), 1.0, 0, 0, 0, 0, 1.0)
        new = np.array([float(v) for v in vals]).reshape(6, 6)
        worst = max(worst, np.abs(new[:4, :4] - ref).max() / np.abs(ref).max())
    print(f'N3 interior matrix (axion off) vs N2 matrix at 5 random states: max relative difference {worst:.1e}  ->', 'OK' if worst < 1e-12 else 'FAIL')
    return worst < 1e-12


if __name__ == '__main__':
    sys.exit(0 if main() else 1)
