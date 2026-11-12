"""SLy4 (Douchin & Haensel 2001): tabulated EoS and core composition, parsed from the authors' files.

- data/src/sly4.dat: n_b [fm^-3], rho [g/cm^3], P [dyn/cm^2] (http://www.ioffe.ru/astro/NSG/NSEOS/sly4.dat).
- paper/dh2001/dh2001.txt: `pdftotext -layout` of arXiv:astro-ph/0111092; Tables 4 and 6 are read from it.
"""
import re
from pathlib import Path

import numpy as np

HERE = Path(__file__).parent
DAT = HERE / 'data' / 'src' / 'sly4.dat'
TXT = HERE / 'paper' / 'dh2001' / 'dh2001.txt'
assert DAT.exists(), f'download http://www.ioffe.ru/astro/NSG/NSEOS/sly4.dat to {DAT}'
assert TXT.exists(), f'run pdftotext -layout on arXiv:astro-ph/0111092 into {TXT}'


def eos():
    rows = [ln.split() for ln in DAT.read_text().splitlines()]
    a = np.array([[float(x) for x in r[1:]] for r in rows if len(r) == 4 and r[0].isdigit()])
    return a[:, 0], a[:, 1], a[:, 2]


def composition():
    """Table 4: n_b [fm^-3] and the fractions x_p, x_e, x_mu (not in percent)."""
    lines = TXT.read_text().splitlines()
    i = next(k for k, ln in enumerate(lines) if ln.startswith('Table 4. Composition of the liquid core'))
    rows = []
    for ln in lines[i + 1:i + 40]:
        v = ln.split()
        if len(v) == 8 and all(re.fullmatch(r'\d+\.\d+', x) for x in v):
            rows += [v[:4], v[4:]]
    a = np.array(sorted([[float(x) for x in r] for r in rows]))
    assert len(a) == 38 and np.allclose(a[:, 1], a[:, 2] + a[:, 3], atol=2e-3), len(a)
    return a[:, 0], a[:, 1] / 100, a[:, 2] / 100, a[:, 3] / 100


def table6():
    """Maximum-mass configuration of the SLy row: M [Msun], R [km], n_c [fm^-3], rho_c [g/cm^3], P_c [dyn/cm^2]."""
    m = re.search(r'^\s+SLy\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)', TXT.read_text(), flags=re.M)
    M, R, nc, rc, Pc = map(float, m.groups())
    return {'m_max': M, 'r_max': R, 'n_c': nc, 'rho_c': rc * 1e15, 'p_c': Pc * 1e36}   # unit of rho_c: see notebook H2


def thresholds_quoted():
    t = TXT.read_text()
    e = float(re.search(r'direct electron Urca process, ([\d.]+) fm', t).group(1))
    mu = float(re.search(r'somewhat higher, ([\d.]+) fm', t).group(1))
    return e, mu
