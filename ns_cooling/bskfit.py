"""BSk22/24/25/26 unified equations of state from Potekhin's analytic fits.

Source: data/src/bskfit18.f (Pearson et al. 2018, MNRAS 481, 2994; http://www.ioffe.ru/astro/NSG/BSk/).
The coefficients are parsed from that file at import; only the formulas of BSKfit, BSkEofN and CORE18
are ported here. Units: n_b in fm^-3, rho in g/cm^3, P in dyn/cm^2.
"""
import hashlib
import re
from pathlib import Path

import numpy as np

SRC = Path(__file__).parent / 'data' / 'src' / 'bskfit18.f'
URL = 'http://www.ioffe.ru/astro/NSG/BSk/bskfit18.f'
SHA256 = '4bd3b716f04e40c69165fa83b5cd8ecb03c30aa1f2b45342fb070cef09a3a99c'
assert SRC.exists(), f'download {URL} to {SRC}'
assert hashlib.sha256(SRC.read_bytes()).hexdigest() == SHA256, 'bskfit18.f differs from the version used here (13.02.23)'
COL = {22: 0, 24: 1, 25: 2, 26: 3}          # column in the BSk22+ arrays
COL_P = {19: 0, 20: 1, 21: 2, 22: 3, 24: 4, 25: 5, 26: 6}


def _routine(name):
    text = SRC.read_text()
    m = re.search(rf'^ {{6}}subroutine {name}\(.*?^ {{6}}end\s*$', text, flags=re.S | re.M | re.I)
    lines = [ln for ln in m.group(0).splitlines() if ln[:1] not in '*cC']
    return ' '.join(ln[6:].split('!')[0] for ln in lines)


def _num(tok):
    return float(tok.lower().replace('d', 'e'))


def _data(name, array, rows):
    """Array `array` of a DATA statement in routine `name`, as (rows, ncols) -> [:, k]."""
    body = _routine(name)
    m = re.search(rf'\b{array}/([^/]*)/', body)
    vals = []
    for tok in m.group(1).replace(' ', '').split(','):
        if '*' in tok:
            n, v = tok.split('*')
            vals += [_num(v)] * int(n)
        elif tok:
            vals.append(_num(tok))
    return np.array(vals).reshape(-1, rows).T


def _const(name, par):
    return _num(re.search(rf'\b{par}=([-+.\dedED]+)', _routine(name)).group(1))


_AP = _data('BSKfit', 'A', 23)
_AE = _data('BSkEofN', 'A', 14)
_AYE = _data('CORE18', 'AYe', 6)
_XNC = _data('CRUCOR', 'AXNC', 4)[:, 0]
DMN_MEV = _const('BSkEofN', 'DMNMeV')
GROUND = _const('BSkEofN', 'GROUND')
FM3MEV = _const('BSkEofN', 'FM3MeV')
COMPTON = _const('CORE18', 'COMPTON')
RATMUE = _const('CORE18', 'RATMUE')
XPOW = {22: 3.674, 24: 4.0, 25: 4.0, 26: 4.6}   # exponents in the CORE18 if-chain
assert _AP.shape == (23, 7) and _AE.shape == (14, 4) and _AYE.shape == (6, 4) and len(_XNC) == 4
for _k, _v in XPOW.items():                      # the exponents above must be the ones in the file
    assert re.search(rf'KEOS\.eq\.{_k}\) then\s+XPOW=XN\*\*{_v:g}\b', _routine('CORE18')), _k


def _fermi(x):
    return 1.0 / (np.exp(np.clip(x, -40, 40)) + 1.0)


def lgP(keos, lgrho):
    a = np.concatenate([[0.0], _AP[:, COL_P[keos]]])
    x = lgrho
    return ((a[1] + a[2] * x + a[3] * x**3) / (1 + a[4] * x) * _fermi(a[5] * (x - a[6]))
            + (a[7] + a[8] * x) * _fermi(a[9] * (a[6] - x))
            + (a[10] + a[11] * x) * _fermi(a[12] * (a[13] - x))
            + (a[14] + a[15] * x) * _fermi(a[16] * (a[17] - x))
            + a[18] / (1 + ((x - a[19]) * a[20])**2) + a[21] / (1 + ((x - a[22]) * a[23])**2))


def rho_of_n(keos, n):
    a = np.concatenate([[0.0], _AE[:, COL[keos]]])
    EH = (a[10] * n)**a[13] / (1 + a[12] * n)
    WH = 1 / (1 + (a[11] * n)**a[14])
    EM = a[6] * n**a[7] * (1 + a[8] * n)
    WM = 1 / (1 + a[9] * n)
    EL = (a[1] * n)**1.16667 / (1 + np.sqrt(a[2] * n))
    EL = EL / (1 + np.sqrt(a[3] * n)) * (1 + np.sqrt(a[4] * n)) / (1 + np.sqrt(a[5] * n))
    return n * (DMN_MEV + EL * WM + EM * (1 - WM) * WH + EH * (1 - WH) + GROUND) * FM3MEV


def n_cc(keos):
    """Baryon density at the crust-core interface."""
    return _XNC[COL[keos]]


def fractions(keos, n):
    """Electron and muon fractions in the core (CORE18). Y_p = Y_e + Y_mu."""
    a = _AYE[:, COL[keos]]
    n = np.asarray(n, float)
    n32, xp = n * np.sqrt(n), n**XPOW[keos]
    ye = np.maximum((a[0] + a[1] * n + a[5] * n32 + a[2] * xp) / (1 + a[3] * n32 + a[4] * xp), 0.0)
    mu_e2 = 1 + (COMPTON * (3 * np.pi**2 * ye * n)**(1 / 3))**2
    xmu = np.sqrt(np.maximum(mu_e2 / RATMUE**2 - 1, 0.0))
    ymu = (RATMUE * xmu / COMPTON)**3 / (3 * np.pi**2 * n)
    return ye, ymu
