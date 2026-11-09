"""Paso 19: ¿cruce evitado entre el modo F y el nivel n = 1 del átomo gravitatorio? (hipótesis H21).

Para cada nu_a de una malla, los dos estados ligados más bajos de la estrella realista (raíces reales de
radial.det_entire bajo el umbral, que es real y sin polos en el eje real con la rama decreciente) y la fracción
fluida de cada uno. Cruce evitado = separación con un mínimo > 0 que no baja al refinar + intercambio de carácter.
"""
import json
import multiprocessing as mp
from pathlib import Path

import numpy as np
from scipy.optimize import brentq

from . import axstar as AX
from . import radial as RD
from .pipeline import ok_field, report

OUT = Path(__file__).resolve().parent.parent / 'results'
NU_A = [3.46, 3.48, 3.50, 3.52, 3.54, 3.57, 3.60, 3.65, 3.70]


def two_lowest(nu_a, lo=3.25, step=0.002):
    st = AX.Star(nu_a)
    f = lambda nu: RD.det_entire(st, AX.omega_of(nu)).real
    grid = np.arange(lo, nu_a - 0.004, step)
    v = [f(x) for x in grid]
    roots = []
    for i in range(len(grid) - 1):
        if v[i] * v[i + 1] < 0:
            roots.append(brentq(f, grid[i], grid[i + 1], xtol=1e-9))
            if len(roots) == 2:
                break
    fr = [RD.fluid_fraction(st, AX.omega_of(r)) for r in roots]
    return {'nu_a': nu_a, 'roots': roots, 'fluid': fr, 'M': st.M}


def run(nus=NU_A):
    with mp.get_context('fork').Pool(6) as pool:
        rows = sorted(pool.map(two_lowest, nus), key=lambda r: r['nu_a'])
    ok_rows = [r for r in rows if len(r['roots']) == 2]
    gap = [(r['nu_a'], r['roots'][1] - r['roots'][0]) for r in ok_rows]
    i = int(np.argmin([g for _, g in gap]))
    # refinamiento alrededor del mínimo
    if 0 < i < len(gap) - 1:
        fine = list(np.linspace(gap[i - 1][0], gap[i + 1][0], 8)[1:-1])
        with mp.get_context('fork').Pool(6) as pool:
            rows = sorted(rows + pool.map(two_lowest, fine), key=lambda r: r['nu_a'])
        ok_rows = [r for r in rows if len(r['roots']) == 2]
        gap = [(r['nu_a'], r['roots'][1] - r['roots'][0]) for r in ok_rows]
        i = int(np.argmin([g for _, g in gap]))
    nu_min, g_min = gap[i]
    interior = 0 < i < len(gap) - 1
    lo_side, hi_side = ok_rows[0], ok_rows[-1]
    swap = (lo_side['fluid'][0] < 0.5 < lo_side['fluid'][1]) and (hi_side['fluid'][0] > 0.5 > hi_side['fluid'][1])
    res = [
        (f'H21(a): separación de los dos ligados más bajos con mínimo {1e3 * g_min:.1f} Hz en nu_a = {nu_min:.3f} kHz '
         f'(interior a la malla, refinada; {len(gap)} valores de nu_a)', interior and g_min > 0, ''),
        (f'H21(b): intercambio de carácter: a nu_a = {lo_side["nu_a"]:.2f} kHz fracción fluida {lo_side["fluid"][0]:.2f} '
         f'(inferior) y {lo_side["fluid"][1]:.2f} (superior); a {hi_side["nu_a"]:.2f} kHz, {hi_side["fluid"][0]:.2f} y '
         f'{hi_side["fluid"][1]:.2f}', swap, ''),
    ]
    out = {'rows': rows, 'gap_min_khz': g_min, 'nu_a_min': nu_min,
           'checks': [{'name': n_, 'ok': ok_field(o_), 'note': t_} for n_, o_, t_ in res]}
    OUT.mkdir(exist_ok=True)
    (OUT / 'radial_crossing.json').write_text(json.dumps(out, indent=1))
    return report('Paso 19: cruce evitado F - átomo gravitatorio', res)


if __name__ == '__main__':
    run()
