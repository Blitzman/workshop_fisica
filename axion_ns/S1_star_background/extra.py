"""S1: the EXTRA stars that N4 (radial modes) needs besides the 12 of Paper II Table II: python -m S1_star_background.run --extra

  * nu_a = 3.4 and 3.5 kHz (steps 15, 16) and the nu_a values of the avoided-crossing study (step 19: its coarse grid and the six refinement points the legacy
    chose around the minimum of the gap at 3.5086 kHz; the committed legacy result lists them);
  * the star with condensate (nu_a = 6 kHz) at p_c = 1000 MeV/fm^3, beyond its maximum mass (step 16, positive control);
  * stars WITHOUT axion for p_c = 700..1100 MeV/fm^3 in steps of 10 (step 15: the turning point of the mass at p_c = 905, where omega_F^2 = 0; step 16: the controls
    at p_c = 800 and 1000). The legacy located the turning point with a bounded minimiser and a root finder on freshly integrated stars; here the scan is a fixed grid
    and N4 interpolates (it cannot call S1: modules exchange files).

Tables are written with fewer nodes than the main ones (they are only used by N4's checks) into outputs_extra/, which is NOT tracked by git (about 20 MB;
regenerate with `make S1-extra`, about 2 minutes). outputs_extra/manifest.json records their content hashes.
"""
import time
from pathlib import Path

from common import artifacts as A
from common.checks import Timer, report

from .export import write_eos, write_star
from .star_logp import SeedTracker, StarLogP

MODULE = 'S1_star_background'
NU_STEP15 = (3.4, 3.5)
NU_STEP19 = (3.46, 3.48, 3.4857142857142858, 3.4914285714285715, 3.4971428571428573, 3.5, 3.5028571428571427, 3.5085714285714285, 3.5142857142857142,
             3.52, 3.54, 3.57, 3.6, 3.65, 3.7)
P0_SCAN = tuple(float(x) for x in range(700, 1101, 10))
AXION_P0 = ((6.0, 1000.0),)
NODES_AXION = dict(n_s=2000, n_r=1000, n_geom=100, n_tail=500)
NODES_NOAXION = dict(n_s=2000, n_r=1000, n_geom=100, n_tail=60)


def name_of(nu, p0=None):
    if nu is None:
        return f'star_noaxion_p0_{p0:g}'
    return f'star_nua{nu:g}' if p0 is None else f'star_nua{nu:g}_p0_{p0:g}'


def all_stars():
    out = [(nu, None) for nu in sorted(set(NU_STEP15) | set(NU_STEP19))]
    out += [(nu, p0) for nu, p0 in AXION_P0]
    out += [(None, p0) for p0 in P0_SCAN]
    return out


def run_extra(out_dir, kw, tol):
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    (out_dir / '.gitignore').write_text('*\n!.gitignore\n')
    tm = Timer()
    arts, rows, bad = [], {}, []
    arts += [out_dir / 'eos_bsk26.csv']
    write_eos(out_dir, MODULE)
    seeds = SeedTracker()
    for nu, p0 in all_stars():
        t0 = time.time()
        seed = seeds.guess(nu) if (nu and p0 is None) else None
        st = StarLogP(nu, p0=p0, seed=seed, **kw)
        if nu and p0 is None:
            seeds.add(nu, st.u_c)
        name = name_of(nu, p0)
        meta, dense = write_star(out_dir, st, name, MODULE, tol=tol, **(NODES_AXION if nu else NODES_NOAXION))
        arts.append(out_dir / f'{name}.csv')
        rows[name] = {'nu_a': nu, 'p0': st.p0, 'M_msun': st.M_stop / 1.476625, 'R_km': st.R_stop, 'theta_c': st.thc, 'interp_err': dense['interp_max_err'],
                      'seconds': round(time.time() - t0, 1)}
        if dense['interp_max_err'] > tol:
            bad.append((name, dense['interp_max_err']))
        print(f'  {name}: M={rows[name]["M_msun"]:.4f} R={st.R_stop:.4f} interp_err={dense["interp_max_err"]:.1e} ({rows[name]["seconds"]} s)', flush=True)
    mass = [(r['p0'], r['M_msun']) for r in rows.values() if r['nu_a'] is None]
    pmax = max(mass, key=lambda t: t[1])
    res = [(f'{len(rows)} stars written (nu_a = 3.4-3.7 kHz, the p_c = 1000 star with condensate, and {len(mass)} stars without axion for p_c = 700-1100 MeV/fm^3)',
            len(rows) == len(all_stars()), ''),
           (f'Tablas de estado: la spline cúbica reproduce la solución dentro de {tol:g} en las {len(rows)} estrellas (peor caso '
            f'{max(r["interp_err"] for r in rows.values()):.1e})', not bad, '; '.join(f'{n}: {e:.1e}' for n, e in bad)),
           (f'La masa sin axión tiene su máximo en el interior de la malla: M = {pmax[1]:.4f} M_sun en p_c = {pmax[0]:g} (el original: 2.1695 M_sun en 904.96)',
            700 < pmax[0] < 1100, '')]
    report('S1 extra: estrellas para N4', res)
    A.write_manifest(out_dir / 'manifest.json', MODULE + '_extra', artifacts=arts, checks=res, seconds=tm.seconds(), extra={'stars': rows})
    return res
