"""S1 entry point:  python -m S1_star_background.run [--nu-a 2.8,3,...] [--witness] [--out DIR] [--tol 1e-6]

Tolerances (defaults, consolidated): background ODE rtol 1e-6, shooting precision 1e-6, Hermite table tolerance 1e-6,
about 4 s per star. The original precision (about 108 s per star) is
    python -m S1_star_background.run --rtol 1e-11 --shoot-tol 0 --tol 1e-9 --k0 8

Reads : nothing (BSk26 analytic fits are in bsk26.py).
Writes: outputs/eos_bsk26.csv(+.meta.json)          BSk26 EOS table
        outputs/star_noaxion.csv(+.meta.json)       reference star without axion (p_c = 100 MeV/fm^3)
        outputs/star_nua<nu>.csv(+.meta.json)       star with condensate for each nu_a (kHz); default: the 12 values of
                                                    Paper II Table II. --witness: only nu_a = 6.0
        outputs/manifest.json                       checks + extra.stars (M, R, theta_c, seconds), extra.condensate_exists
"""
import argparse
import time
from pathlib import Path

from common import artifacts as A
from common.checks import Timer, report

from . import star as S
from .export import write_eos, write_star

MODULE = 'S1_star_background'
HERE = Path(__file__).resolve().parent
NUS = [2.8, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 10.0, 11.0, 12.0, 13.0, 15.0]    # NUS of axion_ns/step9 (T1 u T2 keys)
REF = {'M': 1.53, 'R': 11.7}                                                 # Paper I: 1.53 M_sun, 11.7 km


def tag(nu):
    return f'star_nua{nu:g}'


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument('--nu-a', default=','.join(f'{x:g}' for x in NUS))
    ap.add_argument('--witness', action='store_true', help='only nu_a = 6.0 kHz (fast development point)')
    ap.add_argument('--out', default=str(HERE / 'outputs'))
    ap.add_argument('--tol', type=float, default=1e-6, help='Hermite table tolerance (default 1e-6)')
    ap.add_argument('--rtol', type=float, default=1e-6, help='background ODE rtol (default 1e-6; original code: 1e-11)')
    ap.add_argument('--shoot-tol', type=float, default=1e-6,
                    help='shooting precision in delta (default 1e-6; 0 = the original 64 bisections + refinement to 1e-16)')
    ap.add_argument('--max-step', type=float, default=0.25, help='max integrator step (km)')
    ap.add_argument('--k0', type=int, default=1, help='initial subdivisions per solver step in the table (doubles up to 128)')
    args = ap.parse_args(argv)
    out, tm = Path(args.out), Timer()
    kw = dict(rtol=args.rtol, shoot_tol=(None if args.shoot_tol <= 0 else args.shoot_tol), max_step=args.max_step)
    nus = [6.0] if args.witness else [float(x) for x in args.nu_a.split(',') if x]

    arts, res, stars = [], [], {}
    write_eos(out, MODULE)
    arts.append(out / 'eos_bsk26.csv')

    st0 = S.Star(None, **kw)
    meta0, info0 = write_star(out, st0, 'star_noaxion', MODULE, tol=args.tol, k0=args.k0)
    arts.append(out / 'star_noaxion.csv')
    stars['noaxion'] = {'M_msun': st0.M / S.MSUN, 'R_km': st0.R, **info0}
    res.append((f"Fondo BSk26 sin axión, p_c = 100 MeV/fm^3: M = {st0.M / S.MSUN:.3f} M_sun, R = {st0.R:.2f} km "
                f"(Paper I: 1.53 M_sun, 11.7 km)",
                abs(st0.M / S.MSUN - REF['M']) < 0.01 and abs(st0.R - REF['R']) < 0.05, ''))
    print(f'  noaxion: M={st0.M / S.MSUN:.4f} R={st0.R:.4f} rows={info0["n_rows"]} K={info0["K"]} '
          f'hermite_err={info0["hermite_max_rel_err"]:.1e}', flush=True)

    for nu in nus:
        t0 = time.time()
        st = S.Star(nu, **kw)
        meta, info = write_star(out, st, tag(nu), MODULE, tol=args.tol, k0=args.k0)
        arts.append(out / f'{tag(nu)}.csv')
        stars[f'{nu:g}'] = {'M_msun': st.M / S.MSUN, 'R_km': st.R, 'theta_c': st.thc, 'seconds': round(time.time() - t0, 1),
                            **info}
        print(f'  nu_a={nu:g}: M={st.M / S.MSUN:.4f} R={st.R:.4f} thc={st.thc:.6f} rows={info["n_rows"]} K={info["K"]} '
              f'hermite_err={info["hermite_max_rel_err"]:.1e} ({time.time() - t0:.0f} s)', flush=True)

    worst = max([v['hermite_max_rel_err'] for v in stars.values()])
    res.append((f'Tablas de perfil: el interpolante de Hermite reproduce la solución densa del integrador dentro de '
                f'{args.tol:g} (peor caso {worst:.1e}, {len(stars)} estrellas)', worst <= args.tol, ''))

    # condensate threshold (Paper I: it disappears below ~2.8 kHz)
    cond = {}
    try:
        S.Star(2.5, **kw)
        cond['2.5'] = True
    except ValueError:
        cond['2.5'] = False
    if 2.8 in nus:
        cond['2.8'] = True
    res.append(('Umbral del condensado: no hay a 2.5 kHz y sí a 2.8 kHz (Paper I: desaparece por debajo de ~2.8 kHz)',
                (cond['2.5'] is False and cond.get('2.8', False)) if '2.8' in cond else None,
                '' if '2.8' in cond else 'nu_a = 2.8 no calculado en esta ejecución (descriptivo)'))
    report('S1: fondo BSk26 con condensado (tablas)', res)
    A.write_manifest(out / 'manifest.json', MODULE, artifacts=arts, checks=res, seconds=tm.seconds(),
                     extra={'stars': stars, 'condensate_exists': cond, 'nu_a_list': nus, 'witness': args.witness})
    print(f'S1 done in {tm.seconds()} s', flush=True)
    return res


if __name__ == '__main__':
    main()
