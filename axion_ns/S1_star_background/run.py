"""S1 entry point:  python -m S1_star_background.run [--nu-a 2.8,3,...] [--witness] [--out DIR] [--tol 1e-5]

The star is integrated against ln p (star_logp.py: start at r = 10 m from the centre series, stop at p = 1e-10 p_c, then the vacuum axion
tail in r). Tolerances (defaults): ODE rtol 1e-8, shooting precision 1e-12 (about 5 s per star; with 1e-6 / 1e-6 the valid axion tail is shorter and the damping of the p2 mode just above the threshold moves 4.7e-3); state table of about 7,000 nodes (uniform in ln p
+ uniform in r + geometric near the centre + clustered tail), interpolation checked at 1e-5. The original integration precision is
    python -m S1_star_background.run --rtol 1e-11 --shoot-tol 0

Reads : nothing (BSk26 analytic fits are in bsk26.py).
        --extra: only the extra stars of N4 (extra.py) into outputs_extra/ (untracked)
Writes: outputs/eos_bsk26.csv(+.meta.json)          BSk26 EOS table
        outputs/star_noaxion.csv(+.meta.json)       reference star without axion (p_c = 100 MeV/fm^3), dense state table
        outputs/star_nua<nu>.csv(+.meta.json)       star with condensate for each nu_a (kHz); default: the 12 values of
                                                    Paper II Table II. --witness: only nu_a = 6.0
        outputs/manifest.json                       checks + extra.stars (M, R, theta_c, interpolation errors, seconds),
                                                    extra.condensate_exists
"""
import argparse
import time
from pathlib import Path

from common import artifacts as A
from common.checks import Timer, report

from . import star as S
from .export import write_eos, write_star
from .star_logp import StarLogP

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
    ap.add_argument('--tol', type=float, default=1e-5, help='tolerance of the state-table interpolation check (default 1e-5)')
    ap.add_argument('--rtol', type=float, default=1e-8, help='background ODE rtol (default 1e-8; original code: 1e-11)')
    ap.add_argument('--shoot-tol', type=float, default=1e-12,
                    help='shooting precision in delta (default 1e-12; 0 = the original 64 bisections + refinement to 1e-16)')
    ap.add_argument('--max-step', type=float, default=0.25, help='max integrator step (km)')
    ap.add_argument('--extra', action='store_true', help='only the extra stars N4 needs (extra.py), into outputs_extra/ (not tracked by git)')
    ap.add_argument('--n-s', type=int, default=4000, help='table nodes uniform in ln p')
    ap.add_argument('--n-r', type=int, default=2000, help='table nodes uniform in r')
    ap.add_argument('--n-tail', type=int, default=1000, help='table nodes in the axion tail')
    args = ap.parse_args(argv)
    out, tm = Path(args.out), Timer()
    kw = dict(rtol=args.rtol, shoot_tol=(None if args.shoot_tol <= 0 else args.shoot_tol), max_step=args.max_step)
    nk = dict(n_s=args.n_s, n_r=args.n_r, n_tail=args.n_tail)
    if args.extra:
        from .extra import run_extra
        return run_extra(HERE / 'outputs_extra' if args.out == str(HERE / 'outputs') else args.out, kw, args.tol)
    nus = [6.0] if args.witness else [float(x) for x in args.nu_a.split(',') if x]

    arts, res, stars = [], [], {}
    write_eos(out, MODULE)
    arts.append(out / 'eos_bsk26.csv')

    t0 = time.time()
    st0 = StarLogP(None, **kw)
    meta0, d0 = write_star(out, st0, 'star_noaxion', MODULE, tol=args.tol, **nk)
    arts.append(out / 'star_noaxion.csv')
    stars['noaxion'] = {'M_msun': st0.M / S.MSUN, 'R_km': st0.R_stop, 'seconds': round(time.time() - t0, 1), **d0}
    res.append((f"Fondo BSk26 sin axión, p_c = 100 MeV/fm^3: M = {st0.M / S.MSUN:.3f} M_sun, R = {st0.R:.2f} km "
                f"(Paper I: 1.53 M_sun, 11.7 km)",
                abs(st0.M / S.MSUN - REF['M']) < 0.01 and abs(st0.R - REF['R']) < 0.05, ''))
    print(f'  noaxion: M={st0.M / S.MSUN:.4f} R={st0.R:.4f} rows={d0["n_rows"]} interp_err={d0["interp_max_err"]:.1e}', flush=True)

    for nu in nus:
        t0 = time.time()
        st = StarLogP(nu, **kw)
        meta, d = write_star(out, st, tag(nu), MODULE, tol=args.tol, **nk)
        arts.append(out / f'{tag(nu)}.csv')
        stars[f'{nu:g}'] = {'M_msun': st.M / S.MSUN, 'R_km': st.R_stop, 'theta_c': st.thc, 'seconds': round(time.time() - t0, 1), **d}
        print(f'  nu_a={nu:g}: M={st.M / S.MSUN:.4f} R={st.R:.4f} thc={st.thc:.6f} rows={d["n_rows"]} '
              f'interp_err={d["interp_max_err"]:.1e} ({time.time() - t0:.1f} s)', flush=True)

    worst = max(v['interp_max_err'] for v in stars.values())
    wname = max(stars, key=lambda k: stars[k]['interp_max_err'])
    res.append((f'Tablas de estado: la spline cúbica reproduce la solución del integrador dentro de {args.tol:g} '
                f'(peor caso {worst:.1e}, {len(stars)} estrellas, peor: {wname})', worst <= args.tol, ''))

    # condensate threshold (Paper I: it disappears below ~2.8 kHz)
    cond = {}
    try:
        StarLogP(2.5, **kw)
        cond['2.5'] = True
    except ValueError:
        cond['2.5'] = False
    if 2.8 in nus:
        cond['2.8'] = True
    res.append(('Umbral del condensado: no hay a 2.5 kHz y sí a 2.8 kHz (Paper I: desaparece por debajo de ~2.8 kHz)',
                (cond['2.5'] is False and cond.get('2.8', False)) if '2.8' in cond else None,
                '' if '2.8' in cond else 'nu_a = 2.8 no calculado en esta ejecución (descriptivo)'))
    report('S1: fondo BSk26 con condensado (tablas de estado)', res)
    A.write_manifest(out / 'manifest.json', MODULE, artifacts=arts, checks=res, seconds=tm.seconds(),
                     extra={'stars': stars, 'condensate_exists': cond, 'nu_a_list': nus, 'witness': args.witness})
    print(f'S1 done in {tm.seconds()} s', flush=True)
    return res


if __name__ == '__main__':
    main()
