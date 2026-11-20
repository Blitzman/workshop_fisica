"""N3: modes with the axion and the dynamic metric for the 12 values of nu_a of Paper II Table II (old step 10) and the convergence of the most delicate
ones (old step 25), table-based.

Port of axion_ns/step10_axion_gr.py and step25_exterior_convergence.py with the same algorithms, starting points and (where they exist) check texts. What changed:
the star comes from S1's table, the right-hand sides from the Laurent tables (tables.py: about 25x cheaper per root), the matching radius from an explicit rule
(modes.matching_radius), and the pieces that needed other packages are inputs: the Paper II tables from R, the no-axion Cowling modes and the Cowling shifts
from N1's table. NOT ported: the end-to-end LMI polytrope check of step 10 (it needs a polytrope EOS table in S1; N2 validates the same equations) and the
comparison of the numerical matrix with N2's (tests/compare_n3_n2_matrix.py, cross-module on purpose).
"""
import json
import multiprocessing as mp
import time

import numpy as np

from common.checks import ok_field

from . import modes as MD
from .star_table import StarTable
from .tables import LaurentTables

L = 2
MODES = ('f', 'p1', 'p2')
GR_GUESS = (1.96, 6.59, 9.87)           # kHz: f, p1, p2 of the no-axion star with the dynamic metric (starting points of step 10)
P2_COWLING_NOAXION = (2.43, 7.15, 10.18)   # Paper II, no-axion Cowling f, p1, p2 (kHz): reference of its tables I and II

_W = {}                                  # flat functions, tables of Paper II, paths: built BEFORE the fork


def setup(flats, T1p, T2p, eos_csv, s1_dir):
    _W.update(flat_i=flats[0], flat_e=flats[1], T1p=T1p, T2p=T2p, eos=eos_csv, s1=s1_dir)


def _star(nu):
    name = 'star_noaxion' if nu is None else f'star_nua{nu:g}'
    st = StarTable(_W['eos'], f"{_W['s1']}/{name}.csv")
    return st, LaurentTables(st, _W['flat_i'], _W['flat_e'])


def reference():
    """No-axion star: the three modes with the dynamic metric (the Cowling ones come from N1's table)."""
    st, tabs = _star(None)
    return [MD.find_mode(st, tabs, g, dyn=False) for g in GR_GUESS]


def row(args):
    nu, gr = args
    t0 = time.time()
    st, tabs = _star(nu)
    out = {'nu': nu, 'bg': [], 'cp': []}
    ks = sorted(_W['T1p'])
    for i in range(3):
        g0 = gr[i].real
        # starting points: the no-axion mode shifted as in Tables I and II of Paper II (Cowling)
        gb = g0 * np.interp(nu, ks, [_W['T1p'][k][i] for k in ks]) / P2_COWLING_NOAXION[i]
        gc = g0 * _W['T2p'][nu][i][0] / P2_COWLING_NOAXION[i]
        for key, guess, dyn in (('bg', gb, False), ('cp', gc, True)):
            try:
                out[key].append(MD.find_mode(st, tabs, guess, dyn=dyn))
            except Exception:  # noqa: BLE001
                out[key].append(complex(np.nan, np.nan))
    out['secs'] = time.time() - t0
    return out


def convergence(args):
    """The coupled f-mode at 6 kHz with another tolerance, another matching radius and another extraction distance."""
    nu, guess = args
    st, tabs = _star(nu)
    z0 = MD.find_mode(st, tabs, guess)
    rm = MD.matching_radius(st)
    alts = [MD.find_mode(st, tabs, z0.real, rtol=1e-11), MD.find_mode(st, tabs, z0.real, rm=0.8 * rm), MD.find_mode(st, tabs, z0.real, n_far=80)]
    return z0, alts


def _map(pool, fn, items):
    return pool.map(fn, items) if pool is not None else [fn(x) for x in items]


def compute(nus, procs, quick=False):
    ref_gr = reference()
    pool = mp.get_context('fork').Pool(procs) if (procs > 1 and 'fork' in mp.get_all_start_methods()) else None
    try:
        conv_job = None
        if not quick and 6.0 in nus:
            conv_job = (6.0, ref_gr[0].real * _W['T2p'][6.0][0][0] / P2_COWLING_NOAXION[0])
        res = pool.apply_async(convergence, (conv_job,)) if (pool is not None and conv_job) else None
        rows = {r['nu']: r for r in _map(pool, row, [(nu, ref_gr) for nu in sorted(nus, key=lambda x: -x)])}
        conv = res.get() if res is not None else (convergence(conv_job) if conv_job else None)
    finally:
        if pool is not None:
            pool.close()
            pool.join()
    return ref_gr, rows, conv


# ---------------------------------------------------------------- checks of step 10
def checks(ref_gr, rows, conv, subs, cow_n1=None):
    """cow_n1: dict with 'cow0' (no-axion Cowling f, p1, p2 in kHz), 'f_cow_max' (max Cowling f shift in %), 'p2_cow_im' ({nu_a: Cowling |Im| in Hz})."""
    res = [(f'Sustituciones del axión: {sum(ok for _, ok in subs)}/{len(subs)} cantidades del Paper II (n dm*/dtheta, m_eff^2, eps(n, theta), c_s^2...)',
            all(ok for _, ok in subs), '')]
    gr0 = np.array([z.real for z in ref_gr])
    if conv is not None:
        z0, alts = conv
        dre = max(abs(a.real / z0.real - 1) for a in alts)
        dim = max(abs(a.imag / z0.imag - 1) for a in alts)
        res.append((f'Convergencia (modo f, 6 kHz): tolerancia, radio de empalme y extracción mueven Re nu en {dre:.0e} e Im nu en {100 * dim:.1f} %',
                    dre < 1e-5 and dim < 0.01, ''))
    if cow_n1 is not None:
        cow0 = np.array(cow_n1['cow0'])
        res.append((f"Sin axión, la métrica dinámica baja las notas de Cowling: f, p1, p2 = {', '.join(f'{z.real:.3f}' for z in ref_gr)} kHz frente a "
                    f"{', '.join(f'{c:.3f}' for c in cow0)} ({', '.join(f'{100 * (g / c - 1):+.0f} %' for g, c in zip(gr0, cow0))}); "
                    f"f amortiguado por ondas gravitacionales en {1 / (2 * np.pi * -ref_gr[0].imag * 1e3):.3f} s", all(gr0 < cow0) and ref_gr[0].imag < 0, ''))
    cp = [z for r in rows.values() for z in r['cp']]
    res.append((f'Las {len(cp)} notas acopladas (f, p1, p2 para {len(rows)} masas del axión) encontradas y amortiguadas',
                all(np.isfinite(z.real) and z.imag < 0 for z in cp), ''))
    above = [(nu, i, z) for nu, r in rows.items() for i, z in enumerate(r['cp']) if i > 0]
    res.append(('Umbral: los modos p solo radian axiones por encima de nu_a (|Im nu| > 0.05 Hz si Re nu > nu_a; por debajo, solo ondas gravitacionales, < 0.05 Hz)',
                all((abs(z.imag) * 1e3 > 0.05) == (z.real > nu) for nu, i, z in above), ''))
    nus = sorted(rows)
    sf = [100 * (rows[nu]['cp'][0].real / gr0[0] - 1) for nu in nus]
    last = nus[-1]
    txt = (f'Con la métrica dinámica, el axión sube el modo f hasta un {max(sf):.1f} % '
           + (f'(Cowling, N1: {cow_n1["f_cow_max"]:.1f} %) ' if cow_n1 else '')
           + f'y su amortiguamiento por ondas gravitacionales pasa de {1 / (2 * np.pi * -ref_gr[0].imag * 1e3):.2f} a '
           f'{1 / (2 * np.pi * -rows[last]["cp"][0].imag * 1e3):.2f} s (a {last:g} kHz)')
    res.append((txt, None, 'resultado nuevo'))
    p2 = {nu: -rows[nu]['cp'][2].imag * 1e3 for nu in rows}
    num = max(p2, key=p2.get)
    res.append((f'Amortiguamiento máximo por axiones: p2 a {num:.0f} kHz, {p2[num]:.1f} Hz (tau_a = {1e3 / (2 * np.pi * p2[num]):.1f} ms)'
                + (f"; la suma 1/tau_GW + 1/tau_a con el tau_a de Cowling (p2 a 8 kHz: {p2.get(8.0, float('nan')):.2f} Hz frente a {cow_n1['p2_cow_im'][8.0]:.1f})"
                   if cow_n1 and 8.0 in p2 and 8.0 in cow_n1['p2_cow_im'] else ''), None, 'resultado nuevo'))
    return res


# ---------------------------------------------------------------- step 25: exterior branch and convergence of the delicate modes
BELOW = (2.102 - 0.001j, 2.4 - 0.01j, 5.0 - 0.1j, 5.9 - 0.5j)
ABOVE = (7.4 - 0.001j, 11.2 - 0.01j)
DELICATE = (('p1, justo por debajo del umbral (menor amortiguamiento)', 1, 8.0), ('p2, justo por debajo del umbral', 2, 12.0),
            ('p1, justo por encima del umbral', 1, 7.0), ('p2, justo por encima del umbral', 2, 11.0),
            ('p2, el que radia axiones con menor amortiguamiento', 2, 8.0), ('f, axión más pesado', 0, 15.0))
VARIANTS = {'base': {}, 'rtol': {'rtol': 1e-11}, 'rm': 'rm08', 'n_far': {'n_far': 80}}
NEG_MODE = 2
NEGS = {'rm=R': 'rmR', 'n_far=1': {'n_far': 1}}


def branch_test():
    st, _ = _star(6.0)
    rm = MD.matching_radius(st)
    Lf = lambda nu, **kw: complex(MD.exterior_L(st, L2, MD.omega_of(nu), r_stop=rm, **kw))
    rel = lambda a, b: abs(a / b - 1)
    out = {'below': [], 'above': []}
    for nu in BELOW:
        k = np.sqrt(complex(MD.omega_of(nu)**2 - st.ma**2))
        out['below'].append({'nu': [nu.real, nu.imag], 'k_principal': [k.real, k.imag], 'rel': rel(Lf(nu, branch='principal'), Lf(nu)),
                             'rel_short': rel(Lf(nu, branch='principal', n_decay=2), Lf(nu))})
    for nu in ABOVE:
        out['above'].append({'nu': [nu.real, nu.imag], 'rel': rel(Lf(nu, branch='principal'), Lf(nu)), 'rel_flipped': rel(Lf(nu, branch='flipped'), Lf(nu))})
    return out


L2 = 2


def _task(args):
    i, nu_a, kw, z0 = args
    st, tabs = _star(nu_a)
    if kw == 'rm08':
        kw = {'rm': 0.8 * MD.matching_radius(st)}
    elif kw == 'rmR':
        kw = {'rm': float(st.R)}
    try:
        return MD.find_mode(st, tabs, z0.real, dyn=True, **kw)
    except RuntimeError:
        return complex(np.nan, np.nan)


def step25(step10_rows, procs):
    """step10_rows: {nu_a: {'cp': [z_f, z_p1, z_p2]}} from step 10 (starting points, as the legacy read fig_axgr.json)."""
    res, out = [], {}
    B = branch_test()
    out['branch'] = B
    worst = max(b['rel'] for b in B['below'])
    res.append((f'H31: por debajo del umbral (nu_a = 6 kHz) la rama principal de k tiene Im k < 0 en {sum(b["k_principal"][1] < 0 for b in B["below"])}/{len(B["below"])} '
                f'frecuencias, pero L = psi\'/psi coincide con la de la rama que decae a {worst:.0e} relativo como mucho: '
                + ', '.join(f'{b["nu"][0]:.3g}{b["nu"][1]:+.3g}i kHz: {b["rel"]:.0e}' for b in B['below']), worst < 1e-10,
                'la integración hacia dentro (40 longitudes de decaimiento) filtra la solución creciente'))
    above = max(b['rel'] for b in B['above'])
    res.append((f'H31: por encima del umbral la rama no cambia (diferencia {above:.0e})', above < 1e-12, 'el original pedía == 0 (daba 0 exacto); aquí el redondeo deja ~1e-17'))
    short = min(b['rel_short'] for b in B['below'])
    res.append((f'Control negativo 1: arrancando a 2 longitudes de decaimiento, la rama principal cambia L al menos un {short:.0e} relativo (frente a {worst:.0e} con 40)',
                short > 1e-4, ''))
    flip = min(b['rel_flipped'] for b in B['above'])
    res.append((f'Control negativo 2: por encima del umbral, la onda entrante cambia L al menos un {flip:.1g} relativo', flip > 1e-2, ''))
    nv = len(VARIANTS)
    tasks = [(i, nu_a, kw, step10_rows[nu_a]['cp'][i]) for _, i, nu_a in DELICATE for kw in VARIANTS.values()]
    _, ni, nn = DELICATE[NEG_MODE]
    tasks += [(ni, nn, kw, step10_rows[nn]['cp'][ni]) for kw in NEGS.values()]
    pool = mp.get_context('fork').Pool(procs) if (procs > 1 and 'fork' in mp.get_all_start_methods()) else None
    try:
        Z = _map(pool, _task, tasks)
    finally:
        if pool is not None:
            pool.close()
            pool.join()
    rows = []
    for k, (name, i, nu_a) in enumerate(DELICATE):
        zs = dict(zip(VARIANTS, Z[k * nv:(k + 1) * nv]))
        z0 = zs['base']
        rows.append({'name': name, 'mode': ['f', 'p1', 'p2'][i], 'nu_a': nu_a, 'nu': [z0.real, z0.imag],
                     'dre_rel': {v: abs(zs[v].real / z0.real - 1) for v in VARIANTS if v != 'base'},
                     'dim_Hz': {v: abs(zs[v].imag - z0.imag) * 1e3 for v in VARIANTS if v != 'base'}})
    neg = dict(zip(NEGS, Z[len(DELICATE) * nv:]))
    out['delicate'] = rows
    res_im = lambda r_: max(0.01 * abs(r_['nu'][1]) * 1e3, 1e-9 * abs(complex(*r_['nu'])) * 1e3)
    finite = all(np.isfinite(z) for z in Z[:len(DELICATE) * nv])
    ok_re = finite and all(max(r_['dre_rel'].values()) < 1e-7 for r_ in rows)
    ok_im = finite and all(max(r_['dim_Hz'].values()) < res_im(r_) for r_ in rows)
    res.append((f'H33: en los {len(rows)} modos delicados del paso 10, tolerancia, radio de empalme y distancia de extracción mueven Re nu como mucho '
                f'{max(max(r_["dre_rel"].values()) for r_ in rows):.0e} relativo: ' + '; '.join(f'{r_["mode"]} a {r_["nu_a"]:g} kHz {max(r_["dre_rel"].values()):.0e}' for r_ in rows),
                ok_re, ''))
    res.append(('H33: e Im nu como mucho: ' + '; '.join(f'{r_["mode"]} a {r_["nu_a"]:g} kHz (Im {-r_["nu"][1] * 1e3:.3g} Hz) {max(r_["dim_Hz"].values()):.1e} Hz' for r_ in rows),
                ok_im, 'criterio: 1 % de Im nu o la resolución de Muller (1e-9 |nu|), lo que sea mayor'))
    r0_ = rows[NEG_MODE]
    z0 = complex(*r0_['nu'])
    dneg = {k: (abs(v.imag - z0.imag) * 1e3 if np.isfinite(v) else np.inf) for k, v in neg.items()}
    out['neg'] = {'mode': r0_['name'], 'nu_a': r0_['nu_a'], 'dim_Hz': dneg}
    base = max(r0_['dim_Hz'].values())
    res.append((f'Control negativo 3: empalmando en la superficie (sin la cola del campo), Im nu de {r0_["mode"]} a {r0_["nu_a"]:g} kHz se mueve {dneg["rm=R"]:.2g} Hz, '
                f'frente a {base:.1e} Hz con las variantes de arriba', dneg['rm=R'] > 10 * max(base, res_im(r0_)), ''))
    res.append((f'Control negativo 4: extrayendo la onda a una sola longitud de onda se mueve {dneg["n_far=1"]:.2g} Hz', dneg['n_far=1'] > 10 * max(base, res_im(r0_)), ''))
    out['checks'] = [{'name': n_, 'ok': ok_field(o_), 'note': t_} for n_, o_, t_ in res]
    return res, out
