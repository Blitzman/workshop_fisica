"""N4 against the committed legacy results (results/radial_*.json): the gate for consolidating N4.

    python -m tests.n4_gate [DIR]        DIR defaults to N4_modes_radial/outputs; files of steps that were not run are skipped

Numbers: relative difference below 1e-3 (the tolerance set for the N3 comparison) for frequencies and masses; integer counts must be equal. Check outcomes (ok flags) are compared separately,
by name with the digits masked, and every difference is listed (they do not decide PASS/FAIL: some legacy checks are known failures, see MIGRATION.md).
"""
import json
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
LEG = ROOT.parent / 'results'
NEW = Path(sys.argv[1]) if len(sys.argv) > 1 else ROOT / 'N4_modes_radial' / 'outputs'
TOL = 1e-3
rows, flags = [], []


def load(name):
    pl, pn = LEG / name, NEW / name
    if not pn.exists():
        return None, None
    return json.loads(pl.read_text(encoding='utf-8')), json.loads(pn.read_text(encoding='utf-8'))


def num(fn, q, a, b):
    rows.append((fn, q, 'num', a, b))


def cnt(fn, q, a, b):
    rows.append((fn, q, 'count', a, b))


def cplx(fn, q, a, b):          # [re, im] lists
    num(fn, q + ' Re', a[0], b[0])
    if a[1] or b[1]:
        num(fn, q + ' Im', a[1], b[1])


def norm(n):
    return re.sub(r'[-+]?\d[\d.eE+-]*', '#', n)[:70]


def compare_flags(fn, a, b):
    la = {norm(c['name']): c['ok'] for c in a['checks']}
    lb = {norm(c['name']): c['ok'] for c in b['checks']}
    for k in la:
        if k in lb and la[k] != lb[k]:
            flags.append((fn, k, la[k], lb[k]))
        if k not in lb:
            flags.append((fn, k, la[k], 'missing'))


def main():
    L, N = load('radial_paper1.json')
    if L:
        for k in L['t1']:
            for n, (x, y) in enumerate(zip(L['t1'][k], N['t1'][k])):
                if x and y:
                    num('paper1', f'T1 nu_a={k} f{n}', x, y)
        for k in L['t2']:
            for n, (x, y) in enumerate(zip(L['t2'][k], N['t2'][k])):
                if 'nu' in x and 'nu' in y:
                    num('paper1', f'T2 nu_a={k} f{n} Re', x['nu'][0], y['nu'][0])
                    if x['nu'][1]:
                        num('paper1', f'T2 nu_a={k} f{n} damping', x['nu'][1], y['nu'][1])
        compare_flags('radial_paper1', L, N)
    L, N = load('radial_star.json')
    if L:
        for q in ('p0_Mmax', 'p0_wF2_zero', 'Mmax_msun'):
            num('star', 'turning ' + q, L['turning'][q], N['turning'][q])
        num('star', 'F (6 kHz)', L['nu6']['F'], N['nu6']['F'])
        for a, b in zip(L['nu6']['ladder'], N['nu6']['ladder']):
            num('star', f'ladder n={a["n"]}', a['nu'], b['nu'])
        cplx('star', 'H1 (6 kHz)', L['nu6']['H1'], N['nu6']['H1'])
        for i, (x, y) in enumerate(zip(L['nu35']['bound'], N['nu35']['bound'])):
            num('star', f'3.5 kHz bound {i}', x, y)
        cplx('star', 'F (3.4 kHz)', L['nu34']['F'], N['nu34']['F'])
        compare_flags('radial_star', L, N)
    L, N = load('radial_crossing.json')
    if L:
        num('crossing', 'gap_min_khz', L['gap_min_khz'], N['gap_min_khz'])
        num('crossing', 'nu_a_min', L['nu_a_min'], N['nu_a_min'])
        rn = {round(r['nu_a'], 9): r for r in N['rows']}
        for r in L['rows']:
            m = rn.get(round(r['nu_a'], 9))
            if m:
                for i, (x, y) in enumerate(zip(r['roots'], m['roots'])):
                    num('crossing', f'nu_a={r["nu_a"]:.4f} root {i}', x, y)
        compare_flags('radial_crossing', L, N)
    L, N = load('radial_damped.json')
    if L:
        for k in ('nu6', 'thin', 'deep', 'old_lapse', 'old_lapse_deep'):
            cnt('damped', f'{k} index', L[k]['total'], N[k]['total'])
        for a, b in zip(sorted(L['nu6']['fluid_modes']), sorted(N['nu6']['fluid_modes'])):
            cplx('damped', 'fluid mode', a, b)
        compare_flags('radial_damped', L, N)
    L, N = load('radial_stability.json')
    if L:
        cnt('stability', 'control_axion index', L['control_axion']['total'], N['control_axion']['total'])
        cnt('stability', 'nu6 index', L['nu6.0']['total'], N['nu6.0']['total'])
        cnt('stability', 'nu3.4 index', L['nu3.4']['total'], N['nu3.4']['total'])
        for k in L['control_no_axion']:
            cnt('stability', f'no-axion control p_c={k}', L['control_no_axion'][k], N['control_no_axion'][k])
        compare_flags('radial_stability', L, N)
    L, N = load('radial_contour.json')
    if L:
        compare_flags('radial_contour', L, N)
    bad = 0
    worst = {}
    for fn, q, kind, a, b in rows:
        d = abs(b / a - 1) if kind == 'num' and a else (0.0 if a == b else float('inf'))
        w = worst.setdefault((fn, kind), (0.0, ''))
        if d > w[0]:
            worst[(fn, kind)] = (d, q)
        if d >= (TOL if kind == 'num' else 1e-12):
            bad += 1
            print(f'ABOVE TOLERANCE  {fn:10s} {q:34s} old {a!r:>22} new {b!r:>22}  rel diff {d:.1e}')
    print('\nworst difference by file (numbers: relative; counts: any difference):')
    for (fn, kind), (d, q) in worst.items():
        print(f'   {fn:10s} {kind:6s} {d:9.1e}  ({q})')
    print('\ncheck outcomes that differ between the legacy and N4 (name with digits masked):')
    for fn, k, a, b in flags:
        print(f'   {fn}: legacy {a} -> N4 {b}: {k}')
    if not flags:
        print('   none')
    print(f'\n{len(rows)} comparisons, {bad} above tolerance:', 'PASS' if not bad else 'FAIL')
    return not bad


if __name__ == '__main__':
    sys.exit(0 if main() else 1)
