"""Paper II tables (reference data from R's paper2_tables.ref.json) in the shapes the old steps 6 and 9 used."""
import json
from pathlib import Path


def load_paper2_tables(path):
    """Returns (T1, T2, T3, T4) as in axion_ns/step6: keys nu_a (float; None = no axion in T1), tuple values."""
    p = Path(path)
    if not p.exists():
        raise FileNotFoundError(f'{p} not found. Produce it first: python -m R_paper_refs.run')
    doc = json.loads(p.read_text(encoding='utf-8'))
    key = lambda k: None if k == 'none' else float(k)
    out = []
    for name in ('T1', 'T2', 'T3', 'T4'):
        rows = doc['tables'][name]['nu_a_kHz_rows']
        out.append({key(k): tuple(v) for k, v in rows.items()})
    return tuple(out)


def step9_shapes(T1, T2, T4):
    """The shapes used by axion_ns/step9: T1 {nu: (f, p1, p2)}, T2 {nu: ((re, im) x 3)}, T4M {nu: M}."""
    t1 = {k: v for k, v in T1.items() if k is not None}
    t2 = {k: ((v[0], v[1]), (v[2], v[3]), (v[4], v[5])) for k, v in T2.items()}
    t4m = {k: v[0] for k, v in T4.items()}
    return t1, t2, t4m
