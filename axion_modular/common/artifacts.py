"""Artifact I/O: the only channel between modules.

Three kinds of file:

  *.eq.json      verified equations (SymPy expressions) with provenance and a content hash.
                 Expressions are stored as a de-duplicated node table (a DAG), which is compact
                 and can be decoded without recursion. `text` fields give a readable form; the
                 typeset form lives in the companion PDF (common/pdf.py).
  *.csv + .meta.json   numerical tables: plain CSV (17 significant digits) plus a sidecar with
                 units, parameters, producer and the sha256 of the CSV.
  manifest.json  per-module summary: artifacts with hashes, checks, runtime, environment.

Every artifact records the sha256 of the artifacts it was computed from (`inputs`), so a
consumer can detect that an upstream file changed after the artifact was produced.
"""
from __future__ import annotations

import csv
import datetime as _dt
import hashlib
import json
import platform
import subprocess
import sys
from pathlib import Path

import numpy as np
import sympy as sp

SCHEMA_EQ = 'axion-eq/1'
SCHEMA_TABLE = 'axion-table/1'
SCHEMA_MANIFEST = 'axion-manifest/1'
ROOT = Path(__file__).resolve().parent.parent

sys.setrecursionlimit(max(sys.getrecursionlimit(), 20000))


# ------------------------------------------------------------------ small helpers
def sha256_bytes(b: bytes) -> str:
    return hashlib.sha256(b).hexdigest()


def sha256_file(path) -> str:
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for chunk in iter(lambda: f.read(1 << 20), b''):
            h.update(chunk)
    return h.hexdigest()


def env_info() -> dict:
    out = {'python': platform.python_version(), 'platform': f'{platform.system()} {platform.machine()}'}
    for pkg in ('sympy', 'mpmath', 'numpy', 'scipy', 'matplotlib'):
        try:
            out[pkg] = __import__(pkg).__version__
        except Exception:  # noqa: BLE001
            pass
    return out


def git_info() -> dict:
    try:
        run = lambda *a: subprocess.run(['git', *a], cwd=ROOT, capture_output=True, text=True, timeout=10).stdout.strip()
        return {'commit': run('rev-parse', '--short', 'HEAD'), 'dirty': bool(run('status', '--porcelain', '--', str(ROOT)))}
    except Exception:  # noqa: BLE001
        return {'commit': None, 'dirty': None}


def now_iso() -> str:
    return _dt.datetime.now(_dt.timezone.utc).strftime('%Y-%m-%dT%H:%M:%SZ')


def producer_block(module: str, script: str | None = None) -> dict:
    return {'module': module, 'script': script, 'git': git_info(), 'created': now_iso(), 'env': env_info()}


def fingerprint(path) -> str:
    """Stable identity of an artifact: the content hash for equation artifacts (independent of timestamps and
    provenance fields), the sha256 of the bytes for everything else (CSV, JSON data, manifests)."""
    p = Path(path)
    if p.name.endswith('.eq.json'):
        return json.loads(p.read_text(encoding='utf-8'))['content_sha256']
    return sha256_file(p)


def input_record(path, name: str | None = None) -> dict:
    """Record of an upstream artifact: name, path (relative to the project root when possible) and its fingerprint
    (stored under the key `sha256`; for equation artifacts it is the content hash, so reruns reproduce it)."""
    p = Path(path).resolve()
    try:
        shown = str(p.relative_to(ROOT))
    except ValueError:
        shown = str(p)
    return {'name': name or p.name, 'path': shown, 'sha256': fingerprint(p)}


# ------------------------------------------------------------------ expression <-> node table
_CONST = {'pi': sp.pi, 'E': sp.E, 'I': sp.I, 'oo': sp.oo, '-oo': -sp.oo, 'zoo': sp.zoo, 'nan': sp.nan}


def _assumptions(s) -> dict:
    a = getattr(s, '_assumptions_orig', None)
    if a is None:
        a = {k: v for k, v in s.assumptions0.items() if k in ('real', 'positive', 'negative', 'integer', 'nonzero',
                                                               'nonnegative', 'imaginary', 'complex', 'finite')}
    return {k: bool(v) for k, v in sorted(a.items())}


class _Encoder:
    def __init__(self):
        self.index: dict = {}
        self.nodes: list = []
        self.dummy_ord: dict = {}       # Dummy -> ordinal of first appearance (the index of a Dummy is random per process)

    def visit(self, e) -> int:
        i = self.index.get(e)
        if i is not None:
            return i
        if isinstance(e, sp.Dummy):
            node = ['d', e.name, self.dummy_ord.setdefault(e, len(self.dummy_ord)), _assumptions(e)]
        elif isinstance(e, sp.Symbol):
            node = ['s', e.name, _assumptions(e)]
        elif isinstance(e, sp.Integer):
            node = ['n', int(e)]
        elif isinstance(e, sp.Rational):
            node = ['q', int(e.p), int(e.q)]
        elif isinstance(e, sp.Float):
            node = ['x', str(e), int(e._prec)]
        elif e is sp.pi:
            node = ['c', 'pi']
        elif e is sp.E:
            node = ['c', 'E']
        elif e is sp.I:
            node = ['c', 'I']
        elif e is sp.oo:
            node = ['c', 'oo']
        elif e is -sp.oo:
            node = ['c', '-oo']
        elif e is sp.zoo:
            node = ['c', 'zoo']
        elif e is sp.nan:
            node = ['c', 'nan']
        elif isinstance(e, sp.core.function.AppliedUndef):
            node = ['F', e.func.__name__, [self.visit(a) for a in e.args]]
        elif isinstance(e, sp.Derivative):
            node = ['D', self.visit(e.expr), [[self.visit(v), int(c)] for v, c in e.variable_count]]
        elif isinstance(e, sp.Subs):
            node = ['S', self.visit(e.expr), [self.visit(v) for v in e.variables], [self.visit(p) for p in e.point]]
        else:
            name = type(e).__name__
            if getattr(sp, name, None) is not type(e):
                raise TypeError(f'cannot serialize {name}: {str(e)[:80]}')
            node = ['f', name, [self.visit(a) for a in e.args]]
        self.nodes.append(node)
        self.index[e] = len(self.nodes) - 1
        return self.index[e]


def encode_exprs(exprs: dict) -> dict:
    enc = _Encoder()
    roots = {k: enc.visit(sp.sympify(v)) for k, v in exprs.items()}
    return {'nodes': enc.nodes, 'roots': roots}


def _children(node):
    kind = node[0]
    if kind in ('F', 'f'):
        return list(node[2])
    if kind == 'D':
        return [node[1]] + [v for v, _ in node[2]]
    if kind == 'S':
        return [node[1]] + list(node[2]) + list(node[3])
    return []


def decode_exprs(dag: dict, keys=None) -> dict:
    nodes, roots = dag['nodes'], dag['roots']
    keys = list(roots) if keys is None else list(keys)
    need, stack = set(), [roots[k] for k in keys]
    while stack:                                    # reachable nodes only
        i = stack.pop()
        if i in need:
            continue
        need.add(i)
        stack.extend(_children(nodes[i]))
    objs, dummies = {}, {}
    for i in sorted(need):                          # children always precede parents
        nd = nodes[i]
        kind = nd[0]
        if kind == 's':
            objs[i] = sp.Symbol(nd[1], **nd[2])
        elif kind == 'd':
            if nd[2] not in dummies:
                dummies[nd[2]] = sp.Dummy(nd[1], **nd[3])
            objs[i] = dummies[nd[2]]
        elif kind == 'n':
            objs[i] = sp.Integer(nd[1])
        elif kind == 'q':
            objs[i] = sp.Rational(nd[1], nd[2])
        elif kind == 'x':
            objs[i] = sp.Float(nd[1], precision=nd[2])
        elif kind == 'c':
            objs[i] = _CONST[nd[1]]
        elif kind == 'F':
            objs[i] = sp.Function(nd[1])(*[objs[j] for j in nd[2]])
        elif kind == 'D':
            objs[i] = sp.Derivative(objs[nd[1]], *[(objs[v], c) for v, c in nd[2]])
        elif kind == 'S':
            objs[i] = sp.Subs(objs[nd[1]], tuple(objs[j] for j in nd[2]), tuple(objs[j] for j in nd[3]))
        elif kind == 'f':
            objs[i] = getattr(sp, nd[1])(*[objs[j] for j in nd[2]])
        else:
            raise ValueError(f'unknown node kind {kind!r}')
    return {k: objs[roots[k]] for k in keys}


def _same(a, b) -> bool:
    if a == b:
        return True
    return sp.expand(a - b) == 0


# ------------------------------------------------------------------ equation artifacts
class EqArtifact:
    """Loaded equation artifact: `.eq` maps key -> SymPy expression; `.entries` holds label/ref/note."""

    def __init__(self, path, meta, eq):
        self.path, self.meta, self.eq = Path(path), meta, eq
        self.entries = meta['equations']
        self.sha256 = meta['content_sha256']
        self.name = meta['name']

    def __getitem__(self, k):
        return self.eq[k]

    def keys(self):
        return self.eq.keys()


def _content_hash(name, kind, entries, dag, inputs):
    payload = {'schema': SCHEMA_EQ, 'name': name, 'kind': kind, 'dag': dag,
               'equations': {k: {'root': v['root'], 'label': v.get('label', ''), 'form': v.get('form', '')}
                             for k, v in entries.items()},
               'inputs': [{'name': i['name'], 'sha256': i['sha256']} for i in inputs]}
    return sha256_bytes(json.dumps(payload, sort_keys=True, separators=(',', ':')).encode())


def save_equations(path, name, entries: dict, *, title='', kind='derived', description='', module='?', script=None,
                   inputs=(), checks=None, verify=True, text_limit=4000):
    """Write an equation artifact.

    entries: key -> expression, or key -> dict(expr=..., label=..., ref=..., note=..., form=...).
    kind: 'derived' (from the action), 'transcribed' (typed from a paper), 'candidate' (hand assembled,
    to be verified). Returns the metadata dict (including content_sha256).
    """
    path = Path(path)
    norm = {k: (dict(v) if isinstance(v, dict) and 'expr' in v else {'expr': v}) for k, v in entries.items()}
    dag = encode_exprs({k: v['expr'] for k, v in norm.items()})
    if verify:
        back = decode_exprs(dag)
        bad = [k for k, v in norm.items() if not _same(back[k], sp.sympify(v['expr']))]
        if bad:
            raise ValueError(f'round-trip mismatch in {name}: {bad}')
    eqs = {}
    for k, v in norm.items():
        text = sp.sstr(v['expr'])
        eqs[k] = {'label': v.get('label', ''), 'ref': v.get('ref', ''), 'note': v.get('note', ''),
                  'form': v.get('form', 'expr = 0'), 'root': dag['roots'][k],
                  'n_ops': int(sp.count_ops(v['expr'])),
                  'text': text if len(text) <= text_limit else text[:text_limit] + ' ...[truncated]'}
    inputs = list(inputs)
    meta = {'schema': SCHEMA_EQ, 'name': name, 'title': title or name, 'kind': kind, 'description': description,
            'producer': producer_block(module, script), 'inputs': inputs, 'equations': eqs,
            'content_sha256': _content_hash(name, kind, eqs, dag, inputs)}
    if checks is not None:
        meta['checks'] = checks
    doc = dict(meta)
    doc['dag'] = dag
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(doc, ensure_ascii=False, separators=(',', ':')), encoding='utf-8')
    meta['file_sha256'] = sha256_file(path)
    meta['bytes'] = path.stat().st_size
    return meta


def load_equations(path, keys=None, expect_name=None, expect_sha=None) -> EqArtifact:
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f'{path} not found. Produce it first (see that module\'s README / `make`).')
    doc = json.loads(path.read_text(encoding='utf-8'))
    if doc.get('schema') != SCHEMA_EQ:
        raise ValueError(f'{path}: unexpected schema {doc.get("schema")!r}')
    if expect_name and doc['name'] != expect_name:
        raise ValueError(f'{path}: expected artifact {expect_name!r}, found {doc["name"]!r}')
    dag = {'nodes': doc['dag']['nodes'], 'roots': {k: v['root'] for k, v in doc['equations'].items()}}
    recomputed = _content_hash(doc['name'], doc['kind'], doc['equations'], doc['dag'], doc['inputs'])
    if recomputed != doc['content_sha256']:
        raise ValueError(f'{path}: content hash mismatch (file edited or corrupted)')
    if expect_sha and expect_sha != doc['content_sha256']:
        raise ValueError(f'{path}: content hash differs from the one recorded by the consumer')
    meta = {k: v for k, v in doc.items() if k != 'dag'}
    return EqArtifact(path, meta, decode_exprs(dag, keys))


# ------------------------------------------------------------------ table artifacts
def save_table(csv_path, columns: dict, *, name, units=None, descriptions=None, parameters=None, module='?',
               script=None, inputs=(), title=''):
    """Write columns (equal-length 1-D arrays; numeric or str) as CSV + sidecar <stem>.meta.json."""
    csv_path = Path(csv_path)
    csv_path.parent.mkdir(parents=True, exist_ok=True)
    names = list(columns)
    cols = [np.asarray(columns[c]) for c in names]
    n = len(cols[0])
    if any(len(c) != n for c in cols):
        raise ValueError('columns of different length')
    with open(csv_path, 'w', newline='', encoding='utf-8') as f:
        wr = csv.writer(f, lineterminator='\n')
        wr.writerow(names)
        for i in range(n):
            wr.writerow([repr(float(c[i])) if c.dtype.kind in 'fiub' else str(c[i]) for c in cols])
    meta = {'schema': SCHEMA_TABLE, 'name': name, 'title': title or name, 'n_rows': n,
            'columns': [{'name': c, 'dtype': str(x.dtype), 'unit': (units or {}).get(c, ''),
                         'description': (descriptions or {}).get(c, '')} for c, x in zip(names, cols)],
            'parameters': parameters or {}, 'producer': producer_block(module, script), 'inputs': list(inputs),
            'csv_sha256': sha256_file(csv_path)}
    csv_path.with_suffix('.meta.json').write_text(json.dumps(meta, indent=1, ensure_ascii=False), encoding='utf-8')
    return meta


def load_table(csv_path, verify=True):
    """Returns (columns: dict name -> np.ndarray, meta: dict)."""
    csv_path = Path(csv_path)
    if not csv_path.exists():
        raise FileNotFoundError(f'{csv_path} not found. Produce it first (see that module\'s README / `make`).')
    meta = json.loads(csv_path.with_suffix('.meta.json').read_text(encoding='utf-8'))
    if meta.get('schema') != SCHEMA_TABLE:
        raise ValueError(f'{csv_path}: unexpected schema {meta.get("schema")!r}')
    if verify and sha256_file(csv_path) != meta['csv_sha256']:
        raise ValueError(f'{csv_path}: sha256 differs from its .meta.json (file edited or corrupted)')
    with open(csv_path, newline='', encoding='utf-8') as f:
        rows = list(csv.reader(f))
    head, body = rows[0], rows[1:]
    cols = {}
    for j, c in enumerate(head):
        dt = next(x['dtype'] for x in meta['columns'] if x['name'] == c)
        vals = [row[j] for row in body]
        cols[c] = np.array(vals, dtype=float) if dt[0] in 'fiub' else np.array(vals, dtype=str)
    return cols, meta


# ------------------------------------------------------------------ manifest
def write_manifest(path, module, *, artifacts=(), checks=None, seconds=None, inputs=(), extra=None):
    """Per-module summary. `artifacts` are file paths (hashed here); `checks` is a list of (name, ok, note)."""
    from .checks import summarize, to_records
    path = Path(path)
    arts = [{'path': str(Path(a).relative_to(path.parent)) if Path(a).is_absolute() and path.parent in Path(a).parents
             else str(a), 'sha256': sha256_file(a), 'bytes': Path(a).stat().st_size} for a in artifacts]
    doc = {'schema': SCHEMA_MANIFEST, 'module': module, 'producer': producer_block(module), 'seconds': seconds,
           'inputs': list(inputs), 'artifacts': arts}
    if checks is not None:
        doc['summary'] = summarize(checks)
        doc['checks'] = to_records(checks)
    if extra:
        doc['extra'] = extra
    path.write_text(json.dumps(doc, indent=1, ensure_ascii=False), encoding='utf-8')
    return doc
