"""Cache of the numerical functions that N2, N3 and N4 compile from the equation artifacts (decision 2026-11-21).

The symbolic reduction (solving the algebraic constraints, substituting the EOS and the background) and `lambdify(cse=True)` take
5-30 s per module and give the same functions on every run once the equations are settled. The cache stores the Python source that
lambdify generates, next to the module (`<module>/compiled/<name>.json`, tracked by git), and rebuilds the function from it with the same
namespace lambdify uses. Nothing is recompiled unless asked:

  * cache present and built from the same equations and code  -> loaded (milliseconds);
  * cache present but the equations or the reduction code changed -> StaleCache error: rerun the module with --recompile;
  * no cache yet -> built once and written;
  * `--recompile` (sets RECOMPILE) -> rebuilt and overwritten.

The key is the SHA-256 of the input expressions (srepr, without the per-process Dummy indices) and of the source files of the reduction.
"""
import hashlib
import json
import re
from pathlib import Path

import sympy as sp

RECOMPILE = False          # set by the modules' --recompile flag


class StaleCache(RuntimeError):
    pass


def expr_hash(eqs, keys=None):
    """SHA-256 of the expressions (dict key -> sympy expression), in key order."""
    h = hashlib.sha256()
    for k in sorted(keys if keys is not None else eqs, key=str):
        h.update(str(k).encode())
        h.update(re.sub(r'dummy_index=\d+', '', sp.srepr(eqs[k])).encode())      # Dummy indices are random per process
    return h.hexdigest()


def code_hash(*files):
    h = hashlib.sha256()
    for f in files:
        h.update(Path(f).name.encode())
        h.update(Path(f).read_bytes())
    return h.hexdigest()


_NS = {}


def _restore(source, module):
    if module not in _NS:
        _NS[module] = sp.lambdify((), 0, module).__globals__          # the namespace lambdify itself uses for this module
    ns = dict(_NS[module])
    exec(compile(source, '<compiled cache>', 'exec'), ns)              # noqa: S102 (our own generated code, tracked by git)
    return ns['_lambdifygenerated']


def compiled(path, key, build):
    """Functions {name: callable} from the cache at `path`, or from build() -> {name: (args, expression(s), lambdify module)}."""
    path = Path(path)
    if path.exists() and not RECOMPILE:
        data = json.loads(path.read_text(encoding='utf-8'))
        if data['key'] != key:
            changed = [k for k in key if data['key'].get(k) != key[k]]
            raise StaleCache(f'{path}: compiled from different {" and ".join(changed)}; rerun the module with --recompile to rebuild it')
        return {n: _restore(f['source'], f['module']) for n, f in data['functions'].items()}
    import inspect
    funcs, out = {}, {}
    for name, (args, expr, module) in build().items():
        fn = sp.lambdify(args, expr, module, cse=True)
        funcs[name] = fn
        out[name] = {'module': module, 'source': inspect.getsource(fn)}
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps({'key': key, 'sympy': sp.__version__, 'functions': out}, indent=1), encoding='utf-8')
    print(f'  compiled and cached: {path.name}', flush=True)
    return funcs
