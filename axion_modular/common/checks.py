"""Check records shared by all modules (same semantics as axion_ns/pipeline.report).

A check is a tuple (name, ok, note):
  ok = True / False : a test that can fail;
  ok = None         : a descriptive result (reported, not counted).
"""
import time


def report(title, results):
    print(f'\n== {title}', flush=True)
    for name, ok, note in results:
        tag = 'RES ' if ok is None else ('OK ' if ok else 'FAIL')
        print(f'  [{tag}] {name}' + (f'  -- {note}' if note else ''), flush=True)
    return results


def ok_field(ok):
    return None if ok is None else bool(ok)


def to_records(results):
    return [{'name': n, 'ok': ok_field(ok), 'note': note or ''} for n, ok, note in results]


def summarize(results):
    counted = [c for c in results if c[1] is not None]
    return {'n_ok': sum(bool(c[1]) for c in counted), 'n_fail': sum(not c[1] for c in counted),
            'n_descriptive': len(results) - len(counted)}


class Timer:
    def __init__(self):
        self.t0 = time.time()

    def seconds(self):
        return round(time.time() - self.t0, 1)
