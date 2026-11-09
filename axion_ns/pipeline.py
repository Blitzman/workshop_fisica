"""Derivaciones compartidas por los pasos (cacheadas)."""
from functools import lru_cache

from .compare import Rules
from .model import derive_field_equations, background, linearize


@lru_cache
def field_equations():
    return derive_field_equations()


@lru_cache
def background_equations():
    return {k: background(v) for k, v in field_equations().items()}


@lru_cache
def background_rules():
    return Rules(background_equations())


@lru_cache
def linear_equations():
    return {k: linearize(field_equations()[k]) for k in ('Gtt', 'Gtr', 'Grr', 'KG', 'FLUID_r', 'BAR')}


def report(title, results):
    """ok = True/False: comprobación que puede fallar. ok = None: resultado descriptivo (se informa, no se cuenta)."""
    print(f'\n== {title}')
    for name, ok, note in results:
        tag = 'RES ' if ok is None else ('OK ' if ok else 'FAIL')
        print(f'  [{tag}] {name}' + (f'  -- {note}' if note else ''))
    return results


def ok_field(ok):
    """Valor de 'ok' en los JSON: None (resultado descriptivo) se conserva."""
    return None if ok is None else bool(ok)
