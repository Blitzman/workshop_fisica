"""Datos de las figuras en JSON, para que la página las dibuje en el navegador (report/build.py)."""
import json
from pathlib import Path

import numpy as np

RES = Path(__file__).resolve().parent.parent / 'results'


def r(x, d=4):
    """Redondea números, arrays y listas anidadas a d cifras significativas."""
    if isinstance(x, dict):
        return {k: r(v, d) for k, v in x.items()}
    if isinstance(x, (list, tuple, np.ndarray)):
        return [r(v, d) for v in x]
    if isinstance(x, (int, np.integer)) or x is None or isinstance(x, str):
        return x if not isinstance(x, np.integer) else int(x)
    x = float(x)
    return 0.0 if x == 0 else float(f'{x:.{d}g}')


def save(name, data):
    RES.mkdir(exist_ok=True)
    (RES / f'fig_{name}.json').write_text(json.dumps(data, ensure_ascii=False, separators=(',', ':')))
