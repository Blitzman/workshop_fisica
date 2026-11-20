"""N2: literature mode frequencies from R's literature_modes.ref.json, in the shapes the validation uses."""
import json
from pathlib import Path


def load_literature(path):
    """Returns (LMI97 {l: (Re w, Im w)}, YK97 {M/R: {l: {mode: (exact, Cowling)}}}, R in km of the LMI star); units sqrt(M/R^3)."""
    p = Path(path)
    if not p.exists():
        raise FileNotFoundError(f'{p} not found. Produce it first: python -m R_paper_refs.run')
    d = json.loads(p.read_text(encoding='utf-8'))
    lmi = {int(k): (v['re_w'], v['im_w']) for k, v in d['LMI97']['l'].items()}
    yk = {float(mr): {int(ll): {m: (v['exact'], v['cowling']) for m, v in modes.items()} for ll, modes in t.items()}
          for mr, t in d['YK97']['M_over_R'].items()}
    return lmi, yk, d['LMI97']['star']['R_km']
