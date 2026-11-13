"""File provider of the background: the input file and the EOS star it defines. No plotting."""

from pathlib import Path

import numpy as np

from stellar_structure import STRUCTURE_KEYS, StructureOptions, build_stellar_models
from superfluid import GAP_MODELS
from thermal_fluxes import MODELS as ENVELOPE_MODELS

OPTIONS = ("envelope_model",) + STRUCTURE_KEYS
_TEXT_OPTIONS = ("eos", "mesh")
_LIST_OPTIONS = ("mass", "central_pressure")  # one value per star


def _parse_option(key, value):
    if key == "envelope_model":
        value = value.lower()
        if value not in ENVELOPE_MODELS:
            raise ValueError(f"envelope_model must be one of {', '.join(ENVELOPE_MODELS)}")
        return value
    if key in _TEXT_OPTIONS:
        if not value:
            raise ValueError(f"Option {key} needs a value")
        return value.upper() if key == "mesh" else value
    words = value.replace(",", " ").split() if key in _LIST_OPTIONS else [value]
    try:
        numbers = tuple(float(w) for w in words)
    except ValueError:
        raise ValueError(f"Option {key} needs a number, got '{value}'") from None
    if key in _LIST_OPTIONS:
        if not numbers:
            raise ValueError(f"Option {key} needs a number or a list of numbers")
        return numbers
    return numbers[0]


def read_parameters(path):
    """Read the input file, including the optional ``envelope_model`` selector.

    Layout (0-based line numbers): 1 is not used (``none``), 3 initial T_8
    [1e8 K], 8/10/12 superfluidity switches and gap models for crust neutrons,
    core neutrons and core protons. Line 5 is not used.
    Later lines may hold ``key = value`` options: ``envelope_model`` and ``eos`` with its
    companions, which define the star (see ``stellar_structure``).
    """
    lines = Path(path).read_text().splitlines()

    def pair(i):
        flag, model = map(int, lines[i].split())
        if flag not in (0, 1) or model not in GAP_MODELS:
            raise ValueError(f"Invalid configuration on line {i + 1}")
        return {"enabled": bool(flag), "model": model}

    cfg = dict(
        structure=lines[1].strip(),
        Tinit8=float(lines[3].replace("D", "e").replace("d", "e")),
        crust_n=pair(8),
        core_n=pair(10),
        core_p=pair(12),
    )
    if not np.isfinite(cfg["Tinit8"]) or cfg["Tinit8"] <= 0:
        raise ValueError("Initial temperature must be positive and finite")
    # Only explicit key=value options are parsed; free text lines are skipped.
    cfg["envelope_model"] = "gudmundsson"
    seen = set()
    for line in lines[13:]:
        text = line.split("#", 1)[0].strip()
        if "=" not in text:
            continue
        key, value = (part.strip() for part in text.split("=", 1))
        key = key.lower()
        if key not in OPTIONS:
            raise ValueError(f"Unknown input option: {key}")
        if key in seen:
            raise ValueError(f"Duplicate {key} option")
        seen.add(key)
        cfg[key] = _parse_option(key, value)
    if any(k in cfg for k in STRUCTURE_KEYS):
        StructureOptions.list_from_config(
            cfg
        )  # validates the star options early (no files are read)
    return cfg


def load_background(input_path):
    """The ``Background`` of an input file that defines one star (see ``load_backgrounds``)."""
    cfg = read_parameters(Path(input_path).resolve())
    if "eos" in cfg:  # refuse before building anything
        n_stars = len(StructureOptions.list_from_config(cfg))
        if n_stars != 1:
            raise ValueError(f"The input file defines {n_stars} stars; use load_backgrounds")
    return load_backgrounds(input_path)[0][1]


def load_backgrounds(input_path):
    """Read the input file and return a list of ``(label, Background)``, one per star.

    The stars are built from the equation of state named by the ``eos`` option
    (``stellar_structure``; ``mass`` or ``central_pressure`` may list several stars).
    The critical temperatures come from the pairing gap models selected in the input
    file; they are computed at the faces and averaged to the cells when the background
    is built.
    """
    input_path = Path(input_path).resolve()
    cfg = read_parameters(input_path)
    if "eos" not in cfg:
        raise ValueError(
            "The input file has no eos option: the star is built from an EOS package "
            "(eos = <folder in eos/>, plus mass or central_pressure)"
        )
    models = build_stellar_models(cfg, input_path.parent / "eos")
    return [(m.label, m.background) for m in models]
