# archive_unvalidated (not ported)

Nothing is copied here on purpose. The following stay in the legacy `archive/axion_ns/` and are **not** part of this layout:

| legacy files | why |
|---|---|
| `step11_gmodes.py`, `gmodes.py`, `gmodes_ref.py` | line 1 (g-modes / interface mode): hypothesis H5 RETIRED (constants `CHI_WALL`, `DELTA_COMP` and a fitted factor set by hand) |
| `step12_rotation.py`, `rotation.py`, `rotation_ref.py` | line 2 (slow rotation): H6 retired in its interpretation (rewriting of Hartle's equation, not a new source) |
| `step13_magnetar.py`, `magnetar.py`, `magnetar_ref.py` | line 3 (axion + magnetar): H7 REFUTED (mu_5,eff assigned by hand) |

See `LAB_NOTEBOOK.md` section 1 (H5-H7) and section 3 (the fabricated-report case). They are also excluded from the legacy
`run_all.py`. If one of these lines is ever reopened it should enter this layout as a new module with its own
reference check made WITHOUT the axion first (rule of `archive/axion_ns/AGENTS.md` section 4).
