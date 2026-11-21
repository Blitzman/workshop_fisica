# Each target runs ONE module with the upstream artifacts already on disk; `all` runs them in dependency order.
PY ?= python3
# SymPy simplify is hash-seed dependent (see ../archive/axion_ns/AGENTS.md): fixed for reproducible artifacts. The comment must NOT share the line: make keeps the spaces before a '#' in the value.
export PYTHONHASHSEED = 0

.PHONY: all R D0 D1 D2 S1 S1-extra S1-witness N1 N1-witness N2 N2-quick N3 N3-quick N4 N4-quick T T-fast report test clean help

help:
	@echo "targets: all R D0 D1 D2 S1 S1-extra S1-witness N1 N1-witness N2 N2-quick N3 N3-quick N4 N4-quick T T-fast report test clean"

all: R D0 D1 D2 S1 S1-extra N1 N2 N3 N4 T report

R:           ; $(PY) -m R_paper_refs.run
D0:          ; $(PY) -m D0_derive_action.run
D1:          ; $(PY) -m D1_derive_radial.run
D2:          ; $(PY) -m D2_derive_nonradial.run
S1:          ; $(PY) -m S1_star_background.run
S1-witness:  ; $(PY) -m S1_star_background.run --witness
S1-extra:    ; $(PY) -m S1_star_background.run --extra
N1:          ; $(PY) -m N1_modes_cowling.run
N1-witness:  ; $(PY) -m N1_modes_cowling.run --nu-a 6
N2:          ; $(PY) -m N2_modes_fullgr_noaxion.run
N2-quick:    ; $(PY) -m N2_modes_fullgr_noaxion.run --quick --out /tmp/n2_quick
N3:          ; $(PY) -m N3_modes_fullgr_axion.run
N3-quick:    ; $(PY) -m N3_modes_fullgr_axion.run --quick --out /tmp/n3_quick
N4:          ; $(PY) -m N4_modes_radial.run --part all
N4-quick:    ; $(PY) -m N4_modes_radial.run --part fast --out /tmp/n4_quick
T:           ; $(PY) -m T_toy_model.run
T-fast:      ; $(PY) -m T_toy_model.run --part bridge
report:      ; $(PY) -m report.build

test:        ; $(PY) -m pytest -q tests

clean:       ; rm -f */outputs/manifest.json report/index.html
