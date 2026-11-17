# Each target runs ONE module with the upstream artifacts already on disk; `all` runs them in dependency order.
PY ?= python3
export PYTHONHASHSEED = 0          # SymPy simplify is hash-seed dependent (see axion_ns/AGENTS.md): fixed for reproducible artifacts

.PHONY: all R D0 D1 D2 S1 S1-witness N1 N1-witness T T-fast report test test-slow clean help

help:
	@echo "targets: all R D0 D1 D2 S1 S1-witness N1 N1-witness T T-fast report test test-slow clean"

all: R D0 D1 D2 S1 N1 T report

R:           ; $(PY) -m R_paper_refs.run
D0:          ; $(PY) -m D0_derive_action.run
D1:          ; $(PY) -m D1_derive_radial.run
D2:          ; $(PY) -m D2_derive_nonradial.run
S1:          ; $(PY) -m S1_star_background.run
S1-witness:  ; $(PY) -m S1_star_background.run --witness
N1:          ; $(PY) -m N1_modes_cowling.run
N1-witness:  ; $(PY) -m N1_modes_cowling.run --nu-a 6
T:           ; $(PY) -m T_toy_model.run
T-fast:      ; $(PY) -m T_toy_model.run --part bridge
report:      ; $(PY) -m report.build

test:        ; $(PY) -m pytest -q tests
test-slow:   ; $(PY) -m tests.verify_vs_legacy_n1

clean:       ; rm -f */outputs/manifest.json report/index.html
