# Agent Protocol for `ns_cooling/` (COOL1D)

The `axion_ns/`-specific rules (`axion_ns/AGENTS.md` and `CONVENTIONS.md`) do **not** apply to `ns_cooling/` (decision D-026).

Before modifying code or running calculations in `ns_cooling/`, check `ns_cooling/README.txt` (especially `CURRENT RATE SELECTION`, `WORKING RULES`, and `REFERENCES AND FORMULAS`).

## 1. Working Rules (decision D-026, from `ns_cooling/README.txt`)
- **Hypothesis & outcome tracking:** Before testing a non-trivial idea, write the hypothesis (in the tracker); afterwards record whether it was confirmed or refuted, and why.
- **External validation:** A new result counts as validated only if it first reproduces a published result by other authors. A one-off comparison recorded in the tracker is enough; it does not have to become a test.
- **Group review & preliminary framing:** The physical premise of a new line of work is reviewed by the group before code is written; results are presented as preliminary.
- **Free numbers allowed with provenance:** Free numbers are allowed. A number that is not a physical constant or a cited value (a fit parameter, a default, an approximation) carries a comment that says where it comes from, or that its source is not recorded.

## 2. Testing (decision D-028, from `Makefile` and `pytest.ini`)
- Run only the test file for the code area you changed, or the full suite on request (`nothing runs by itself`):
  - `make test-structure`, `make test-emissivities`, `make test-heat-conduction`, `make test-envelope`, `make test-integration`, `make test-output`
  - `make test-all` (skips `slow` by default via `pytest.ini`) or `make test-slow` (for mesh convergence).
- When using the repository virtualenv from `ns_cooling/`, pass `PYTHON=../.venv/bin/python` (and note `pytest` is listed in `ns_cooling/requirements-dev.txt`).

## 3. Key Paths & Context
- **Code & physics documentation:** `ns_cooling/README.txt`
- **Entry point:** `ns_cooling/cooling_driver.py` (use `--no-animation` for headless runs)
- **Historical benchmark & comparison notes:** `docs/ns_cooling_lab_notebook.md` and `docs/cooling_comparison.html`
- **Observational data:** `observations/` (`marino2024_sources.json`, `ns_thermal_sample.json`, `ns_thermal_catalogue.html`)
