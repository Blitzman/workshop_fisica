PYTHON ?= python3
.PHONY: cooling test test-all test-slow test-structure test-emissivities test-heat-conduction test-envelope test-integration test-output

cooling:
	$(PYTHON) cooling_driver.py

# Tests (decision D-028): run the area you changed, or everything on request. Nothing runs by itself.
test test-all:
	$(PYTHON) -m pytest tests -q
test-slow:
	$(PYTHON) -m pytest tests -q -m slow
test-structure:
	$(PYTHON) -m pytest tests/test_structure.py -q
test-emissivities:
	$(PYTHON) -m pytest tests/test_emissivities.py -q
test-heat-conduction:
	$(PYTHON) -m pytest tests/test_heat_conduction.py -q
test-envelope:
	$(PYTHON) -m pytest tests/test_envelope.py -q
test-integration:
	$(PYTHON) -m pytest tests/test_integration.py -q
test-output:
	$(PYTHON) -m pytest tests/test_output.py -q
