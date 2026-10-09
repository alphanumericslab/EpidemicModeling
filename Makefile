.PHONY: test notebooks parity package

test:
	python -m pytest

notebooks:
	python scripts/execute_notebooks.py

parity:
	python scripts/compare_matlab_python.py reports/matlab_results.mat

package:
	python scripts/package_release.py
