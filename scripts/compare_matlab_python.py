"""Compare export_parity_results.m output with the Python implementation."""
import sys
from epidemic_modeling.validation import compare_matlab_python
compare_matlab_python(sys.argv[1] if len(sys.argv)>1 else 'reports/matlab_results.mat')
