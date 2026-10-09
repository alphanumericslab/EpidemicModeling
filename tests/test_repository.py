"""Check conversion coverage, imports, naming and notebook structure."""
import ast
import importlib
import json
from pathlib import Path
import re
import nbformat
from epidemic_modeling.validation import parity_results

ROOT=Path(__file__).resolve().parents[1]


def test_every_original_function_has_python_port():
    names=json.loads((ROOT/'docs/name_map.json').read_text())
    original_public=['ForecastQualityAssessment','GenericExtendedKalmanFilter','MyTanhLayer','NPICost','NewCaseEKFEstimatorWithOptimalNPI',
        'PrescribeNPI','ReadCOVID19Data','Rt_ExpFitEKF','Rt_ExpFitGenRatios','Rt_ExpFitLogLinReg','Rt_ExpFitNonlinLS',
        'SEIRP','SEIRPSaturatedResource','SIAlphaModelBackwardEKF','SIAlphaModelBackwardEKFOptControlled','SIAlphaModelEKF',
        'SIAlphaModelEKFOptControlled','SI_Controlled','SIalpha_Controlled','TrainNPIPrescriptor','TrainPredictPrescribeNPI','expLayer']
    package=importlib.import_module('epidemic_modeling')
    for old in original_public:
        assert callable(getattr(package,names[old]))
        assert (ROOT/'matlab/functions'/f'{names[old]}.m').is_file()


def test_python_numerical_functions_are_documented_snake_case():
    for path in (ROOT/'src').rglob('*.py'):
        for node in ast.walk(ast.parse(path.read_text())):
            if isinstance(node,(ast.FunctionDef,ast.AsyncFunctionDef)):
                assert re.fullmatch(r'[a-z_][a-z0-9_]*',node.name), (path,node.name)
                assert ast.get_docstring(node), (path,node.name)


def test_matlab_function_and_filename_names():
    for path in (ROOT/'matlab').rglob('*.m'):
        for name in re.findall(r'^\s*function\s+(?:\[[^\]]*\]\s*=\s*|\w+\s*=\s*)?(\w+)\s*\(',path.read_text(),re.M):
            assert re.fullmatch(r'[a-z][a-z0-9_]*',name),(path,name)


def test_notebooks_import_no_model_definitions_and_credit_author():
    notebooks=sorted((ROOT/'notebooks').glob('*.ipynb')); assert len(notebooks)==7
    for path in notebooks:
        nb=nbformat.read(path,as_version=4); nbformat.validate(nb)
        text='\n'.join(c.source for c in nb.cells)
        assert 'Reza Sameni' in text and 'Emory University' in text
        assert '2003.11371' in text and '3129118' in text
        for cell in nb.cells:
            if cell.cell_type=='code':
                assert not any(isinstance(n,(ast.FunctionDef,ast.ClassDef)) for n in ast.walk(ast.parse(cell.source)))


def test_parity_case_suite_is_finite():
    import numpy as np
    results=parity_results(); assert len(results)>100
    assert all(np.all(np.isfinite(x)) for x in results.values())
