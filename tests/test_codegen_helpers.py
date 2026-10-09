"""Standalone Coder helpers agree with the six-state analytic model."""

from types import SimpleNamespace
import numpy as np
from epidemic_modeling import codegen, default_si_params
from epidemic_modeling.kalman import _si_handles


def test_standalone_helper_protocol():
    """Check the standalone Coder callback shapes and conventions."""

    p = default_si_params([3, 3])
    p.update(a=np.array([0.06, 0.04]), b=0.02)
    handles = _si_handles(controlled=True, nonnegative_tie=True)
    np.testing.assert_array_equal(codegen.obs_hard_margins(np.array([-0.1]), p), [-0.1])
    s = np.array([0.8, 0.1, 0.3, 0.05, 0.02, 0.03])

    u = np.array([1.0, 2.0])
    w = np.zeros(6)

    for name, args in [
        ("state_hard_margins", (s,)),
        ("nlin_state_update", (u, s, w)),
        ("nlin_obs_update", (u, s, 0)),
        ("state_jacobians", (u, s, w)),
        ("obs_jacobian", (u, s, 0)),
        ("state_hessian_terms", (u, s, np.eye(6), w, np.eye(6))),
        ("obs_hessian_terms", (u, s, np.eye(6), 0, 1)),
    ]:
        result = getattr(codegen, name)(*args, p)
        reference = getattr(handles, name)(*args, SimpleNamespace(**p), 0)

        if isinstance(result, tuple):
            for actual, expected in zip(result, reference):
                np.testing.assert_allclose(actual, expected)
        else:
            np.testing.assert_allclose(result, reference)
