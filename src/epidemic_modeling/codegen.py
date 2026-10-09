"""Python counterparts to the standalone MATLAB Coder SI-alpha helpers.

These use the six-state model, phi >= 0 switching, and zero Hessian terms.
Author: Reza Sameni, Emory University.
"""

from types import SimpleNamespace
import numpy as np
from .kalman import _si_handles


def _call(name, *args):
    """Dispatch a Coder helper with a dict or attribute-based parameter object."""
    args = list(args)

    if isinstance(args[-1], dict):
        args[-1] = SimpleNamespace(**args[-1])

    return getattr(_si_handles(controlled=True, nonnegative_tie=True), name)(*args, 0)


def state_hard_margins(s, params):
    """Clip s[0:2] to [0,1] and contact rate to its configured bounds.

    u is an intervention vector; s contains [s, i, alpha] followed by
    three costates. params is a dict or object with the model fields from
    default_si_params plus epsilon, intervention costs w, and sigma.
    Returns a length-6 state vector; costates remain unchanged.
    """

    return _call("state_hard_margins", s, params)


def obs_hard_margins(x, params):
    """Leave observations unchanged, as in the original standalone Coder helper.

    x is the observation vector. params is accepted for callback
    compatibility. Returns x unchanged.
    """

    return np.asarray(x)


def nlin_state_update(u, s, w_bar, params):
    """Return optimal control and the next six-state SI-alpha/costate vector.

    u is an intervention vector; s contains [s, i, alpha] followed by
    three costates. params is a dict or object with the model fields from
    default_si_params plus epsilon, intervention costs w, and sigma.
    Returns (selected_controls, next_state). NaN controls request a bound;
    zero switching values select the lower bound. SI process-noise means
    are accepted but not added to this historical transition.
    """

    return _call("nlin_state_update", u, s, w_bar, params)


def nlin_obs_update(u, s, v_bar, params):
    """Evaluate incidence plus v_bar; the original Coder model observes NEWCASES only.

    u is an intervention vector; s contains [s, i, alpha] followed by
    three costates. params is a dict or object with the model fields from
    default_si_params plus epsilon, intervention costs w, and sigma.
    v_bar is the observation-noise mean. Returns a length-1 incidence
    fraction; multiply by population to convert it to case counts.
    """

    return np.atleast_1d(np.prod(np.asarray(s)[:3]) + v_bar)


def state_jacobians(u, s, w_bar, params):
    """Return analytic state and process-noise Jacobians (6x6 each).

    u is an intervention vector; s contains [s, i, alpha] followed by
    three costates. params is a dict or object with the model fields from
    default_si_params plus epsilon, intervention costs w, and sigma.
    w_bar is accepted for callback compatibility. Returns A and B with
    shape (6, 6); B is the identity process-noise map.
    """

    return _call("state_jacobians", u, s, w_bar, params)


def obs_jacobian(u, s, v_bar, params):
    """Return incidence and noise Jacobians (1x6 and 1x1), as in the Coder model.

    u is an intervention vector; s contains [s, i, alpha] followed by
    three costates. params is a dict or object with the model fields from
    default_si_params plus epsilon, intervention costs w, and sigma.
    v_bar is accepted for callback compatibility. Returns C with shape
    (1, 6) and D with shape (1, 1).
    """

    return np.array([[s[1] * s[2], s[0] * s[2], s[0] * s[1], 0, 0, 0]]), np.ones((1, 1))


def state_hessian_terms(u, s, covariance, w_bar, q_w, params):
    """Return the four original zero-valued state/noise correction arrays.

    u is an intervention vector; s contains [s, i, alpha] followed by
    three costates. params is a dict or object with the model fields from
    default_si_params plus epsilon, intervention costs w, and sigma.
    covariance and q_w are state and process covariances, each (6, 6).
    Returns two length-6 mean corrections and two (6, 6) covariance
    corrections in (fs, Cs, fw, Cw) order; all entries are zero.
    """

    return _call("state_hessian_terms", u, s, covariance, w_bar, q_w, params)


def obs_hessian_terms(u, s, covariance, v_bar, r_v, params):
    """Return the four original zero-valued observation correction arrays.

    u is an intervention vector; s contains [s, i, alpha] followed by
    three costates. params is a dict or object with the model fields from
    default_si_params plus epsilon, intervention costs w, and sigma.
    covariance is (6, 6); r_v is the observation variance. Returns four
    (1,) or (1, 1) zero arrays in (gs, Gsp, gv, Gvp) order.
    """

    return _call("obs_hessian_terms", u, s, covariance, v_bar, r_v, params)


def new_case_ekf_estimator_with_optimal_npi(
    u,
    x,
    params,
    s_init,
    ps_init,
    s_final,
    ps_final,
    w_bar,
    v_bar,
    q_w,
    r_v,
    beta=1,
    gamma=1,
    inv_monitor_len=21,
    order=1,
):
    """Run the standalone Coder six-state estimator in its alternate output order.

    The original Coder model observes NEWCASES only and leaves predicted
    observations unclipped. Returns u_opt, s_minus, s_plus, p_minus, p_plus,
    k_gain, s_smooth, p_smooth, innovations, rho. Other argument contracts
    match the dedicated estimator in kalman.py.
    """
    from .kalman import generic_extended_kalman_filter

    p = SimpleNamespace(**(params if isinstance(params, dict) else vars(params)))
    p.obs_type = "NEWCASES"
    handles = _si_handles(controlled=True, nonnegative_tie=True)

    handles.obs_hard_margins = lambda x, p, k: x
    result = generic_extended_kalman_filter(
        u,
        x,
        handles,
        p,
        s_init,
        ps_init,
        s_final,
        ps_final,
        w_bar,
        v_bar,
        q_w,
        r_v,
        beta,
        gamma,
        inv_monitor_len,
        order,
        covariance_update="simple",
        switching_rho_epsilon=False,
    )

    return (
        result.u_opt,
        result.s_minus,
        result.s_plus,
        result.p_minus,
        result.p_plus,
        result.k_gain,
        result.s_smooth,
        result.p_smooth,
        result.innovations,
        result.rho,
    )
