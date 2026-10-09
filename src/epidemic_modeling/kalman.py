"""Extended Kalman filtering, smoothing, and SI-alpha model callbacks.

Arrays follow MATLAB's variables-by-time convention. Author: Reza Sameni,
Emory University. Reference: doi:10.1109/JSTSP.2021.3129118.
"""

from collections import namedtuple
from types import SimpleNamespace
import numpy as np

filter_result = namedtuple(
    "filter_result",
    "u_opt u_opt_smooth s_minus s_plus s_smooth p_minus p_plus p_smooth k_gain innovations rho",
)


def _covariance_series(value, count):
    """Normalize a scalar, square covariance, or time-varying covariance stack."""
    a = np.asarray(value, float)

    if a.ndim == 0:
        a = a.reshape(1, 1)

    if a.ndim == 1:
        if a.size != count:
            raise ValueError("time-dependent scalar covariance must have N samples")

        return a.reshape(1, 1, count).copy(), False

    if a.ndim == 2 and a.shape[0] == a.shape[1]:
        return np.repeat(a[:, :, None], count, axis=2), True

    if a.ndim == 3 and a.shape[0] == a.shape[1] and a.shape[2] == count:
        return a.copy(), False

    raise ValueError("covariance must be scalar, square, or square-by-N")


def _right_divide(a, b):
    """Compute MATLAB-style matrix right division for nonsingular matrices."""

    return np.linalg.solve(b.T, a.T).T


def generic_extended_kalman_filter(
    u,
    x,
    handles,
    params,
    s_init,
    ps_init,
    s_final,
    ps_final,
    w_bar,
    v_bar,
    q_w,
    r_v,
    beta=1.0,
    gamma=1.0,
    inv_monitor_len=21,
    order=1,
    *,
    covariance_update="joseph",
    switching_rho_epsilon=True
):
    """Run an EKF and fixed-interval Rauch--Tung--Striebel smoother.

    Parameters
    ----------
    u, x : array_like, shape (controls, N), (observations, N)
        Inputs and observations. A column with any NaN observation skips
        correction. Model callbacks may replace NaN controls with optima.
    handles : mapping or object
        Eight snake-case callbacks: state_hard_margins, obs_hard_margins,
        nlin_state_update, nlin_obs_update, state_jacobians, obs_jacobian,
        state_hessian_terms, obs_hessian_terms. Signatures match MATLAB,
        except time index k is zero-based. See docs/api.md.
    params : object
        Passed unchanged to each callback.
    s_init, ps_init : array_like
        Initial mean (m,) and covariance (m,m).
    s_final, ps_final : array_like
        Final smoothing boundary values; NaN entries remain unconstrained.
    w_bar, v_bar : array_like
        Process and observation noise means.
    q_w, r_v : array_like
        Fixed covariance matrices, scalar variance series, or covariance
        stacks with shape (dimension, dimension, N).
    beta, gamma : float
        Observation-noise adaptation factor [0,1] and covariance stability
        factor (0,1], respectively. beta=gamma=1 disables these adjustments.
    inv_monitor_len : int
        Positive length of the rolling innovation monitor.
    order : {1, 2}
        Linearized EKF or callback-supplied second-order corrections.

    Returns
    -------
    filter_result
        Tuple of eleven arrays in the original generic MATLAB output order.
        Covariances and gains retain a final time axis; rho is squeezed.
        The final smoothed control is zero, preserving the legacy convention.

    Notes
    -----
    Fixed R may adapt; time-dependent R does not. Joseph covariance updates
    and symmetrization follow the generic original. The optional simple
    update is used only by the two older dedicated filter ports.
    """

    if not 0 <= beta <= 1 or not 0 < gamma <= 1 or order not in (1, 2):
        raise ValueError("beta must be in [0,1], gamma in (0,1], order in {1,2}")

    if inv_monitor_len < 1 or int(inv_monitor_len) != inv_monitor_len:
        raise ValueError("inv_monitor_len must be a positive integer")

    if isinstance(handles, dict):
        handles = SimpleNamespace(**handles)

    x = np.asarray(x, float)
    u = np.asarray(u, float)

    if x.ndim == 1:
        x = x.reshape(1, -1)

    if u.ndim == 1:
        u = u.reshape(1, -1)

    if x.ndim != 2 or u.ndim != 2 or not x.shape[1] or u.shape[1] != x.shape[1]:
        raise ValueError("x and u must be matrices with the same nonzero time length")

    n, count = x.shape
    s = np.asarray(s_init, float).reshape(-1).copy()
    m = len(s)
    p = np.asarray(ps_init, float).copy()

    if p.shape != (m, m) or not np.all(np.isfinite(s)) or not np.all(np.isfinite(p)):
        raise ValueError(
            "initial states/covariance must be finite with matching dimensions"
        )

    q, _ = _covariance_series(q_w, count)
    r, fixed_r = _covariance_series(r_v, count)
    sm = np.zeros((m, count))
    sp = sm.copy()
    pm = np.zeros((m, m, count))

    pp = pm.copy()
    kg = np.zeros((m, n, count))
    innovations = np.zeros((n, count))
    rho = np.zeros((n, n, count))
    u_opt = np.zeros_like(u)

    u_smooth = np.zeros_like(u)
    history = np.zeros((n, inv_monitor_len))
    covhist = np.zeros((n, n, inv_monitor_len))
    normhist = covhist.copy()

    for k in range(count):
        sm[:, k], pm[:, :, k] = s, p
        gs, gsp, gv, gvp = (
            np.zeros(n),
            np.zeros((n, n)),
            np.zeros(n),
            np.zeros((n, n)),
        )

        if order == 2:
            gs, gsp, gv, gvp = handles.obs_hessian_terms(
                u[:, k], s, p, v_bar, r[:, :, k], params, k
            )

        c, d = handles.obs_jacobian(u[:, k], s, v_bar, params, k)
        c = np.atleast_2d(c)
        d = np.atleast_2d(d)
        predicted = handles.obs_hard_margins(
            np.asarray(handles.nlin_obs_update(u[:, k], s, v_bar, params, k)).reshape(
                -1
            )
            + gs
            + gv,
            params,
            k,
        )
        valid = np.all(np.isfinite(x[:, k]))

        if valid:
            innovations[:, k] = x[:, k] - predicted
            gain = _right_divide(
                p @ c.T, c @ p @ c.T + gamma * d @ r[:, :, k] @ d.T + gsp + gvp
            )
            identity = np.eye(m) - gain @ c
            posterior_p = (
                (identity @ p @ identity.T + gain @ d @ r[:, :, k] @ d.T @ gain.T)
                / gamma
                if covariance_update == "joseph"
                else identity @ p / gamma
            )
            posterior_s = s + gain @ innovations[:, k]
        else:
            gain = np.zeros((m, n))
            posterior_s = s.copy()
            posterior_p = p.copy()

        if covariance_update == "joseph":
            posterior_p = (posterior_p + posterior_p.T) / 2

        posterior_s = handles.state_hard_margins(posterior_s, params, k)
        fs, fsp, fw, fwp = np.zeros(m), np.zeros((m, m)), np.zeros(m), np.zeros((m, m))

        if order == 2:
            fs, fsp, fw, fwp = handles.state_hessian_terms(
                u[:, k], posterior_s, posterior_p, w_bar, q[:, :, k], params, k
            )

        u_opt[:, k], s = handles.nlin_state_update(
            u[:, k].copy(), posterior_s, w_bar, params, k
        )
        s = s + fs + fw
        a, b = handles.state_jacobians(u[:, k], posterior_s, w_bar, params, k)
        p = a @ posterior_p @ a.T + b @ q[:, :, k] @ b.T + fsp + fwp

        if covariance_update == "joseph":
            p = (p + p.T) / 2

        s = handles.state_hard_margins(s, params, k)
        sp[:, k], pp[:, :, k], kg[:, :, k] = posterior_s, posterior_p, gain
        stats = min(k + 1, inv_monitor_len)

        # Track recent residuals for diagnostics and optional noise adaptation.
        history = np.concatenate((innovations[:, k, None], history[:, :-1]), axis=1)
        centered = innovations[:, k] - history.sum(axis=1) / stats

        cc = np.outer(centered, centered)
        covhist = np.concatenate((cc[:, :, None], covhist[:, :, :-1]), axis=2)
        denom = r[:, :, k] + (np.finfo(float).eps if switching_rho_epsilon else 0)
        normhist = np.concatenate(
            (_right_divide(cc, denom)[:, :, None], normhist[:, :, :-1]), axis=2
        )
        rho[:, :, k] = normhist.sum(axis=2) / stats

        if beta != 1 and valid and fixed_r and k + 1 < count:
            r[:, :, k + 1] = (
                beta * r[:, :, k] + (1 - beta) * covhist.sum(axis=2) / stats
            )

    # Apply only finite boundary entries, then smooth backwards in time.
    ss = sp.copy()
    ps = pp.copy()
    end_s = np.asarray(s_final, float).reshape(m)
    end_p = np.asarray(ps_final, float).reshape(m, m)
    ss[np.isfinite(end_s), -1] = end_s[np.isfinite(end_s)]

    ps[:, :, -1][np.isfinite(end_p)] = end_p[np.isfinite(end_p)]

    for k in range(count - 2, -1, -1):
        a, _ = handles.state_jacobians(u[:, k], sp[:, k], w_bar, params, k)
        j = (
            pp[:, :, k] @ a.T @ np.linalg.pinv(pm[:, :, k + 1])
            if np.all(np.isfinite(pm[:, :, k + 1]))
            else np.zeros((m, m))
        )
        ss[:, k] = handles.state_hard_margins(
            sp[:, k] + j @ (ss[:, k + 1] - sm[:, k + 1]), params, k
        )
        ps[:, :, k] = pp[:, :, k] - j @ (pm[:, :, k + 1] - ps[:, :, k + 1]) @ j.T

        if covariance_update == "joseph":
            ps[:, :, k] = (ps[:, :, k] + ps[:, :, k].T) / 2

        u_smooth[:, k], _ = handles.nlin_state_update(
            u[:, k].copy(), ss[:, k], w_bar, params, k
        )

    return filter_result(
        u_opt, u_smooth, sm, sp, ss, pm, pp, ps, kg, innovations, np.squeeze(rho)
    )


def _si_handles(controlled=False, backward=False, nonnegative_tie=False):
    """Construct analytic callbacks for the original three/six-state SI model."""
    sign = -1 if backward else 1
    m = 6 if controlled else 3

    def state_hard_margins(s, p, k):
        """Clip fractions and contact rate, leaving costates unrestricted."""
        s = s.copy()
        s[0] = np.clip(
            s[0], getattr(p, "s_min", 0) if not controlled and not backward else 0, 1
        )
        s[1] = np.clip(
            s[1], getattr(p, "i_min", 0) if not controlled and not backward else 0, 1
        )
        s[2] = np.clip(s[2], p.alpha_min, p.alpha_max)

        return s

    def obs_hard_margins(x, p, k):
        """Constrain predicted cases to be nonnegative."""

        return np.maximum(0, x)

    def nlin_state_update(u, s, w, p, k):
        """Advance states/costates and replace missing controls by optima."""
        u = u.copy()

        if controlled:

            # The switching function selects the lower or upper control bound.
            phi = p.epsilon * np.asarray(p.w) - p.gamma * s[5] * np.asarray(p.a)
            lower = phi >= 0 if nonnegative_tie else phi > 0
            u = np.where(np.isnan(u), np.where(lower, p.u_min, p.u_max), u)

        s0, i, alpha = s[:3]
        inc = alpha * s0 * i
        derivative = [
            -inc,
            inc - p.beta * i,
            -p.gamma * alpha
            + p.gamma * p.b
            + p.gamma * np.asarray(p.a) @ (np.asarray(p.u_max) - u),
        ]

        if controlled:
            rr = s[3] - s[4] - (1 - p.epsilon)
            derivative += [
                rr * alpha * i,
                rr * alpha * s0 + p.beta * s[4],
                rr * s0 * i + p.gamma * s[5],
            ]

        return u, state_hard_margins(s + sign * p.dt * np.asarray(derivative), p, k)

    def nlin_obs_update(u, s, v, p, k):
        """Observe incidence or cumulative cases according to obs_type."""

        if p.obs_type == "NEWCASES":
            return np.atleast_1d(np.prod(s[:3]) + v)

        if p.obs_type == "TOTALCASES":
            return np.atleast_1d(1 - s[0] + v)

        raise ValueError("obs_type must be NEWCASES or TOTALCASES")

    def state_jacobians(u, s, w, p, k):
        """Evaluate the historical analytic state and process-noise Jacobians."""
        s0, i, alpha = s[:3]
        dt = sign * p.dt
        a = np.zeros((m, m))
        a[0, :3] = [-alpha * i, -alpha * s0, -s0 * i]

        a[1, :3] = [i * alpha, s0 * alpha - p.beta, s0 * i]
        a[2, 2] = -p.gamma

        if controlled:

            # The switching function selects the lower or upper control bound.
            phi = p.epsilon * np.asarray(p.w) - p.gamma * s[5] * np.asarray(p.a)
            active = np.isnan(u) & (phi > -1 / p.sigma) & (phi < 1 / p.sigma)
            a[2, 5] = (
                -p.gamma
                * (p.sigma / 2)
                * np.sum(
                    np.asarray(p.a)[active]
                    * (np.asarray(p.u_max) - np.asarray(p.u_min))[active]
                )
            )
            rr = s[3] - s[4] - (1 - p.epsilon)
            a[3, 1:5] = [alpha * rr, i * rr, i * alpha, -i * alpha]

            a[4, [0, 2, 3, 4]] = [alpha * rr, s0 * rr, s0 * alpha, -s0 * alpha + p.beta]
            a[5, [0, 1, 3, 4, 5]] = [i * rr, s0 * rr, s0 * i, -s0 * i, p.gamma]

        return np.eye(m) + dt * a, np.eye(m)

    def obs_jacobian(u, s, v, p, k):
        """Differentiate incidence or cumulative-case observations."""
        c = np.zeros((1, m))

        if p.obs_type == "NEWCASES":
            c[0, :3] = [s[1] * s[2], s[0] * s[2], s[0] * s[1]]
        elif p.obs_type == "TOTALCASES":
            c[0, 0] = -1
        else:
            raise ValueError("obs_type must be NEWCASES or TOTALCASES")

        return c, np.ones((1, 1))

    def state_hessian_terms(u, s, pk, w, q, p, k):
        """Return zero corrections, matching the SI routines' original order-2 implementation."""

        return np.zeros(m), np.zeros((m, m)), np.zeros(m), np.zeros((m, m))

    def obs_hessian_terms(u, s, pk, v, r, p, k):
        """Return the original zero observation Hessian corrections."""

        return np.zeros(1), np.zeros((1, 1)), np.zeros(1), np.zeros((1, 1))

    return SimpleNamespace(**{n: v for n, v in locals().items() if callable(v)})


def _si_filter(
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
    beta,
    gamma,
    inv_monitor_len,
    order,
    controlled=False,
    backward=False,
    tie=False,
    simple=False,
):
    """Share the SI model setup while retaining forward/backward output order."""
    p = SimpleNamespace(**params) if isinstance(params, dict) else params
    args = [
        np.asarray(u, float),
        np.asarray(x, float),
        _si_handles(controlled, backward, tie),
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
    ]

    if backward:
        args[0] = args[0][..., ::-1]
        args[1] = args[1][..., ::-1]
        args[4:8] = [s_final, ps_final, s_init, ps_init]

        for i in (10, 11):
            v = np.asarray(args[i])

            if v.ndim in (1, 3):
                args[i] = v[..., ::-1]

    result = generic_extended_kalman_filter(
        *args,
        covariance_update="simple" if simple else "joseph",
        switching_rho_epsilon=not simple
    )

    if backward:
        result = filter_result(*(a[..., ::-1] if a.ndim else a for a in result))

    return result


def si_alpha_model_ekf(
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
    """Estimate three SI-alpha states; inputs/outputs follow the generic EKF contract.

    The state is [susceptible fraction, infected fraction, contact rate].
    u is interventions-by-time; x is observations-by-time and contains
    population fractions, not raw case counts. params is a dict or object
    with the fields from default_si_params. Select NEWCASES for incidence
    or TOTALCASES for cumulative cases using params.obs_type.

    Initial/final means have length 3 and covariances have shape (3, 3).
    NaN observations skip correction; NaN final entries impose no boundary.
    Returns filter_result with filtered and smoothed states, covariances,
    controls, gains, and innovations. See generic_extended_kalman_filter
    for noise inputs and adaptation settings.
    """

    return _si_filter(
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
        beta,
        gamma,
        inv_monitor_len,
        order,
    )


def si_alpha_model_ekf_opt_controlled(
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
    """Estimate six SI-alpha/costate states and fill NaN controls using phi > 0.

    The six states are susceptible fraction, infected fraction, contact
    rate, and three costates. Means have length 6; covariances are (6, 6).
    u is interventions-by-time; NaN controls request a bound selected by
    the switching function. Finite controls are used as supplied.

    params adds epsilon, intervention costs w, and switching slope sigma
    to default_si_params. Returns the same eleven-array filter_result as
    generic_extended_kalman_filter. Observations use population fractions.
    NaN final mean/covariance entries leave that boundary unconstrained.
    """

    return _si_filter(
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
        beta,
        gamma,
        inv_monitor_len,
        order,
        controlled=True,
    )


def si_alpha_model_backward_ekf(
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
    """Run the three-state reverse-time SI filter and return chronological arrays.

    Inputs use the same units, shapes, and parameter fields as
    si_alpha_model_ekf. Supply samples in chronological order; this function
    reverses them internally. The initial mean/covariance starts the
    reverse-time pass at the latest sample. Returned arrays are chronological.
    """

    return _si_filter(
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
        beta,
        gamma,
        inv_monitor_len,
        order,
        backward=True,
    )


def si_alpha_model_backward_ekf_opt_controlled(
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
    """Run six-state reverse-time filtering with bounded optimal NPI controls.

    Inputs follow si_alpha_model_ekf_opt_controlled, with six states and
    NaN controls requesting optimal bounds. Supply chronological data.
    The initial mean/covariance starts at the latest sample; the final
    boundary applies at the earliest sample. Outputs are chronological.
    """

    return _si_filter(
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
        beta,
        gamma,
        inv_monitor_len,
        order,
        controlled=True,
        backward=True,
    )


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
    """Port the older six-state estimator (simple covariance; switching phi >= 0).

    Returns ten arrays: u_opt, s_minus, s_plus, s_smooth, p_minus, p_plus,
    p_smooth, k_gain, innovations, rho, as in Tools' original function.
    """
    r = _si_filter(
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
        beta,
        gamma,
        inv_monitor_len,
        order,
        controlled=True,
        tie=True,
        simple=True,
    )

    return (r.u_opt,) + tuple(r[2:])


def _hessian_moments(cov, hessians):
    """Compute Gaussian second-order mean and covariance correction terms."""
    means = np.array([np.trace(cov @ h) / 2 for h in hessians])
    variances = np.array(
        [[np.trace(cov @ hi @ cov @ hj) / 2 for hj in hessians] for hi in hessians]
    )

    return means, variances


def rt_exp_fit_ekf(
    x,
    s_init,
    params,
    w_bar,
    v_bar,
    ps_init,
    q_w,
    r_v,
    beta=1,
    gamma=1,
    inv_monitor_len=21,
    order=1,
):
    """Estimate cases and bounded growth with the original two-state EKF/EKS.

    params = [time_scale, growth_memory, growth_saturation]; saturation > 0.
    Returns nine arrays in the MATLAB order: s_minus, s_plus, p_minus,
    p_plus, k_gain, s_smooth, p_smooth, innovations, rho. Missing observations
    skip correction. Both first- and second-order filters are supported.

    x is a length-N case series. The two states are case amplitude and
    growth; s_init has length 2 and ps_init is (2, 2). q_w and r_v are
    process and observation covariances. Noise means are w_bar and v_bar.
    All state arrays have shape (2, N); covariance arrays are (2, 2, N).
    """
    time_scale, memory, sigma = np.asarray(params, float)

    if sigma <= 0:
        raise ValueError("growth saturation must be positive")

    w = np.asarray(w_bar, float).reshape(2)

    def state_hard_margins(s, p, k):
        """Leave exponential state estimates unconstrained, matching MATLAB."""

        return s

    def obs_hard_margins(x, p, k):
        """Leave the linear observation unchanged."""

        return x

    def nlin_state_update(u, s, w, p, k):
        """Propagate amplitude and tanh-bounded growth."""

        return u, np.array(
            [
                s[0] * np.exp(time_scale * s[1]) + w[0],
                sigma * np.tanh((memory * s[1] + w[1]) / sigma),
            ]
        )

    def nlin_obs_update(u, s, v, p, k):
        """Observe the case amplitude."""

        return np.atleast_1d(s[0] + v)

    def state_jacobians(u, s, w, p, k):
        """Differentiate the exponential and saturated state equations."""
        e = np.exp(time_scale * s[1])
        t = np.tanh((memory * s[1] + w[1]) / sigma)

        return np.array(
            [[e, time_scale * s[0] * e], [0, memory * (1 - t * t)]]
        ), np.diag([1, 1 - t * t])

    def obs_jacobian(u, s, v, p, k):
        """Return the linear case observation Jacobian."""

        return np.array([[1, 0]]), np.ones((1, 1))

    def state_hessian_terms(u, s, pk, w, q, p, k):
        """Compute analytic second-order Gaussian corrections."""
        e = np.exp(time_scale * s[1])
        t = np.tanh((memory * s[1] + w[1]) / sigma)
        hs = [
            np.array([[0, time_scale * e], [time_scale * e, time_scale**2 * s[0] * e]]),
            np.diag([0, -2 * memory**2 / sigma * t * (1 - t * t)]),
        ]
        hw = [np.zeros((2, 2)), np.diag([0, -2 / sigma * t * (1 - t * t)])]

        fs, cs = _hessian_moments(pk, hs)
        fw, cw = _hessian_moments(q, hw)

        return fs, cs, fw, cw

    def obs_hessian_terms(u, s, pk, v, r, p, k):
        """Return zero corrections for the linear observation."""

        return np.zeros(1), np.zeros((1, 1)), np.zeros(1), np.zeros((1, 1))

    handles = SimpleNamespace(**{n: v for n, v in locals().items() if callable(v)})

    count = np.asarray(x).shape[-1]
    result = generic_extended_kalman_filter(
        np.zeros((1, count)),
        x,
        handles,
        params,
        s_init,
        ps_init,
        np.full(2, np.nan),
        np.full((2, 2), np.nan),
        w,
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
