"""Forward-Euler compartment models, matching the original MATLAB equations.

Author: Reza Sameni, Emory University.
Original MATLAB solver notices: released under the GNU General Public License.
Reference: Sameni (2020), https://arxiv.org/abs/2003.11371.
"""
import numpy as np


def _rate_series(value, count, name):
    """Expand a scalar or validate a finite time-dependent rate sequence."""
    a = np.asarray(value, dtype=float).reshape(-1)
    if a.size == 1:
        a = np.full(count, a[0])
    if a.size < count or not np.all(np.isfinite(a)) or np.any(a < 0):
        raise ValueError(f"{name} must be nonnegative, finite, and contain at least {count} samples")
    return a


def _sample_count(duration, dt):
    """Use MATLAB's positive-number round convention for the time grid."""
    if not np.isfinite(dt) or dt <= 0 or not np.isfinite(duration) or duration <= 0:
        raise ValueError("duration and dt must be finite and positive")
    k = int(np.floor(duration / dt + 0.5))
    if k < 1:
        raise ValueError("duration/dt must round to at least one sample")
    return k


def _seirp_step(state, alpha_e, alpha_i, kappa, rho, beta, mu, gamma):
    """Evaluate the five-compartment vector field at one time point."""
    s, e, i, r, p = state
    incidence = alpha_e * s * e + alpha_i * s * i
    return np.array([-incidence + gamma*r, incidence-(kappa+rho)*e,
                     kappa*e-(beta+mu)*i, beta*i+rho*e-gamma*r, mu*i])


def seirp(alpha_e, alpha_i, kappa, rho, beta, mu, gamma,
          s0, e0, i0, r0, p0, duration, dt):
    """Solve the SEIRP model using forward Euler.

    Parameters
    ----------
    alpha_e, alpha_i, kappa, rho, beta, mu, gamma : float or array_like
        Nonnegative rates per time unit, scalar or at least K-1 samples.
    s0, e0, i0, r0, p0 : float
        Initial population fractions (susceptible, exposed, infected,
        recovered, and passed/deceased).
    duration, dt : float
        Requested duration and step size in the same time unit.

    Returns
    -------
    s, e, i, r, p : ndarray, shape (K,)
        K = round(duration/dt) samples, including the initial conditions.
        The final time is (K-1)*dt, preserving the original MATLAB convention.

    Notes
    -----
    No clipping is applied: reduce dt if Euler produces negative fractions.
    The sum of the compartments is conserved up to floating-point precision.
    """
    k = _sample_count(duration, dt)
    rates = [_rate_series(v, k-1, n) for v, n in zip(
        (alpha_e, alpha_i, kappa, rho, beta, mu, gamma),
        ('alpha_e', 'alpha_i', 'kappa', 'rho', 'beta', 'mu', 'gamma'))]
    states = np.empty((5, k)); states[:, 0] = [s0, e0, i0, r0, p0]
    for t in range(k-1):
        states[:, t+1] = states[:, t] + dt*_seirp_step(states[:, t], *(v[t] for v in rates))
    return tuple(states)


def seirp_saturated_resource(alpha_e, alpha_i, kappa, rho, gamma,
                             s0, e0, i0, r0, p0, duration, dt,
                             beta_0, beta_s, mu_0, mu_s, sigma, i_0):
    """Solve SEIRP with a smooth transition to saturated healthcare rates.

    ``h = (tanh((i-i_0)/sigma)+1)/2`` interpolates recovery from beta_0 to
    beta_s and mortality from mu_0 to mu_s. sigma must be positive; i_0 is
    the prevalence threshold. Other inputs and outputs follow :func:`seirp`.
    """
    if sigma <= 0:
        raise ValueError("sigma must be positive")
    k = _sample_count(duration, dt)
    rates = [_rate_series(v, k-1, n) for v, n in zip(
        (alpha_e, alpha_i, kappa, rho, gamma), ('alpha_e','alpha_i','kappa','rho','gamma'))]
    states = np.empty((5, k)); states[:, 0] = [s0, e0, i0, r0, p0]
    for t in range(k-1):
        h = (np.tanh((states[2, t]-i_0)/sigma)+1)/2
        ae, ai, ka, ro, ga = (v[t] for v in rates)
        states[:, t+1] = states[:, t] + dt*_seirp_step(
            states[:, t], ae, ai, ka, ro, beta_0+(beta_s-beta_0)*h,
            mu_0+(mu_s-mu_0)*h, ga)
    return tuple(states)


def si_controlled(alpha, beta, s0, i0, count, dt):
    """Integrate bounded SI fractions with a supplied contact-rate schedule.

    alpha is scalar or length >= count-1; beta is the removal rate. Returns
    two length-count arrays including the initial state. Each Euler update
    is clipped to [0, 1], as in the original MATLAB solver.
    """
    if count < 1 or int(count) != count or dt <= 0:
        raise ValueError("count must be a positive integer and dt positive")
    alpha = _rate_series(alpha, count-1, 'alpha')
    states = np.empty((2, count)); states[:, 0] = [s0, i0]
    for t in range(count-1):
        s, i = states[:, t]; incidence = alpha[t]*s*i
        states[:, t+1] = np.clip([s-dt*incidence, i+dt*(incidence-beta*i)], 0, 1)
    return tuple(states)


def si_alpha_controlled(u, s0, i0, alpha0, u_max, alpha_min, alpha_max,
                        gamma, a, b, beta, s_noise_std, i_noise_std,
                        alpha_noise_std, count, dt, *, rng=None, noise=None):
    """Integrate NPI-driven SI-alpha dynamics, returning post-update samples.

    u has shape (interventions, count); a and u_max have one element per
    intervention. gamma is the contact-rate response speed. The equilibrium
    contact rate is b + a @ (u_max-u). States and contact rate are clipped.
    The returned three length-count arrays exclude the initial condition.

    Pass ``noise`` of shape (3, count) containing standard-normal draws for
    exact MATLAB/Python comparison. Otherwise ``rng`` is a NumPy Generator.
    Zero noise standard deviations yield deterministic trajectories.
    """
    u = np.asarray(u, float)
    if u.ndim == 1: u = u.reshape(1, -1)
    a = np.asarray(a, float).reshape(-1); u_max = np.asarray(u_max, float).reshape(-1)
    if u.shape != (a.size, count) or a.size != u_max.size or dt <= 0:
        raise ValueError("u must have shape (len(a), count); u_max must match a; dt must be positive")
    if noise is None: noise = (rng or np.random.default_rng()).standard_normal((3, count))
    noise = np.asarray(noise, float)
    if noise.shape != (3, count): raise ValueError("noise must have shape (3, count)")
    states = np.empty((3, count+1)); states[:, 0] = [s0, i0, alpha0]
    for t in range(count):
        s, i, alpha = states[:, t]; inc = alpha*s*i
        states[:, t+1] = [np.clip(s-dt*(inc+noise[0,t]*s_noise_std),0,1),
            np.clip(i+dt*(inc-beta*i+noise[1,t]*i_noise_std),0,1),
            np.clip(alpha+dt*(-gamma*alpha+gamma*b+gamma*a@(u_max-u[:,t])+
                             noise[2,t]*alpha_noise_std),alpha_min,alpha_max)]
    return tuple(states[:, 1:])
