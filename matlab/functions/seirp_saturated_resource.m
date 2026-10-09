function [s, e, i, r, p] = seirp_saturated_resource(alpha_e, alpha_i, kappa, rho, gamma, s0, e0, ...
    i0, r0, p0, T, dt, beta_0, beta_s, mu_0, mu_s, sigma, i_0)
    % SEIRP_SATURATED_RESOURCE Integrate SEIRP with recovery and mortality changing as
    % healthcare saturates.
    %
    % alpha_e/alpha_i, kappa, rho, gamma and initial fractions follow seirp.
    % T is the duration and dt is the step size, in one consistent time unit.
    % beta_0/mu_0 are recovery/mortality rates before saturation; beta_s/mu_s
    % are their saturated values. i_0 is the infected-fraction threshold and
    % sigma > 0 sets the transition width. Rates use inverse time units.
    % Returns five 1-by-round(T/dt) row vectors, including the initial state.
    % The transition is h = (tanh((i-i_0)/sigma)+1)/2. Reduce dt if fractions
    % become negative; this Euler solver does not clip the output.
    %
    % Author: Reza Sameni | Emory University
    % References: Sameni (2020), arXiv:2003.11371;
    % Sameni (2022), doi:10.1109/JSTSP.2021.3129118.
    % The Open Source Electrophysiological Toolbox, version 3.14, March 2020
    % Released under the GNU General Public License
    % https://gitlab.com/rsameni/OSET/

    assert(isfinite(dt) && dt > 0 && isfinite(T) && T > 0, "Positive finite T and dt required.");

    assert(sigma > 0, "sigma must be positive.");
    K = round(T / dt);

    assert(K >= 1, "At least one sample is required.");

    if isscalar(alpha_e)
        alpha_e = repmat(alpha_e, 1, max(K - 1, 1));
    end

    assert(numel(alpha_e) >= K - 1 && all(isfinite(alpha_e)) && all(alpha_e >= 0), ...
        "Invalid alpha_e schedule.");

    if isscalar(alpha_i)
        alpha_i = repmat(alpha_i, 1, max(K - 1, 1));
    end

    assert(numel(alpha_i) >= K - 1 && all(isfinite(alpha_i)) && all(alpha_i >= 0), ...
        "Invalid alpha_i schedule.");

    if isscalar(kappa)
        kappa = repmat(kappa, 1, max(K - 1, 1));
    end

    assert(numel(kappa) >= K - 1 && all(isfinite(kappa)) && all(kappa >= 0), ...
        "Invalid kappa schedule.");

    if isscalar(rho)
        rho = repmat(rho, 1, max(K - 1, 1));
    end

    assert(numel(rho) >= K - 1 && all(isfinite(rho)) && all(rho >= 0), "Invalid rho schedule.");

    if isscalar(gamma)
        gamma = repmat(gamma, 1, max(K - 1, 1));
    end

    assert(numel(gamma) >= K - 1 && all(isfinite(gamma)) && all(gamma >= 0), ...
        "Invalid gamma schedule.");
    s = zeros(1, K);
    e = zeros(1, K);
    i = zeros(1, K);
    r = zeros(1, K);

    p = zeros(1, K);

    s(1) = s0;
    e(1) = e0;
    i(1) = i0;
    r(1) = r0;
    p(1) = p0;

    for t = 1:K - 1
        h = (tanh((i(t) - i_0) / sigma) + 1) / 2;
        beta = (beta_s - beta_0) * h + beta_0;
        mu = (mu_s - mu_0) * h + mu_0;

        s(t + 1) = (-alpha_e(t) * s(t) * e(t) - alpha_i(t) * s(t) * i(t) + gamma(t) * ...
            r(t)) * dt + s(t);
        e(t + 1) = (alpha_e(t) * s(t) * e(t) + alpha_i(t) * s(t) * i(t) - kappa(t) * ...
            e(t) - rho(t) * e(t)) * dt + e(t);
        i(t + 1) = (kappa(t) * e(t) - beta * i(t) - mu * i(t)) * dt + i(t);
        r(t + 1) = (beta * i(t) + rho(t) * e(t) - gamma(t) * r(t)) * dt + r(t);
        p(t + 1) = (mu * i(t)) * dt + p(t);
    end
