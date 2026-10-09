function [s, e, i, r, p] = seirp(alpha_e, alpha_i, kappa, rho, beta, mu, gamma, s0, e0, i0, r0, p0, T, dt)
% SEIRP Integrate the five-compartment SEIRP model by forward Euler.
%
% Syntax: [s, e, i, r, p] = seirp(alpha_e, alpha_i, kappa, rho, beta, mu, gamma, s0, e0, i0, r0, p0, T, dt)
% Author: Reza Sameni | Emory University
% Reference: Sameni (2020), arXiv:2003.11371; Sameni (2022),
% doi:10.1109/JSTSP.2021.3129118.
%
% Inputs, outputs, units, shapes and endpoint conventions: docs/api.md.
% All rates use a consistent inverse time unit; dt uses that time unit.
% Outputs are row vectors. SEIRP and SI include their initial conditions;
% SI-alpha returns post-update states and accepts optional 3-by-K noise.
% Scalar or time-dependent transmission schedules are accepted.
% The Open Source Electrophysiological Toolbox, version 3.14, March 2020
% Released under the GNU General Public License
% https://gitlab.com/rsameni/OSET/

assert(isfinite(dt) && dt > 0 && isfinite(T) && T > 0, "Positive finite T and dt required.");
K = round(T/dt);
assert(K >= 1, "At least one sample is required.");
if isscalar(alpha_e), alpha_e = repmat(alpha_e, 1, max(K-1, 1)); end
assert(numel(alpha_e) >= K-1 && all(isfinite(alpha_e)) && all(alpha_e >= 0), "Invalid alpha_e schedule.");
if isscalar(alpha_i), alpha_i = repmat(alpha_i, 1, max(K-1, 1)); end
assert(numel(alpha_i) >= K-1 && all(isfinite(alpha_i)) && all(alpha_i >= 0), "Invalid alpha_i schedule.");
if isscalar(kappa), kappa = repmat(kappa, 1, max(K-1, 1)); end
assert(numel(kappa) >= K-1 && all(isfinite(kappa)) && all(kappa >= 0), "Invalid kappa schedule.");
if isscalar(rho), rho = repmat(rho, 1, max(K-1, 1)); end
assert(numel(rho) >= K-1 && all(isfinite(rho)) && all(rho >= 0), "Invalid rho schedule.");
if isscalar(beta), beta = repmat(beta, 1, max(K-1, 1)); end
assert(numel(beta) >= K-1 && all(isfinite(beta)) && all(beta >= 0), "Invalid beta schedule.");
if isscalar(mu), mu = repmat(mu, 1, max(K-1, 1)); end
assert(numel(mu) >= K-1 && all(isfinite(mu)) && all(mu >= 0), "Invalid mu schedule.");
if isscalar(gamma), gamma = repmat(gamma, 1, max(K-1, 1)); end
assert(numel(gamma) >= K-1 && all(isfinite(gamma)) && all(gamma >= 0), "Invalid gamma schedule.");
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

for t = 1 : K - 1
    s(t + 1) = (-alpha_e(t) * s(t) * e(t) - alpha_i(t) * s(t) * i(t) + gamma(t) * r(t)) * dt + s(t);
    e(t + 1) = (alpha_e(t) * s(t) * e(t) + alpha_i(t) * s(t) * i(t) - kappa(t) * e(t) - rho(t) * e(t)) * dt + e(t);
    i(t + 1) = (kappa(t) * e(t) - beta(t) * i(t) - mu(t) * i(t))* dt + i(t);
    r(t + 1) = (beta(t) * i(t) + rho(t) * e(t) - gamma(t) * r(t)) * dt + r(t);
    p(t + 1) = (mu(t) * i(t)) * dt + p(t);
end


