function [s, i] = si_controlled(alpha, beta, s0, i0, K, dt)
% SI_CONTROLLED Integrate bounded SI fractions with a supplied contact-rate schedule.
%
% Syntax: [s, i] = si_controlled(alpha, beta, s0, i0, K, dt)
% Author: Reza Sameni | Emory University
% Reference: Sameni (2020), arXiv:2003.11371; Sameni (2022),
% doi:10.1109/JSTSP.2021.3129118.
%
% Inputs, outputs, units, shapes and endpoint conventions: docs/api.md.
% All rates use a consistent inverse time unit; dt uses that time unit.
% Outputs are row vectors. SEIRP and SI include their initial conditions;
% SI-alpha returns post-update states and accepts optional 3-by-K noise.
% Scalar or time-dependent transmission schedules are accepted.
% The Open Source Electrophysiological Toolbox, version 3.14, January 2021
% Released under the GNU General Public License

assert(K >= 1 && K == fix(K) && dt > 0, "Positive integer K and dt required.");
if isscalar(alpha), alpha = repmat(alpha, 1, max(K-1, 1)); end
assert(numel(alpha) >= K-1 && all(isfinite(alpha)) && all(alpha >= 0), "Invalid alpha schedule.");
s = zeros(1, K);
i = zeros(1, K);

s(1) = s0;
i(1) = i0;

% State equations
for t = 1 : K - 1
    s(t + 1) = max(0.0, min(1.0, s(t) - dt * alpha(t) * s(t) * i(t)));
    i(t + 1) = max(0.0, min(1.0, i(t) + dt * (alpha(t) * s(t) * i(t) - beta * i(t))));
end
