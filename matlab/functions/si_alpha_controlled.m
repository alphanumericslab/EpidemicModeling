function [s, i, alpha] = si_alpha_controlled(u, s0, i0, alpha0, u_max, alpha_min, alpha_max, ...
    gamma, a, b, beta, s_noise_std, i_noise_std, alpha_noise_std, K, dt, noise)
    % SI_ALPHA_CONTROLLED Integrate susceptible fraction, infected fraction, and contact rate.
    %
    % u is interventions-by-K; a and u_max have one entry per intervention.
    % s0/i0/alpha0 are initial values. gamma sets the contact-response rate;
    % b + a'*(u_max-u) is the target contact rate. beta is the removal rate.
    % alpha_min/alpha_max bound the contact rate; dt is the time step.
    % s_noise_std/i_noise_std/alpha_noise_std scale additive process noise.
    % An optional 3-by-K noise array supplies standard-normal draws; otherwise
    % MATLAB generates them. Supply the same array in Python for exact parity.
    % Returns three 1-by-K vectors AFTER each update, excluding the initial state.
    % Population fractions and contact rate are clipped to their bounds.
    %
    % Author: Reza Sameni | Emory University
    % References: Sameni (2020), arXiv:2003.11371;
    % Sameni (2022), doi:10.1109/JSTSP.2021.3129118.
    % The Open Source Electrophysiological Toolbox, version 3.14, January 2021
    % Released under the GNU General Public License

    if nargin < 17
        noise = randn(3, K);
    end

    assert(isequal(size(noise), [3 K]), "noise must be 3-by-K.");

    assert(isequal(size(u), [numel(a) K]) && numel(u_max) == numel(a) && dt > 0, "Invalid inputs.");
    s = zeros(1, K + 1);
    i = zeros(1, K + 1);
    alpha = zeros(1, K + 1);

    s(1) = s0;
    i(1) = i0;
    alpha(1) = alpha0;

    % State equations

    for t = 1:K
        s(t + 1) = max(0.0, min(1.0, s(t) - dt * (alpha(t) * s(t) * i(t) + noise(1, t) * ...
            s_noise_std)));
        i(t + 1) = max(0.0, min(1.0, i(t) + dt * (alpha(t) * s(t) * i(t) - beta * i(t) + ...
            noise(2, t) * i_noise_std)));
        alpha(t + 1) = max(alpha_min, min(alpha_max, alpha(t) + dt * (-gamma * ...
            alpha(t) + gamma * b + gamma * a' * (u_max - u(:, t)) + noise(3, t) * alpha_noise_std)));
    end

    s = s(2:end);
    i = i(2:end);
    alpha = alpha(2:end);
