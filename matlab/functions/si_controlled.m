function [s, i] = si_controlled(alpha, beta, s0, i0, K, dt)
    % SI_CONTROLLED Integrate susceptible and infected fractions under a given contact schedule.
    %
    % alpha is a scalar or a vector of at least K-1 contact rates; beta is the
    % removal rate. s0/i0 are initial fractions; K is the number of samples and
    % dt is the time step. Rates and dt use consistent inverse-time/time units.
    % Returns two 1-by-K row vectors including the initial fractions.
    % Each Euler update is clipped to [0,1], following the original model.
    %
    % Author: Reza Sameni | Emory University
    % References: Sameni (2020), arXiv:2003.11371;
    % Sameni (2022), doi:10.1109/JSTSP.2021.3129118.
    % The Open Source Electrophysiological Toolbox, version 3.14, January 2021
    % Released under the GNU General Public License

    assert(K >= 1 && K == fix(K) && dt > 0, "Positive integer K and dt required.");

    if isscalar(alpha)
        alpha = repmat(alpha, 1, max(K - 1, 1));
    end

    assert(numel(alpha) >= K - 1 && all(isfinite(alpha)) && all(alpha >= 0), ...
        "Invalid alpha schedule.");
    s = zeros(1, K);
    i = zeros(1, K);

    s(1) = s0;
    i(1) = i0;

    % State equations

    for t = 1:K - 1
        s(t + 1) = max(0.0, min(1.0, s(t) - dt * alpha(t) * s(t) * i(t)));
        i(t + 1) = max(0.0, min(1.0, i(t) + dt * (alpha(t) * s(t) * i(t) - beta * i(t))));
    end
