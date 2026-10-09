function [daily, smoothed] = prepare_cases(cumulative, window)
    % PREPARE_CASES Clean cumulative counts and apply a causal moving average.
    % [daily,smoothed] = prepare_cases(cumulative,window)
    % cumulative: cumulative counts; window: positive integer (default 7).
    % Outputs are row vectors. First difference is zero; negative revisions and
    % nonfinite differences become zero. Smoothing has zero-padded startup.
    % Author: Reza Sameni | Emory University

    if nargin < 2
        window = 7;
    end

    assert(window >= 1 && window == fix(window), 'Invalid window.');
    daily = [0 diff(cumulative(:)')];
    daily(~isfinite(daily)) = 0;
    daily = max(0, daily);
    smoothed = filter(ones(1, window) / window, 1, daily);
end
