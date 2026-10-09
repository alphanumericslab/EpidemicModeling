function [Rt, Lambda, RtSmoothed, LambdaSmoothed] = rt_exp_fit_gen_ratios(NewCases, wlen, ...
    generation_period, time_unit)
    % RT_EXP_FIT_GEN_RATIOS Estimate growth from case counts separated by a fixed lag.
    %
    % NewCases is a nonnegative length-N vector. generation_period is an integer
    % lag in samples; wlen >= 2 is the trailing averaging window in samples.
    % time_unit > 0 is the interval used in Rt = exp(Lambda*time_unit).
    % Returns four 1-by-N vectors: growth factor Rt, log growth Lambda, and
    % their smoothed versions. Lambda starts with generation_period zeros.
    % Smoothing uses a causal zero-padded average. Zero case counts can produce
    % NaN or Inf. Rt is the historical growth factor, not a mechanistic R0.
    %
    % Author: Reza Sameni | Emory University
    % References: Sameni (2020), arXiv:2003.11371;
    % Sameni (2022), doi:10.1109/JSTSP.2021.3129118.

    NewCases = NewCases(:)';

    assert(wlen >= 2 && wlen == fix(wlen) && time_unit > 0, "Invalid fitting grid.");

    assert(all(isfinite(NewCases)) && all(NewCases >= 0), "Cases must be finite and nonnegative.");

    assert(generation_period >= 1 && generation_period <= numel(NewCases) && ...
        generation_period == fix(generation_period), "Invalid generation period.");
    Lambda = [zeros(1, generation_period), log(NewCases(1 + generation_period:end) ./ ...
        NewCases(1:end - generation_period))] / ...
        generation_period; % The reproduction eigenvalue (inverse time unit)
    LambdaSmoothed = filter(ones(1, wlen), wlen, Lambda);

    Rt = exp(Lambda * time_unit); % The reproduction rate
    RtSmoothed = exp(LambdaSmoothed * time_unit); % The smoothed reproduction rate
