function [Rt, Lambda, RtSmoothed, LambdaSmoothed] = rt_exp_fit_gen_ratios(NewCases, wlen, generation_period, time_unit)
% RT_EXP_FIT_GEN_RATIOS Estimate growth using lagged ratios and causal moving averages.
%
% Syntax: [Rt, Lambda, RtSmoothed, LambdaSmoothed] = rt_exp_fit_gen_ratios(NewCases, wlen, generation_period, time_unit)
% Author: Reza Sameni | Emory University
% Reference: Sameni (2020), arXiv:2003.11371; Sameni (2022),
% doi:10.1109/JSTSP.2021.3129118.
%
% Inputs, outputs, units, shapes and endpoint conventions: docs/api.md.
% All rates use a consistent inverse time unit; dt uses that time unit.
% NewCases is a nonnegative count vector; wlen is an integer >=2;
% time_unit is positive. Outputs are row vectors with input length.
% rt is a legacy growth factor; ExpFit is a one-step forecast.
% See docs/api.md for zero handling and nonlinear time scaling.

NewCases = NewCases(:)';
assert(wlen >= 2 && wlen == fix(wlen) && time_unit > 0, "Invalid fitting grid.");
assert(all(isfinite(NewCases)) && all(NewCases >= 0), "Cases must be finite and nonnegative.");
assert(generation_period >= 1 && generation_period <= numel(NewCases) && generation_period == fix(generation_period), "Invalid generation period.");
Lambda = [zeros(1, generation_period) , log(NewCases(1 + generation_period : end)./NewCases(1 : end - generation_period))]/generation_period; % The reproduction eigenvalue (inverse time unit)
LambdaSmoothed = filter(ones(1, wlen), wlen, Lambda);

Rt = exp(Lambda * time_unit); % The reproduction rate
RtSmoothed = exp(LambdaSmoothed * time_unit); % The smoothed reproduction rate


