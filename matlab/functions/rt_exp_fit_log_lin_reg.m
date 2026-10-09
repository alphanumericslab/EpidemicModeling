function [Rt, A, Lambda, ExpFit] = rt_exp_fit_log_lin_reg(NewCases, wlen, time_unit, varargin)
% RT_EXP_FIT_LOG_LIN_REG Fit rolling exponential growth by regression of strictly positive log cases.
%
% Syntax: [Rt, A, Lambda, ExpFit] = rt_exp_fit_log_lin_reg(NewCases, wlen, time_unit, varargin)
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

if(nargin > 3)
    causal = varargin{1};
else
    causal = 1;
end

NewCases = NewCases(:)';
assert(wlen >= 2 && wlen == fix(wlen) && time_unit > 0, "Invalid fitting grid.");
assert(all(isfinite(NewCases)) && all(NewCases >= 0), "Cases must be finite and nonnegative.");
assert(all(NewCases > 0), "Log fitting requires strictly positive cases.");
L = length(NewCases); % The input signal length
NewCasesLog = log(NewCases); % The log of the NewCases numbers
ALog = zeros(1, L); % The log amplitude of the exp fit
r = zeros(1, L); % The growth rate (scalar)
if(isequal(causal, 1))
    n = -wlen + 1 : 0; % time sequence of the last wlen samples
    En = mean(n); % E{n}
    En2 = mean(n.^2); % E{n^2}
    Det = En2 - En^2;
    for mm = wlen : L
        segment = NewCasesLog(mm - wlen + 1: mm); % a segment of wlen samples
        ALog(mm) = (mean(segment)*En2 - mean(n .* segment)*En)/Det;
        r(mm) = (mean(n .* segment) - mean(segment)* En)/Det;
    end
else
    wlen_half = floor(wlen/2);
    n = -wlen_half : wlen_half; % time sequence of the wlen_half previous and next samples
    En = mean(n); % E{n}
    En2 = mean(n.^2); % E{n^2}
    Det = En2 - En^2;
    for mm = wlen_half + 1 : L - wlen_half
        segment = NewCasesLog(mm - wlen_half : mm + wlen_half); % a segment of wlen samples
        ALog(mm) = (mean(segment)*En2 - mean(n .* segment)*En)/Det;
        r(mm) = (mean(n .* segment) - mean(segment)* En)/Det;
    end
end

A = exp(ALog); % The exponential amplitudes
Rt = exp(r); % The reproduction rate
ExpFit = A .* Rt; % The exponential fit
Lambda = r/time_unit; % The reproduction eigenvalue (inverse time unit)
