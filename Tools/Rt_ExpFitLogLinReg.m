function [Rt, A, Lambda, ExpFit] = Rt_ExpFitLogLinReg(NewCases, wlen, time_unit, varargin)
% Estimates the parameters of an exponential fit. The function performs a
% linear regression over the log of the past wlen number of new cases
% acquired over constant time intervals time_unit
%
% Reza Sameni
% Dec 2020
% Email: reza.sameni@gmail.com

if(nargin > 3)
    causal = varargin{1};
else
    causal = 1;
end

NewCases = NewCases(:)';
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
