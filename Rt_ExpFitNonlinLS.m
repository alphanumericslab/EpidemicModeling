function [Rt, A, Lambda, ExpFit] = Rt_ExpFitNonlinLS(NewCases, wlen, time_unit, varargin)
% Estimates the parameters of an exponential fit using nonlinear least
% squares over the past wlen number of new cases acquired over constant
% time intervals time_unit
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
r = zeros(1, L); % The growth rate (scalar)
if(causal)
    % The amplitude of the exp fit:
    A = filter([zeros(1, wlen-1), 1], 1, NewCases); A = A(:)'; % Just a causal lag to fill in the end points with the raw input
    n = -wlen + 1 : 0; % time sequence of the last wlen samples
    options = optimset('TolX', 1e-6, 'TolFun', 1e-6, 'MaxIter', 250);
    for mm = wlen : L
        segment = NewCases(mm - wlen + 1: mm); % a segment of wlen samples
        if(length(find(segment ~= 0)) < wlen)
            A(mm) = NewCases(mm);
            r(mm) = 0;
        else
            InitParams = [NewCases(mm) 0];
            EstParams = nlinfit(n/time_unit, segment, @ExpModel, InitParams, options);
            %     EstParams = lsqcurvefit(@ExpModel, InitParams, n/time_unit, segment);
            A(mm) = EstParams(1);
            r(mm) = EstParams(2);
        end
    end
else
    A = NewCases(:)'; %zeros(1, L); % The amplitude of the exp fit is equal to the input at its end-points
    wlen_half = floor(wlen/2);
    n = -wlen_half : wlen_half; % time sequence of the wlen_half previous and next samples
    options = optimset('TolX', 1e-6, 'TolFun', 1e-6, 'MaxIter', 250);
    for mm = wlen_half + 1 : L - wlen_half
        segment = NewCases(mm - wlen_half : mm + wlen_half); % a segment of wlen samples
        if(length(find(segment ~= 0)) < wlen)
            A(mm) = NewCases(mm);
            r(mm) = 0;
        else
            InitParams = [NewCases(mm) 0];
            EstParams = nlinfit(n/time_unit, segment, @ExpModel, InitParams, options);
            %     EstParams = lsqcurvefit(@ExpModel, InitParams, n/time_unit, segment);
            A(mm) = EstParams(1);
            r(mm) = EstParams(2);
        end
    end
end

Rt = exp(r); % The reproduction rate
ExpFit = A .* Rt; % The exponential fit
Lambda = r/time_unit; % The reproduction eigenvalue (inverse time unit)

end

function y = ExpModel(params, t)
A = params(1);
lambda = params(2);
y = A*exp(lambda*t);
end
