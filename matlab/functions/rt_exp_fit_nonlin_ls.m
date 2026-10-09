function [Rt, A, Lambda, ExpFit] = rt_exp_fit_nonlin_ls(NewCases, wlen, time_unit, varargin)
    % RT_EXP_FIT_NONLIN_LS Fit exponential curves directly to rolling windows of case counts.
    %
    % NewCases is a nonnegative length-N vector. wlen >= 2, time_unit > 0 and
    % the optional fourth argument causal follow rt_exp_fit_log_lin_reg.
    % causal defaults to true; false selects centered windows.
    % Returns 1-by-N vectors Rt, A, Lambda, ExpFit. The historical fitting grid
    % uses sample offsets/time_unit, with Lambda = fitted_rate/time_unit.
    % Windows containing zeros use the current count and zero growth.
    % Unfitted amplitudes use delayed counts for causal windows and raw counts
    % for centered windows. Requires nlinfit from Statistics and Machine Learning Toolbox.
    %
    % Author: Reza Sameni | Emory University
    % References: Sameni (2020), arXiv:2003.11371;
    % Sameni (2022), doi:10.1109/JSTSP.2021.3129118.

    if nargin > 3
        causal = varargin{1};
    else
        causal = 1;
    end

    NewCases = NewCases(:)';

    assert(wlen >= 2 && wlen == fix(wlen) && time_unit > 0, "Invalid fitting grid.");

    assert(all(isfinite(NewCases)) && all(NewCases >= 0), "Cases must be finite and nonnegative.");
    L = length(NewCases); % The input signal length
    r = zeros(1, L); % The growth rate (scalar)

    if causal
        % The amplitude of the exp fit:
        A = filter([zeros(1, wlen - 1), 1], 1, NewCases);
        A = A(:)'; % Just a causal lag to fill in the end points with the raw input
        n = -wlen + 1:0; % time sequence of the last wlen samples
        options = optimset('TolX', 1e-6, 'TolFun', 1e-6, 'MaxIter', 250);

        for mm = wlen:L
            segment = NewCases(mm - wlen + 1:mm); % a segment of wlen samples

            if any(segment == 0)
                A(mm) = NewCases(mm);
                r(mm) = 0;
            else
                InitParams = [NewCases(mm) 0];
                EstParams = nlinfit(n / time_unit, segment, @exp_model, InitParams, options);
                A(mm) = EstParams(1);
                r(mm) = EstParams(2);
            end
        end
    else
        A = NewCases(:)'; % zeros(1, L); % The amplitude of the exp fit is equal to the input at its end-points
        wlen_half = floor(wlen / 2);
        n = -wlen_half:wlen_half; % time sequence of the wlen_half previous and next samples
        options = optimset('TolX', 1e-6, 'TolFun', 1e-6, 'MaxIter', 250);

        for mm = wlen_half + 1:L - wlen_half
            segment = NewCases(mm - wlen_half:mm + wlen_half); % a segment of wlen samples

            if any(segment == 0)
                A(mm) = NewCases(mm);
                r(mm) = 0;
            else
                InitParams = [NewCases(mm) 0];
                EstParams = nlinfit(n / time_unit, segment, @exp_model, InitParams, options);
                A(mm) = EstParams(1);
                r(mm) = EstParams(2);
            end
        end
    end

    Rt = exp(r); % The reproduction rate
    ExpFit = A .* Rt; % The exponential fit
    Lambda = r / time_unit; % The reproduction eigenvalue (inverse time unit)

end

function y = exp_model(params, t)
    % EXP_MODEL Evaluate amplitude times the exponential of growth times time.
    A = params(1);
    lambda = params(2);
    y = A * exp(lambda * t);
end
