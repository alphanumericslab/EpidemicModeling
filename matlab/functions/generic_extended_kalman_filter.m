function [u_opt, u_opt_smooth, S_MINUS, S_PLUS, S_SMOOTH, P_MINUS, P_PLUS, P_SMOOTH, K_GAIN, ...
    innovations, rho] = generic_extended_kalman_filter(u, x, handles, params, s_init, Ps_init, ...
    s_final, Ps_final, w_bar, v_bar, Q_w, R_v, beta, gamma, inv_monitor_len, order)
    % GENERIC_EXTENDED_KALMAN_FILTER Run nonlinear filtering and fixed-interval smoothing
    % using model callbacks.
    %
    % handles supplies state/observation updates, bounds, Jacobians and Hessian
    % terms; callback signatures are documented in docs/api.md.
    % u and x are controls-by-time and observations-by-time. Their units
    % are defined by the callbacks. A missing observation skips correction. params is passed unchanged to each model callback.
    % s_init/Ps_init are initial mean/covariance; s_final/Ps_final set the
    % smoothing boundary. NaN final entries leave that boundary unconstrained.
    % w_bar/v_bar are noise means; Q_w/R_v are fixed covariance matrices or
    % covariance stacks with time on the last axis. beta in [0,1] controls
    % measurement-noise adaptation; gamma in (0,1] stabilizes covariance.
    % inv_monitor_len is a positive window length; order selects 1 or 2.
    % Outputs *_MINUS, *_PLUS, *_SMOOTH are predicted, corrected, and smoothed
    % states/covariances. State columns and the last covariance axis are time.
    % K_GAIN is the correction gain; innovations are observation residuals;
    % rho monitors normalized innovation covariance. u_opt/u_opt_smooth store
    % supplied or inferred controls. See docs/api.md for exact output order.
    %
    % Author: Reza Sameni | Emory University
    % References: Sameni (2020), arXiv:2003.11371;
    % Sameni (2022), doi:10.1109/JSTSP.2021.3129118.
    % (C) Reza Sameni, 2021

    T = size(x, 2); % number of time samples
    n = size(x, 1); % number of observations
    m = length(s_init); % number of state variables
    % l = length(w_bar); % number of process noises
    % p = length(v_bar); % number of observation noises

    % //////////////////////////////////////////////////////////////////////////
    S_MINUS = zeros(m, T);
    S_PLUS = zeros(m, T);
    P_MINUS = zeros(m, m, T);
    P_PLUS = zeros(m, m, T);
    K_GAIN = zeros(m, n, T);

    innovations = zeros(n, T);
    rho = zeros(n, n, T);
    InnovationsMean = zeros(n, inv_monitor_len);
    InnovationsCovNormalized = zeros(n, n, inv_monitor_len);
    InnovationsCov = zeros(n, n, inv_monitor_len);

    % Initialization
    sk_minus = s_init(:);
    Pk_minus = Ps_init;

    % Normalize fixed and time-dependent covariance inputs.
    [Q, ~] = covariance_series(Q_w, T);
    [R, fixed_R] = covariance_series(R_v, T);

    assert(order == 1 || order == 2, "order must be 1 or 2.");

    assert(beta >= 0 && beta <= 1 && gamma > 0 && gamma <= 1, "Invalid filter factors.");

    assert(inv_monitor_len >= 1 && inv_monitor_len == fix(inv_monitor_len), ...
        "Invalid monitor length.");

    % equal to control input whenever available, otherwise equal to optimal control
    u_opt = zeros(size(u));
    u_opt_smooth = zeros(size(u));

    % Forward Kalman Filtering Stage

    for k = 1:T
        % Store results of s_minus from previous iteration
        S_MINUS(:, k) = sk_minus;
        P_MINUS(:, :, k) = Pk_minus;

        if order == 1
            gs = zeros(n, 1);
            Gsp = zeros(n);
            gv = zeros(n, 1);
            Gvp = zeros(n);
        elseif order == 2
            [gs, Gsp, gv, Gvp] = handles.obs_hessian_terms(u(:, k), sk_minus, ...
                Pk_minus, v_bar, R(:, :, k), params, k);
        else
            error('Undefined order');
        end

        % Calculate s(k|k) and P(k|k)
        [Ck_minus, Dk_minus] = handles.obs_jacobian(u(:, k), sk_minus, v_bar, params, k);
        xk_minus = handles.nlin_obs_update(u(:, k), sk_minus, v_bar, params, k) + gs + gv;

        % Apply hard margins on observations
        xk_minus = handles.obs_hard_margins(xk_minus, params, k);

        % time update if observation is valid

        if all(isfinite(x(:, k)))
            innovations(:, k) = x(:, k) - xk_minus;
            Kgain = Pk_minus * Ck_minus' / (Ck_minus * Pk_minus * Ck_minus' + ...
                gamma * (Dk_minus * R(:, :, k) * Dk_minus') + Gsp + Gvp); % Kalman gain

            Pk_plus = ((eye(m) - Kgain * Ck_minus) * Pk_minus * (eye(m) - Kgain * ...
                Ck_minus)' + Kgain * (Dk_minus * R(:, :, k) * Dk_minus') * Kgain') / ...
                gamma; % Stabilized Kalman cov. matrix

            sk_plus = sk_minus + Kgain * innovations(:, k); % As posteriori state estimate
        else
            innovations(:, k) = 0;
            Kgain = zeros(m, n); % Kalman gain
            Pk_plus = Pk_minus;
            sk_plus = sk_minus; % As posteriori state estimate
        end

        % Trivial condition added to guarantee numerical stability
        Pk_plus = (Pk_plus + Pk_plus') / 2.0;

        % Apply hard margins on states
        sk_plus = handles.state_hard_margins(sk_plus, params, k);

        if order == 1
            fs = zeros(m, 1);
            Fsp = zeros(m);
            fw = zeros(m, 1);
            Fwp = zeros(m);
        elseif order == 2
            [fs, Fsp, fw, Fwp] = handles.state_hessian_terms(u(:, k), sk_plus, ...
                Pk_plus, w_bar, Q(:, :, k), params, k);
        else
            error('Undefined order');
        end

        % Calculate s(k+1|k) and P(k+1|k) for k+1
        [u_opt(:, k), sk_minus] = handles.nlin_state_update(u(:, k), sk_plus, w_bar, ...
            params, k); % State update
        sk_minus = sk_minus + fs + fw; % Add second order terms (is available)
        [Ak_plus, Bk_plus] = handles.state_jacobians(u(:, k), sk_plus, w_bar, params, k);
        Pk_minus = (Ak_plus * Pk_plus * Ak_plus') + (Bk_plus * Q(:, :, k) * Bk_plus') + ...
            Fsp + Fwp; % Cov. matrix update

        % Trivial condition added to guarantee numerical stability
        Pk_minus = (Pk_minus + Pk_minus') / 2.0;

        % Apply hard margins on states
        sk_minus = handles.state_hard_margins(sk_minus, params, k);

        % Store results of s_plus
        S_PLUS(:, k) = sk_plus;
        P_PLUS(:, :, k) = Pk_plus;
        K_GAIN(:, :, k) = Kgain;

        % Monitoring the innovation variance and update the observation noise
        stats_counter = min(k, inv_monitor_len);
        InnovationsMean = cat(2, innovations(:, k), InnovationsMean(:, 1:inv_monitor_len - 1));
        mu_k = sum(InnovationsMean, 2) / stats_counter;
        cc = (innovations(:, k) - mu_k) * (innovations(:, k) - mu_k)'; % with mean cancellation
        InnovationsCov = cat(3, cc, InnovationsCov(:, :, 1:inv_monitor_len - 1));

        InnovationsCovNormalized = cat(3, cc / (R(:, :, k) + eps), ...
            InnovationsCovNormalized(:, :, 1:inv_monitor_len - 1));
        rho(:, :, k) = sum(InnovationsCovNormalized, 3) / stats_counter;

        if ~isequal(beta, 1) && all(isfinite(x(:, k))) && fixed_R == true && k < T
            R_estim = sum(InnovationsCov, 3) / stats_counter;
            R(:, :, k + 1) = beta * R(:, :, k) + (1 - beta) * R_estim;
        end
    end

    % Backward Kalman Smoothing Stage
    S_SMOOTH = zeros(size(S_PLUS));
    S_SMOOTH(:, T) = S_PLUS(:, T);
    P_SMOOTH = zeros(size(P_PLUS));
    P_SMOOTH(:, :, T) = P_PLUS(:, :, T);

    % Replace estimates with boundary conditions, if available
    fixed_end_state = find(~isnan(s_final));
    S_SMOOTH(fixed_end_state, T) = s_final(fixed_end_state);

    fixed_end_covs = find(~isnan(Ps_final));
    [row, col] = ind2sub(size(Ps_final), fixed_end_covs);

    for kk = 1:length(row)
        P_SMOOTH(row(kk), col(kk), T) = Ps_final(row(kk), col(kk));
    end

    for k = T - 1:-1:1
        sk_plus = S_PLUS(:, k);
        Ak_plus = handles.state_jacobians(u(:, k), sk_plus, w_bar, params, k);

        % Check to make sure that P_MINUS is not ill-conditioned
        pmns = P_MINUS(:, :, k + 1);

        if sum(isnan(pmns(:))) > 0 || sum(isinf(pmns(:))) > 0 % || rcond(pmns) < 2.0*eps)
            J = zeros(m);
        else
            J = (P_PLUS(:, :, k) * Ak_plus') * pinv(P_MINUS(:, :, k + 1));
        end

        S_SMOOTH(:, k) = S_PLUS(:, k) + J * (S_SMOOTH(:, k + 1) - S_MINUS(:, k + 1));

        % Apply hard margins on states
        S_SMOOTH(:, k) = handles.state_hard_margins(S_SMOOTH(:, k), params, k);

        P_SMOOTH(:, :, k) = P_PLUS(:, :, k) - J * (P_MINUS(:, :, k + 1) - P_SMOOTH(:, :, ...
            k + 1)) * J';

        % Trivial condition added to guarantee numerical stability
        P_SMOOTH(:, :, k) = (P_SMOOTH(:, :, k) + P_SMOOTH(:, :, k)') / 2.0;

        % rerun the state equation to find the optimal input
        [u_opt_smooth(:, k), ~] = handles.nlin_state_update(u(:, k), S_SMOOTH(:, k), ...
            w_bar, params, k);
    end

    % Squeeze excess dimensions if applicable
    rho = squeeze(rho);
end

function [stack, fixed] = covariance_series(value, count)
    % COVARIANCE_SERIES Normalize fixed matrices and time-dependent covariances.

    if ndims(value) == 3 && size(value, 1) == size(value, 2) && size(value, 3) == count
        stack = value;
        fixed = false;
    elseif ismatrix(value) && size(value, 1) == size(value, 2)
        stack = repmat(value, 1, 1, count);
        fixed = true;
    elseif isvector(value) && numel(value) == count
        stack = reshape(value, 1, 1, count);
        fixed = false;
    else
        error('Covariance must be square or square-by-time.');
    end
end
