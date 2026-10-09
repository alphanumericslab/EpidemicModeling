function [u_opt, S_MINUS, S_PLUS, P_MINUS, P_PLUS, K_GAIN, S_SMOOTH, P_SMOOTH, innovations, ...
    rho] = new_case_ekf_estimator_with_optimal_npi(u, x, params, s_init, Ps_init, s_final, ...
    Ps_final, w_bar, v_bar, Q_w, R_v, beta, gamma, inv_monitor_len, order)
    % NEW_CASE_EKF_ESTIMATOR_WITH_OPTIMAL_NPI Run the six-state SI-alpha/costate filter and
    % smoother for code generation.
    %
    % u is interventions-by-time; x contains incidence fractions. The state is
    % [s; i; alpha; costate_s; costate_i; costate_alpha]. Initial/final means
    % have length 6; covariance matrices are 6-by-6. NaN observations skip
    % correction; NaN final entries leave that boundary unconstrained.
    % NaN controls select intervention bounds using the switching function.
    % params contains the fields from default_si_params plus epsilon, w and sigma.
    % Q_w/R_v are noise covariances; w_bar/v_bar are noise means. beta/gamma,
    % inv_monitor_len and order follow generic_extended_kalman_filter.
    % Outputs follow this function signature, which differs from the generic
    % filter: gains precede smoothed states. Uses simple covariance updates and
    % the nonnegative switching tie. See docs/api.md for callback details.
    %
    % Author: Reza Sameni | Emory University

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
    Q = Q_w;
    R = R_v;

    % equal to control input whenever available, otherwise equal to optimal control
    u_opt = zeros(size(u));

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
            [gs, Gsp, gv, Gvp] = obs_hessian_terms(u(:, k), sk_minus, Pk_minus, v_bar, R, params);
        else
            error('Undefined order');
        end

        % Calculate s(k|k) and P(k|k)
        [Ck_minus, Dk_minus] = obs_jacobian(u(:, k), sk_minus, v_bar, params);
        xk_minus = nlin_obs_update(u(:, k), sk_minus, v_bar, params) + gs + gv;

        % Apply hard margins on observations
        xk_minus = obs_hard_margins(xk_minus, params);

        % time update if observation is valid

        if ~isnan(x(:, k))
            innovations(:, k) = x(:, k) - xk_minus;
            Kgain = Pk_minus * Ck_minus' / (Ck_minus * Pk_minus * Ck_minus' + ...
                gamma * (Dk_minus * R * Dk_minus') + Gsp + Gvp); % Kalman gain
            Pk_plus = (eye(m) - Kgain * Ck_minus) * Pk_minus / gamma;
            sk_plus = sk_minus + Kgain * innovations(:, ...
                k);                 % As posteriori state estimate
        else
            innovations(:, k) = 0;
            Kgain = zeros(m, n); % Kalman gain
            Pk_plus = Pk_minus;
            sk_plus = sk_minus; % As posteriori state estimate
        end

        % Apply hard margins on states
        sk_plus = state_hard_margins(sk_plus, params);

        if order == 1
            fs = zeros(m, 1);
            Fsp = zeros(m);
            fw = zeros(m, 1);
            Fwp = zeros(m);
        elseif order == 2
            [fs, Fsp, fw, Fwp] = state_hessian_terms(u(:, k), sk_plus, Pk_plus, w_bar, Q, params);
        else
            error('Undefined order');
        end

        % Calculate s(k+1|k) and P(k+1|k) for k+1
        [u_opt(:, k), sk_minus] = nlin_state_update(u(:, k), sk_plus, w_bar, params); % State update
        sk_minus = sk_minus + fs + fw; % Add second order terms (is available)
        [Ak_plus, Bk_plus] = state_jacobians(u(:, k), sk_plus, w_bar, params);
        Pk_minus = (Ak_plus * Pk_plus * Ak_plus') + (Bk_plus * Q * Bk_plus') + Fsp + ...
            Fwp; % Cov. matrix update

        % Apply hard margins on states
        sk_minus = state_hard_margins(sk_minus, params);

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

        InnovationsCovNormalized = cat(3, cc / R, InnovationsCovNormalized(:, :, ...
            1:inv_monitor_len - 1));
        rho(:, :, k) = sum(InnovationsCovNormalized, 3) / stats_counter;

        if ~isequal(beta, 1) && ~isnan(x(:, k))
            R = beta * R + (1 - beta) * sum(InnovationsCov, 3) / stats_counter;
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

    for kk = 1:numel(row)
        P_SMOOTH(row(kk), col(kk), T) = Ps_final(row(kk), col(kk));
    end

    for k = T - 1:-1:1
        sk_plus = S_PLUS(:, k);
        Ak_plus = state_jacobians(u(:, k), sk_plus, w_bar, params);
        J = (P_PLUS(:, :, k) * Ak_plus') * pinv(P_MINUS(:, :, k + 1));
        S_SMOOTH(:, k) = S_PLUS(:, k) + J * (S_SMOOTH(:, k + 1) - S_MINUS(:, k + 1));

        % Apply hard margins on states
        S_SMOOTH(:, k) = state_hard_margins(S_SMOOTH(:, k), params);

        P_SMOOTH(:, :, k) = P_PLUS(:, :, k) - J * (P_MINUS(:, :, k + 1) - P_SMOOTH(:, :, ...
            k + 1)) * J';
    end

    % Squeeze excess dimensions if applicable
    rho = squeeze(rho);
end
