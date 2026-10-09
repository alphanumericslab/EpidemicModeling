function [u_opt, u_opt_smooth, S_MINUS, S_PLUS, S_SMOOTH, P_MINUS, P_PLUS, P_SMOOTH, K_GAIN, ...
    innovations, rho] = si_alpha_model_backward_ekf(u, x, params, s_init, Ps_init, s_final, ...
    Ps_final, w_bar, v_bar, Q_w, R_v, beta, gamma, inv_monitor_len, order)
    % SI_ALPHA_MODEL_BACKWARD_EKF Filter three SI-alpha states and smooth their trajectories.
    %
    % The first three state rows are susceptible fraction, infected fraction,
    % and contact rate. Means have length 3; covariances are 3-by-3.
    % u and x are controls-by-time and observations-by-time. Observations use
    % population fractions, not raw case counts. A missing observation skips
    % correction. params.obs_type selects NEWCASES (s*i*alpha) or TOTALCASES (1-s).
    % params contains dt, beta, gamma, a, b, u_min/u_max, alpha_min/alpha_max
    % and s_min/i_min; default_si_params creates these fields.
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
    % Supply chronological samples. The initial boundary starts at the latest
    % sample; the final boundary applies at the earliest. Outputs are chronological.
    %
    % Author: Reza Sameni | Emory University
    % References: Sameni (2020), arXiv:2003.11371;
    % Sameni (2022), doi:10.1109/JSTSP.2021.3129118.

    handles.state_hard_margins = @state_hard_margins_flipped;
    handles.obs_hard_margins = @obs_hard_margins_flipped;
    handles.nlin_state_update = @nlin_state_update_flipped;
    handles.nlin_obs_update = @nlin_obs_update_flipped;
    handles.state_jacobians = @state_jacobians_flipped;

    handles.obs_jacobian = @obs_jacobian_flipped;
    handles.state_hessian_terms = @state_hessian_terms_flipped;
    handles.obs_hessian_terms = @obs_hessian_terms_flipped;

    % time flip the samples
    u_flipped = u(:, end:-1:1);
    x_flipped = x(end:-1:1);
    s_init_flipped = s_final;
    s_final_flipped = s_init;
    Ps_init_flipped = Ps_final;

    Ps_final_flipped = Ps_init;

    if ndims(Q_w) == 3
        Q_w = flip(Q_w, 3);
    elseif isvector(Q_w) && ~isscalar(Q_w)
        Q_w = flip(Q_w);
    end

    if ndims(R_v) == 3
        R_v = flip(R_v, 3);
    elseif isvector(R_v) && ~isscalar(R_v)
        R_v = flip(R_v);
    end

    % Call the generic EKF function with the desired function handles
    [u_opt_flipped, u_opt_smooth_flipped, S_MINUS_flipped, S_PLUS_flipped, S_SMOOTH_flipped, ...
        P_MINUS_flipped, P_PLUS_flipped, P_SMOOTH_flipped, K_GAIN_flipped, innovations_flipped, ...
        rho_flipped] = generic_extended_kalman_filter(u_flipped, x_flipped, handles, params, ...
        s_init_flipped, Ps_init_flipped, s_final_flipped, Ps_final_flipped, w_bar, v_bar, Q_w, ...
        R_v, beta, gamma, inv_monitor_len, order);

    % flip back the results in time
    u_opt = u_opt_flipped(:, end:-1:1);
    u_opt_smooth = u_opt_smooth_flipped(:, end:-1:1);
    S_MINUS = S_MINUS_flipped(:, end:-1:1);
    S_PLUS = S_PLUS_flipped(:, end:-1:1);
    S_SMOOTH = S_SMOOTH_flipped(:, end:-1:1);

    P_MINUS = P_MINUS_flipped(:, :, end:-1:1);
    P_PLUS = P_PLUS_flipped(:, :, end:-1:1);
    P_SMOOTH = P_SMOOTH_flipped(:, :, end:-1:1);
    K_GAIN = K_GAIN_flipped(:, :, end:-1:1);
    innovations = innovations_flipped(:, end:-1:1);

    if isvector(rho_flipped)
        rho = flip(rho_flipped);
    else
        rho = flip(rho_flipped, ndims(rho_flipped));
    end

end

% THE SYSTEM EQUATIONS

% Hard margins on state vectors
function s_k = state_hard_margins_flipped(s_k, params, k)
    % STATE_HARD_MARGINS_FLIPPED Constrain population fractions and contact rate; preserve costates.
    s_k(1) = min(1.0, max(0, s_k(1)));
    s_k(2) = min(1.0, max(0, s_k(2)));
    s_k(3) = min(params.alpha_max, max(params.alpha_min, s_k(3)));
end

% Hard margins on observations
function x_k = obs_hard_margins_flipped(x_k, params, k)
    % OBS_HARD_MARGINS_FLIPPED Constrain predicted observations to nonnegative values.
    x_k = max(0, x_k);
end

% Nonlinear state update
function [u, s_k_plus_one] = nlin_state_update_flipped(u, s_k, w_bar, params, k)
    % NLIN_STATE_UPDATE_FLIPPED Evaluate one Euler transition and fill missing controls when
    % applicable.

    s_k_plus_one = zeros(3, 1);

    % State equations
    s_k_plus_one(1) = max(0.0, min(1.0, s_k(1) + params.dt * s_k(3) * s_k(1) * s_k(2)));
    s_k_plus_one(2) = max(0.0, min(1.0, s_k(2) - params.dt * (s_k(3) * s_k(1) * s_k(2) - ...
        params.beta * s_k(2))));
    s_k_plus_one(3) = max(params.alpha_min, min(params.alpha_max, s_k(3) - params.dt * (- ...
        params.gamma * s_k(3) + params.gamma * params.b + params.gamma * params.a' * ...
        (params.u_max - u))));

end

% Nonlinear observation update
function x_k = nlin_obs_update_flipped(u, s_k, v_bar, params, k)
    % NLIN_OBS_UPDATE_FLIPPED Evaluate the configured new-case or total-case observation.

    if isequal(params.obs_type, 'NEWCASES')
        x_k = s_k(1) * s_k(2) * s_k(3) + v_bar;
    elseif isequal(params.obs_type, 'TOTALCASES')
        x_k = 1 - s_k(1) + v_bar; % following the revised model that takes total cases as input
    else
        error('unknown observation type');
    end
end

% State equation Jacobian
function [A, B] = state_jacobians_flipped(u, s_k, w_bar, params, k)
    % STATE_JACOBIANS_FLIPPED Evaluate analytic state-transition and process-noise Jacobians.

    A = zeros(3);
    A(1, 1) = 1 + params.dt * s_k(3) * s_k(2);
    A(1, 2) = params.dt * s_k(3) * s_k(1);
    A(1, 3) = params.dt * s_k(1) * s_k(2);

    A(2, 1) = -params.dt * s_k(2) * s_k(3);
    A(2, 2) = 1 - params.dt * (s_k(1) * s_k(3) - params.beta);
    A(2, 3) = -params.dt * s_k(1) * s_k(2);

    A(3, 3) = 1 + params.dt * params.gamma;

    B = eye(3);
end

% Observation equation Jacobian
function [C, D] = obs_jacobian_flipped(u, s_k, v_bar, params, k)
    % OBS_JACOBIAN_FLIPPED Evaluate analytic observation and measurement-noise Jacobians.

    if isequal(params.obs_type, 'NEWCASES')
        C = [s_k(2) * s_k(3), s_k(1) * s_k(3), s_k(1) * s_k(2)];
        D = 1;
    elseif isequal(params.obs_type, 'TOTALCASES')
        C = [-1, 0, 0]; % following the revised model that takes total cases as input
        D = 1;
    else
        error('unknown observation type');
    end
end

% State equation Hessian terms
function [fs, Cs, fw, Cw] = state_hessian_terms_flipped(u, s_k, Pk, w_bar, Qk, params, k)
    % STATE_HESSIAN_TERMS_FLIPPED Return the model-specific Gaussian second-order corrections.
    fs = zeros(3, 1);
    Cs = zeros(3);

    fw = zeros(3, 1);
    Cw = zeros(3);

end

% Observation equation Hessian terms
function [gs, Gsp, gv, Gvp] = obs_hessian_terms_flipped(u, s_k, Pk, v_bar, Rk, params, k)
    % OBS_HESSIAN_TERMS_FLIPPED Return observation mean and covariance second-order corrections.
    gs = zeros(1);
    Gsp = zeros(1);

    gv = zeros(1);
    Gvp = zeros(1);

end
