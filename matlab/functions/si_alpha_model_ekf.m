function [u_opt, u_opt_smooth, S_MINUS, S_PLUS, S_SMOOTH, P_MINUS, P_PLUS, P_SMOOTH, K_GAIN, innovations, rho] = si_alpha_model_ekf(u, x, params, s_init, Ps_init, s_final, Ps_final, w_bar, v_bar, Q_w, R_v, beta, gamma, inv_monitor_len, order)
% SI_ALPHA_MODEL_EKF Filter and smooth three SI-alpha states from new or total case observations.
%
% Syntax: [u_opt, u_opt_smooth, S_MINUS, S_PLUS, S_SMOOTH, P_MINUS, P_PLUS, P_SMOOTH, K_GAIN, innovations, rho] = si_alpha_model_ekf(u, x, params, s_init, Ps_init, s_final, Ps_final, w_bar, v_bar, Q_w, R_v, beta, gamma, inv_monitor_len, order)
% Author: Reza Sameni | Emory University
% Reference: Sameni (2020), arXiv:2003.11371; Sameni (2022),
% doi:10.1109/JSTSP.2021.3129118.
%
% Inputs, outputs, units, shapes and endpoint conventions: docs/api.md.
% All rates use a consistent inverse time unit; dt uses that time unit.
% u/x: controls/observations-by-time; state means: states-by-time;
% covariances/gains: dimensions-by-dimensions-by-time. Missing observations
% skip correction. NaN final mean/covariance entries remain unconstrained.
% SI-alpha params fields and callback signatures are listed in docs/api.md.
% beta in [0,1] controls R adaptation; gamma in (0,1] stabilizes covariance;
% inv_monitor_len is a positive integer; order is 1 or 2.

handles.state_hard_margins = @state_hard_margins;
handles.obs_hard_margins = @obs_hard_margins;
handles.nlin_state_update = @nlin_state_update;
handles.nlin_obs_update = @nlin_obs_update;
handles.state_jacobians = @state_jacobians;
handles.obs_jacobian = @obs_jacobian;
handles.state_hessian_terms = @state_hessian_terms;
handles.obs_hessian_terms = @obs_hessian_terms;

% Call the generic EKF function with the desired function handles
[u_opt, u_opt_smooth, S_MINUS, S_PLUS, S_SMOOTH, P_MINUS, P_PLUS, P_SMOOTH, K_GAIN, innovations, rho] = generic_extended_kalman_filter(u, x, handles, params, s_init, Ps_init, s_final, Ps_final, w_bar, v_bar, Q_w, R_v, beta, gamma, inv_monitor_len, order);

end
% THE SYSTEM EQUATIONS

% Hard margins on state vectors
function s_k = state_hard_margins(s_k, params, k)
% STATE_HARD_MARGINS Constrain population fractions and contact rate; preserve costates.
s_k(1) = min(1.0, max(params.s_min, s_k(1)));
s_k(2) = min(1.0, max(params.i_min, s_k(2)));
s_k(3) = min(params.alpha_max, max(params.alpha_min, s_k(3)));
end

% Hard margins on observations
function x_k = obs_hard_margins(x_k, params, k)
% OBS_HARD_MARGINS Constrain predicted observations to nonnegative values.
    x_k = max(0, x_k);
end

% Nonlinear state update
function [u, s_k_plus_one] = nlin_state_update(u, s_k, w_bar, params, k)
% NLIN_STATE_UPDATE Evaluate one Euler transition and fill missing controls when applicable.

s_k_plus_one = zeros(3, 1);

% State equations
s_k_plus_one(1) = max(params.s_min, min(1.0, s_k(1) - params.dt * s_k(3) * s_k(1) * s_k(2)));
s_k_plus_one(2) = max(params.i_min, min(1.0, s_k(2) + params.dt * (s_k(3) * s_k(1) * s_k(2) - params.beta * s_k(2))));
s_k_plus_one(3) = max(params.alpha_min, min(params.alpha_max, s_k(3) + params.dt * (-params.gamma * s_k(3) + params.gamma * params.b + params.gamma * params.a'*(params.u_max - u))));

end

% Nonlinear observation update
function x_k = nlin_obs_update(u, s_k, v_bar, params, k)
% NLIN_OBS_UPDATE Evaluate the configured new-case or total-case observation.
    if(isequal(params.obs_type, 'NEWCASES'))
        x_k = s_k(1) * s_k(2) * s_k(3) + v_bar;
    elseif(isequal(params.obs_type, 'TOTALCASES'))
        x_k = 1 - s_k(1) + v_bar; % following the revised model that takes total cases as input
    else
        error('unknown observation type');
    end
end

% State equation Jacobian
function [A, B] = state_jacobians(u, s_k, w_bar, params, k)
% STATE_JACOBIANS Evaluate analytic state-transition and process-noise Jacobians.

A = zeros(3);
A(1, 1) = 1 - params.dt * s_k(3) * s_k(2);
A(1, 2) = - params.dt * s_k(3) * s_k(1);
A(1, 3) = - params.dt * s_k(1) * s_k(2);

A(2, 1) = params.dt * s_k(2) * s_k(3);
A(2, 2) = 1 + params.dt * (s_k(1) * s_k(3) - params.beta);
A(2, 3) = params.dt * s_k(1) * s_k(2);

A(3, 3) = 1 - params.dt * params.gamma;

B = eye(3);
end

% Observation equation Jacobian
function [C, D] = obs_jacobian(u, s_k, v_bar, params, k)
% OBS_JACOBIAN Evaluate analytic observation and measurement-noise Jacobians.
    if(isequal(params.obs_type, 'NEWCASES'))
        C = [s_k(2)*s_k(3), s_k(1)*s_k(3), s_k(1)*s_k(2)];
        D = 1;
    elseif(isequal(params.obs_type, 'TOTALCASES'))
        C = [-1, 0, 0]; % following the revised model that takes total cases as input
        D = 1;
    else
        error('unknown observation type');
    end
end

% State equation Hessian terms
function [fs, Cs, fw, Cw] = state_hessian_terms(u, s_k, Pk, w_bar, Qk, params, k)
% STATE_HESSIAN_TERMS Return the model-specific Gaussian second-order corrections.
fs = zeros(3, 1);
Cs = zeros(3);

fw = zeros(3, 1);
Cw = zeros(3);

end

% Observation equation Hessian terms
function [gs, Gsp, gv, Gvp] = obs_hessian_terms(u, s_k, Pk, v_bar, Rk, params, k)
% OBS_HESSIAN_TERMS Return observation mean and covariance second-order corrections.
gs = zeros(1);
Gsp = zeros(1);

gv = zeros(1);
Gvp = zeros(1);

end
