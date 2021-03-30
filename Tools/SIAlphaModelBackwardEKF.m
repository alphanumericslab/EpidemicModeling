function [u_opt, u_opt_smooth, S_MINUS, S_PLUS, S_SMOOTH, P_MINUS, P_PLUS, P_SMOOTH, K_GAIN, innovations, rho] = SIAlphaModelBackwardEKF(u, x, params, s_init, Ps_init, s_final, Ps_final, w_bar, v_bar, Q_w, R_v, beta, gamma, inv_monitor_len, order)
% Estimates the parameters of an exponential fit using an Extended Kalman
% Filter (EKF) and Extended Kalman Smoother (EKS) over the number of new cases
%
% Reza Sameni
% Jan 2021
% Email: reza.sameni@gmail.com

handles.StateHardMargins = @StateHardMarginsFlipped;
handles.ObsHardMargins = @ObsHardMarginsFlipped;
handles.NlinStateUpdate = @NlinStateUpdateFlipped;
handles.NlinObsUpdate = @NlinObsUpdateFlipped;
handles.StateJacobians = @StateJacobiansFlipped;
handles.ObsJacobian = @ObsJacobianFlipped;
handles.StateHessianTerms = @StateHessianTermsFlipped;
handles.ObsHessianTerms = @ObsHessianTermsFlipped;

% time flip the samples
u_flipped = u(:, end : -1 : 1);
x_flipped = x(end : -1 : 1);
s_init_flipped = s_final;
s_final_flipped = s_init;
Ps_init_flipped = Ps_final;
Ps_final_flipped = Ps_init;

% Call the generic EKF function with the desired function handles
[u_opt_flipped, u_opt_smooth_flipped, S_MINUS_flipped, S_PLUS_flipped, S_SMOOTH_flipped, P_MINUS_flipped, P_PLUS_flipped, P_SMOOTH_flipped, K_GAIN_flipped, innovations_flipped, rho_flipped] = GenericExtendedKalmanFilter(u_flipped, x_flipped, handles, params, s_init_flipped, Ps_init_flipped, s_final_flipped, Ps_final_flipped, w_bar, v_bar, Q_w, R_v, beta, gamma, inv_monitor_len, order);

% flip back the results in time
u_opt = u_opt_flipped(:, end : -1 : 1);
u_opt_smooth = u_opt_smooth_flipped(:, end : -1 : 1);
S_MINUS = S_MINUS_flipped(:, end : -1 : 1);
S_PLUS = S_PLUS_flipped(:, end : -1 : 1);
S_SMOOTH = S_SMOOTH_flipped(:, end : -1 : 1);
P_MINUS = P_MINUS_flipped(:, :, end : -1 : 1);
P_PLUS = P_PLUS_flipped(:, :, end : -1 : 1);
P_SMOOTH = P_SMOOTH_flipped(:, :, end : -1 : 1);
K_GAIN = K_GAIN_flipped(:, :, end : -1 : 1);
innovations = innovations_flipped(:, end : -1 : 1);
rho = rho_flipped(:, :, end : -1 : 1);

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% THE SYSTEM EQUATIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Hard margins on state vectors
function s_k = StateHardMarginsFlipped(s_k, params, k)
s_k(1) = min(1.0, max(0, s_k(1)));
s_k(2) = min(1.0, max(0, s_k(2)));
s_k(3) = min(params.alpha_max, max(params.alpha_min, s_k(3)));
end

% Hard margins on observations
function x_k = ObsHardMarginsFlipped(x_k, params, k)
    x_k = max(0, x_k);
end

% Nonlinear state update
function [u, s_k_plus_one] = NlinStateUpdateFlipped(u, s_k, w_bar, params, k)

s_k_plus_one = zeros(3, 1);

% State equations
s_k_plus_one(1) = max(0.0, min(1.0, s_k(1) + params.dt * s_k(3) * s_k(1) * s_k(2)));
s_k_plus_one(2) = max(0.0, min(1.0, s_k(2) - params.dt * (s_k(3) * s_k(1) * s_k(2) - params.beta * s_k(2))));
s_k_plus_one(3) = max(params.alpha_min, min(params.alpha_max, s_k(3) - params.dt * (-params.gamma * s_k(3) + params.gamma * params.b + params.gamma * params.a'*(params.u_max - u))));

end

% Nonlinear observation update
function x_k = NlinObsUpdateFlipped(u, s_k, v_bar, params, k)
    if(isequal(params.obs_type, 'NEWCASES'))
        x_k = s_k(1) * s_k(2) * s_k(3) + v_bar;
    elseif(isequal(params.obs_type, 'TOTALCASES'))
        x_k = 1 - s_k(1) + v_bar; % following the revised model that takes total cases as input
    else
        error('unknown observation type');
    end
end

% State equation Jacobian
function [A, B] = StateJacobiansFlipped(u, s_k, w_bar, params, k)

A = zeros(3);
A(1, 1) = 1 + params.dt * s_k(3) * s_k(2);
A(1, 2) = params.dt * s_k(3) * s_k(1);
A(1, 3) = params.dt * s_k(1) * s_k(2);

A(2, 1) = - params.dt * s_k(2) * s_k(3);
A(2, 2) = 1 - params.dt * (s_k(1) * s_k(3) - params.beta);
A(2, 3) = - params.dt * s_k(1) * s_k(2);

A(3, 3) = 1 + params.dt * params.gamma;

B = eye(3);
end

% Observation equation Jacobian
function [C, D] = ObsJacobianFlipped(u, s_k, v_bar, params, k)
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
function [fs, Cs, fw, Cw] = StateHessianTermsFlipped(u, s_k, Pk, w_bar, Qk, params, k)
fs = zeros(3, 1);
Cs = zeros(3);

fw = zeros(3, 1);
Cw = zeros(3);

end

% Observation equation Hessian terms
function [gs, Gsp, gv, Gvp] = ObsHessianTermsFlipped(u, s_k, Pk, v_bar, Rk, params, k)
gs = zeros(1);
Gsp = zeros(1);

gv = zeros(1);
Gvp = zeros(1);

end
