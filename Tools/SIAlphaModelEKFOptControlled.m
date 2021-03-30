function [u_opt, u_opt_smooth, S_MINUS, S_PLUS, S_SMOOTH, P_MINUS, P_PLUS, P_SMOOTH, K_GAIN, innovations, rho] = SIAlphaModelEKFOptControlled(u, x, params, s_init, Ps_init, s_final, Ps_final, w_bar, v_bar, Q_w, R_v, beta, gamma, inv_monitor_len, order)
% Estimates the parameters of an exponential fit using an Extended Kalman
% Filter (EKF) and Extended Kalman Smoother (EKS) over the number of new cases
%
% Reza Sameni
% Jan 2021
% Email: reza.sameni@gmail.com

handles.StateHardMargins = @StateHardMargins;
handles.ObsHardMargins = @ObsHardMargins;
handles.NlinStateUpdate = @NlinStateUpdate;
handles.NlinObsUpdate = @NlinObsUpdate;
handles.StateJacobians = @StateJacobians;
handles.ObsJacobian = @ObsJacobian;
handles.StateHessianTerms = @StateHessianTerms;
handles.ObsHessianTerms = @ObsHessianTerms;

% Call the generic EKF function with the desired function handles
[u_opt, u_opt_smooth, S_MINUS, S_PLUS, S_SMOOTH, P_MINUS, P_PLUS, P_SMOOTH, K_GAIN, innovations, rho] = GenericExtendedKalmanFilter(u, x, handles, params, s_init, Ps_init, s_final, Ps_final, w_bar, v_bar, Q_w, R_v, beta, gamma, inv_monitor_len, order);

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% THE SYSTEM EQUATIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Hard margins on state vectors
function s_k = StateHardMargins(s_k, params, k)
s_k(1) = min(1.0, max(0, s_k(1)));
s_k(2) = min(1.0, max(0, s_k(2)));
s_k(3) = min(params.alpha_max, max(params.alpha_min, s_k(3)));
end

% Hard margins on observations
function x_k = ObsHardMargins(x_k, params, k)
    x_k = max(0, x_k);
end

% Nonlinear state update
function [u, s_k_plus_one] = NlinStateUpdate(u, s_k, w_bar, params, k)

% if(isnan(u)) % Optimal control
%     u = zeros(length(params.u_max), 1);
%     phi = params.epsilon * params.w - params.gamma * s_k(6) * params.a;
%     u(phi >= 0) = params.u_min(phi >= 0);
%     u(phi < 0) = params.u_max(phi < 0);
% end

% Optimal control replaces nans (by design feature)
phi = params.epsilon * params.w - params.gamma * s_k(6) * params.a;
for kk = 1 : length(u)
    if(isnan(u(kk)))
        if(phi(kk) > 0) % or >= 0 ?
            u(kk) = params.u_min(kk);
        else
            u(kk) = params.u_max(kk);
        end
    end
end

rho = s_k(4) - s_k(5) - (1 - params.epsilon);

s_k_plus_one = zeros(6, 1);

% State equations
s_k_plus_one(1) = max(0.0, min(1.0, s_k(1) - params.dt * s_k(3) * s_k(1) * s_k(2)));
s_k_plus_one(2) = max(0.0, min(1.0, s_k(2) + params.dt * (s_k(3) * s_k(1) * s_k(2) - params.beta * s_k(2))));
s_k_plus_one(3) = max(params.alpha_min, min(params.alpha_max, s_k(3) + params.dt * (-params.gamma * s_k(3) + params.gamma * params.b + params.gamma * params.a'*(params.u_max - u))));

% Costate equations
s_k_plus_one(4) = s_k(4) + params.dt * rho * s_k(3) * s_k(2);
s_k_plus_one(5) = s_k(5) + params.dt * (rho * s_k(3) * s_k(1) + params.beta * s_k(5));
s_k_plus_one(6) = s_k(6) + params.dt * (rho * s_k(1) * s_k(2) + params.gamma * s_k(6));

end

% Nonlinear observation update
function x_k = NlinObsUpdate(u, s_k, v_bar, params, k)
    if(isequal(params.obs_type, 'NEWCASES'))
        x_k = s_k(1) * s_k(2) * s_k(3) + v_bar;
    elseif(isequal(params.obs_type, 'TOTALCASES'))
        x_k = 1 - s_k(1) + v_bar; % following the revised model that takes total cases as input
    else
        error('unknown observation type');
    end
end

% State equation Jacobian
function [A, B] = StateJacobians(u, s_k, w_bar, params, k)

A = zeros(6);
A(1, 1) = 1 - params.dt * s_k(3) * s_k(2);
A(1, 2) = - params.dt * s_k(3) * s_k(1);
A(1, 3) = - params.dt * s_k(1) * s_k(2);

A(2, 1) = params.dt * s_k(2) * s_k(3);
A(2, 2) = 1 + params.dt * (s_k(1) * s_k(3) - params.beta);
A(2, 3) = params.dt * s_k(1) * s_k(2);

A(3, 3) = 1 - params.dt * params.gamma;
% Sigmoid function:
% % % if(isnan(u)) % entry non-zero only for optimal control (that is costate dependent)
% % %     x = -params.sigma * (s_k(6) - (params.epsilon * params.w)./(params.gamma .* params.a) );
% % %     A(3, 6) = - params.gamma * params.dt * params.sigma * params.a' * ((params.u_max - params.u_min) .* exp(x)./(1 + exp(x)).^2);
% % % end

% Linear slope:
phi = params.epsilon * params.w - params.gamma * s_k(6) * params.a;
for kk = 1 : length(u)
    if(isnan(u(kk)))
        if(phi(kk) > -1.0/params.sigma && phi(kk) < 1.0/params.sigma)
            A(3, 6) = A(3, 6) - params.gamma * params.dt * (params.sigma/2) * params.a(kk) * (params.u_max(kk) - params.u_min(kk));
        end
    end
end


rho = s_k(4) - s_k(5) - (1 - params.epsilon);
A(4, 2) = params.dt * s_k(3) * rho;
A(4, 3) = params.dt * s_k(2) * rho;
A(4, 4) = 1 + params.dt * s_k(2) * s_k(3);
A(4, 5) = - params.dt * s_k(2) * s_k(3);

A(5, 1) = params.dt * s_k(3) * rho;
A(5, 3) = params.dt * s_k(1) * rho;
A(5, 4) = params.dt * s_k(1) * s_k(3);
A(5, 5) = 1 - params.dt * (s_k(1) * s_k(3) - params.beta);

A(6, 1) = params.dt * s_k(2) * rho;
A(6, 2) = params.dt * s_k(1) * rho;
A(6, 4) = params.dt * s_k(1) * s_k(2);
A(6, 5) = -params.dt * s_k(1) * s_k(2);
A(6, 6) = 1 + params.dt * params.gamma;

B = eye(6);
end

% Observation equation Jacobian
function [C, D] = ObsJacobian(u, s_k, v_bar, params, k)
    if(isequal(params.obs_type, 'NEWCASES'))
        C = [s_k(2)*s_k(3), s_k(1)*s_k(3), s_k(1)*s_k(2), 0 , 0, 0];
        D = 1;
    elseif(isequal(params.obs_type, 'TOTALCASES'))
        C = [-1, 0, 0, 0 , 0, 0]; % following the revised model that takes total cases as input
        D = 1;
    else
        error('unknown observation type');
    end
end

% State equation Hessian terms
function [fs, Cs, fw, Cw] = StateHessianTerms(u, s_k, Pk, w_bar, Qk, params, k)
fs = zeros(6, 1);
Cs = zeros(6);

fw = zeros(6, 1);
Cw = zeros(6);

end

% Observation equation Hessian terms
function [gs, Gsp, gv, Gvp] = ObsHessianTerms(u, s_k, Pk, v_bar, Rk, params, k)
gs = zeros(1);
Gsp = zeros(1);

gv = zeros(1);
Gvp = zeros(1);

end
