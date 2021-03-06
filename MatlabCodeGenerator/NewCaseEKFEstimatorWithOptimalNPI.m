function [u_opt, S_MINUS, S_PLUS, P_MINUS, P_PLUS, K_GAIN, S_SMOOTH, P_SMOOTH, innovations, rho] = NewCaseEKFEstimatorWithOptimalNPI(u, x, params, s_init, Ps_init, s_final, Ps_final, w_bar, v_bar, Q_w, R_v, beta, gamma, inv_monitor_len, order)
% Estimates the parameters of an exponential fit using an Extended Kalman
% Filter (EKF) and Extended Kalman Smoother (EKS) over the number of new cases
%
% Reza Sameni
% Jan 2021
% Email: reza.sameni@gmail.com

T = size(x, 2); % number of time samples
n = size(x, 1); % number of observations
m = length(s_init); % number of state variables
% l = length(w_bar); % number of process noises
% p = length(v_bar); % number of observation noises

%//////////////////////////////////////////////////////////////////////////
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
for k = 1 : T
    % Store results of s_minus from previous iteration
    S_MINUS(:, k) = sk_minus;
    P_MINUS(:, :, k) = Pk_minus;
    
    if(order == 1)
        gs = zeros(n, 1);
        Gsp = zeros(n);
        gv = zeros(n, 1);
        Gvp = zeros(n);
    elseif(order == 2)
        [gs, Gsp, gv, Gvp] = ObsHessianTerms(u(:, k), sk_minus, Pk_minus, v_bar, R, params);
    else
        error('Undefined order');
    end
    
    % Calculate s(k|k) and P(k|k)
    [Ck_minus, Dk_minus] = ObsJacobian(u(:, k), sk_minus, v_bar, params);
    xk_minus = NlinObsUpdate(u(:, k), sk_minus, v_bar, params) + gs + gv;

    % Apply hard margins on observations
    xk_minus = ObsHardMargins(xk_minus, params);
    
    % time update if observation is valid
    if(~isnan(x(:, k)))
        innovations(:, k) = x(:, k) - xk_minus;
        Kgain = Pk_minus * Ck_minus' / (Ck_minus * Pk_minus * Ck_minus' + gamma * (Dk_minus * R * Dk_minus') + Gsp + Gvp); % Kalman gain
        Pk_plus = (eye(m) - Kgain * Ck_minus) * Pk_minus / gamma;
        sk_plus = sk_minus + Kgain * innovations(:, k);                 % As posteriori state estimate
    else
        innovations(:, k) = 0;
        Kgain = zeros(m, n); % Kalman gain
        Pk_plus = Pk_minus;
        sk_plus = sk_minus; % As posteriori state estimate
    end
    
    % Apply hard margins on states
    sk_plus = StateHardMargins(sk_plus, params);
    
    if(order == 1)
        fs = zeros(m, 1);
        Fsp = zeros(m);
        fw = zeros(m, 1);
        Fwp = zeros(m);
    elseif(order == 2)
        [fs, Fsp, fw, Fwp] = StateHessianTerms(u(:, k), sk_plus, Pk_plus, w_bar, Q, params);
    else
        error('Undefined order');
    end
    
    % Calculate s(k+1|k) and P(k+1|k) for k+1
    [u_opt(:, k), sk_minus] = NlinStateUpdate(u(:, k), sk_plus, w_bar, params); % State update
    sk_minus = sk_minus + fs + fw; % Add second order terms (is available)
    [Ak_plus, Bk_plus] = StateJacobians(u(:, k), sk_plus, w_bar, params);
    Pk_minus = (Ak_plus * Pk_plus * Ak_plus') + (Bk_plus * Q * Bk_plus') + Fsp + Fwp; % Cov. matrix update
    
    % Apply hard margins on states
    sk_minus = StateHardMargins(sk_minus, params);

    % Store results of s_plus
    S_PLUS(:, k) = sk_plus;
    P_PLUS(:, :, k) = Pk_plus;
    K_GAIN(:, :, k) = Kgain;
    
    % Monitoring the innovation variance and update the observation noise
    stats_counter = min(k, inv_monitor_len);
    InnovationsMean = cat(2, innovations(:, k), InnovationsMean(:, 1 : inv_monitor_len - 1));
    mu_k = sum(InnovationsMean, 2)/stats_counter;
    %     cc = innovations(:, k)*innovations(:, k)'; % without mean cancellation
    cc = (innovations(:, k) - mu_k) * (innovations(:, k) - mu_k)'; % with mean cancellation
    InnovationsCov = cat(3, cc, InnovationsCov(:, :, 1 : inv_monitor_len - 1));
    InnovationsCovNormalized = cat(3, cc / R, InnovationsCovNormalized(:, :, 1 : inv_monitor_len - 1));
    rho(:, :, k) = sum(InnovationsCovNormalized, 3) / stats_counter;
    if(~isequal(beta, 1) && ~isnan(x(:, k)))
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
[row,col] = ind2sub(size(Ps_final), fixed_end_covs);
P_SMOOTH(row, col, T) = Ps_final(row,col);

for k = T - 1 : -1 : 1
    sk_plus = S_PLUS(:, k);
    Ak_plus = StateJacobians(u(:, k), sk_plus, w_bar, params);
    J = (P_PLUS(:, :, k) * Ak_plus') / (P_MINUS(:, :, k + 1));
    S_SMOOTH(:, k) = S_PLUS(:, k) + J * (S_SMOOTH(:, k + 1) - S_MINUS(:, k + 1));
    
    % Apply hard margins on states
    S_SMOOTH(:, k) = StateHardMargins(S_SMOOTH(:, k), params);
    
    P_SMOOTH(:, :, k) = P_PLUS(:, :, k) - J * (P_MINUS(:, :, k+1) - P_SMOOTH(:, :, k+1)) * J';
end

% Squeeze excess dimensions if applicable
rho = squeeze(rho);
end
