function [controls, cases, states] = optimal_npi(model, horizon, weights, epsilon, iterations)
% OPTIMAL_NPI Solve bounded SI-alpha control by repeated EKF/EKS shooting.
% [controls,cases,states] = optimal_npi(model,horizon,weights,epsilon,iterations)
% model: fitted model; horizon: integer days; weights: nonnegative NPI vector;
% epsilon: cost mix [0,1] (default .5); iterations: shooting passes (default 8).
% Controls are interventions-by-days; cases are daily counts; states are 3-by-days.
% Zero terminal costates are imposed. This finite-iteration solver does not guarantee
% global optimality. Reference: Sameni (2022), doi:10.1109/JSTSP.2021.3129118.
% Author: Reza Sameni | Emory University
if nargin < 4, epsilon = .5; end
if nargin < 5, iterations = 8; end
assert(epsilon >= 0 && epsilon <= 1 && horizon >= 1 && horizon == fix(horizon) && iterations >= 1, 'Invalid control settings.');
p = model.params; p.epsilon = epsilon; p.w = weights(:); p.sigma = 10000;
p.a = p.a(:); p.u_min = p.u_min(:); p.u_max = p.u_max(:);
assert(numel(p.w) == numel(p.a) && all(p.w >= 0), 'Invalid intervention costs.');
initial = [model.state(:); zeros(3,1)]; cov = zeros(6); cov(1:3,1:3) = model.covariance; cov(4:6,4:6) = eye(3);
terminal = [nan(3,1); zeros(3,1)]; terminal_cov = nan(6); terminal_cov(4:6,4:6) = 0;
u = nan(numel(p.a),horizon);
for k = 1:iterations
    [u_opt,u_smooth,~,~,smoothed] = si_alpha_model_ekf_opt_controlled(u,nan(1,horizon),p,initial,cov,terminal,terminal_cov,zeros(6,1),0,diag([1e-10 1e-10 1e-6 1e-4 1e-4 1e-4]),1e-6,1,1,21,1);
    initial(4:6) = smoothed(4:6,1);
end
controls = u_smooth; controls(:,end) = u_opt(:,end);
[cases,states] = forecast_npi(model,controls);
end
