function params = default_si_params(u_max, epsilon, weights)
% DEFAULT_SI_PARAMS Reproducible daily SI-alpha parameters for example exercises.
% params = default_si_params(u_max,epsilon,weights)
% u_max: intervention bounds; epsilon: human/NPI cost mix (default .5);
% weights: nonnegative intervention costs (default ones). Returns a struct.
% beta=-log(.01)/21 and gamma=1/7 are historical educational assumptions.
% Author: Reza Sameni | Emory University
if nargin < 2, epsilon = .5; end
if nargin < 3, weights = ones(numel(u_max),1); end
params.dt = 1; params.beta = -log(.01)/21; params.gamma = 1/7;
params.a = zeros(numel(u_max),1); params.b = 0;
params.u_min = zeros(numel(u_max),1); params.u_max = u_max(:);
params.alpha_min = 0; params.alpha_max = 5;
params.epsilon = epsilon; params.w = weights(:); params.sigma = 10000;
params.obs_type = 'NEWCASES'; params.s_min = 0; params.i_min = 0;
end
