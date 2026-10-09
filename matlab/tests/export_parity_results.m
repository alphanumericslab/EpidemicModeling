function export_parity_results(output_file)
% EXPORT_PARITY_RESULTS Save actual MATLAB results for cross-language comparison.
% export_parity_results(output_file); default: reports/matlab_results.mat.
% Compare with python scripts/compare_matlab_python.py reports/matlab_results.mat.
% Requires base MATLAB; nonlinear fitting is included only when nlinfit exists.
% Author: Reza Sameni | Emory University
root = fileparts(fileparts(fileparts(mfilename('fullpath')))); addpath(fullfile(root,'matlab','functions'));
if nargin < 1, output_file = fullfile(root,'reports','matlab_results.mat'); end
out = struct();
[a,b,c,d,e] = seirp(.4,.3,.12,.04,.09,.002,0,.998,.001,.001,0,0,10,.1); out.seirp = [a;b;c;d;e];
[a,b,c,d,e] = seirp_saturated_resource(.4,.3,.12,.04,0,.998,.001,.001,0,0,10,.1,.09,.03,.002,.015,.01,.04); out.saturated = [a;b;c;d;e];
[a,b] = si_controlled(.35,.1,.999,.001,100,.1); out.si = [a;b];
u = [linspace(0,3,40);linspace(3,0,40)]; noise = reshape(sin(0:119),3,40);
[a,b,c] = si_alpha_controlled(u,.999,.001,.3,[3;3],0,5,1/7,[.06;.04],.02,.1,.001,.001,.01,40,.1,noise); out.si_alpha = [a;b;c];
cases = 25*exp(.035*(0:39));
[a,b,c,d] = rt_exp_fit_gen_ratios(cases,7,5,1); out.ratios = [a;b;c;d];
[a,b,c,d] = rt_exp_fit_log_lin_reg(cases,7,1); out.log_fit = [a;b;c;d];
[a,b,c,d] = rt_exp_fit_log_lin_reg(cases,8,1,false); out.log_fit_centered = [a;b;c;d];
if exist('nlinfit','file')==2
    [a,b,c,d] = rt_exp_fit_nonlin_ls(cases,7,1); out.nonlinear_fit = [a;b;c;d];
    [a,b,c,d] = rt_exp_fit_nonlin_ls(cases,8,1,false); out.nonlinear_fit_centered = [a;b;c;d];
else
    out.optional_nonlinear_missing = true;
end
x = cases+2*sin(0:39); x(16:18) = nan;
exp_names = {'s_minus','s_plus','p_minus','p_plus','gain','s_smooth','p_smooth','innovations','rho'};
for order = 1:2
    values = cell(1,9);
    [values{:}] = rt_exp_fit_ekf(x,[25;.03],[1 1 .2],[0;0],0,diag([4 .001]),diag([1 1e-5]),4,1,1,7,order);
    for j = 1:9, out.(sprintf('exp_%d_%s',order,exp_names{j})) = values{j}; end
end
p = default_si_params([3;3]); p.a = [.06;.04]; p.b = .02; y = .0003*ones(1,40); y(16:18) = nan;
common = {p,[.999;.001;.3],diag([1e-6 1e-6 .01]),nan(3,1),nan(3),zeros(3,1),0,diag([1e-8 1e-8 1e-5]),1e-8,1,1,7,1};
names = {'u_opt','u_opt_smooth','s_minus','s_plus','s_smooth','p_minus','p_plus','p_smooth','k_gain','innovations','rho'}; values = cell(1,11);
[values{:}] = si_alpha_model_ekf(u,y,common{:}); out = add_filter_fields(out,'si_filter_',names,values);
q = repmat(diag([1e-8 1e-8 1e-5]),1,1,40).*reshape(linspace(1,2,40),1,1,[]); variance = linspace(1e-8,2e-8,40);
variable = common; variable(8:9) = {q,variance};
[values{:}] = si_alpha_model_ekf(u,y,variable{:}); out = add_filter_fields(out,'si_variable_',names,values);
back = common; back(4:5) = {[.995;.002;.3],diag([1e-6 1e-6 .01])};
[values{:}] = si_alpha_model_backward_ekf(u,y,back{:}); out = add_filter_fields(out,'si_backward_',names,values);
control_common = {p,[.999;.001;.3;0;0;0],diag([1e-6 1e-6 .01 1 1 1]),[nan(3,1);zeros(3,1)],nan(6),zeros(6,1),0,diag([1e-8 1e-8 1e-5 1e-4 1e-4 1e-4]),1e-8,1,1,7,1};
[values{:}] = si_alpha_model_ekf_opt_controlled(nan(size(u)),nan(1,40),control_common{:}); out = add_filter_fields(out,'si_control_',names,values);
back = control_common; back(4:5) = {[.995;.002;.3;0;0;0],diag([1e-6 1e-6 .01 1 1 1])};
[values{:}] = si_alpha_model_backward_ekf_opt_controlled(u,y,back{:}); out = add_filter_fields(out,'si_backward_control_',names,values);
legacy_names = {'u_opt','s_minus','s_plus','s_smooth','p_minus','p_plus','p_smooth','k_gain','innovations','rho'};
legacy_values = cell(1,10); [legacy_values{:}] = new_case_ekf_estimator_with_optimal_npi(u,y,control_common{:});
out = add_filter_fields(out,'legacy_',legacy_names,legacy_values);
[a,b] = npi_cost(cases,u,[1;1.5]); out.costs = [a b];
xx = linspace(-2,2,21); out.exp_layer = exp(.7*xx); out.tanh_layer = .7*tanh(xx/.7);
initial = zeros(11); initial(6,6) = 1; out.diffusion = diffusion_2d(initial,1,.2,1,10);
out.motion = population_motion_2d([.2 .3;.7 .8],[.14 .09;-.11 .07],.1,30);
design = [ones(40,1) (0:39)'/40 sin((0:39)').^2]; out.nnls = nonnegative_least_squares(design,design*[.1;.2;.3]);
day = 0:49; train_u = [2*(day>=25);2*(day>=35)];
truth.population = 1e6; truth.params = p; truth.state = [.9999;.0001;.3]; truth.covariance = diag([1e-8 1e-8 .01]);
[train_cases,~] = forecast_npi(truth,train_u); model = fit_npi_model(cumsum(train_cases),train_u,1e6,[3;3]);
out.trained_coef = model.coefficients; out.trained_refined = model.refined_coefficients; out.trained_state = model.state; out.trained_cov = model.covariance;
[out.forecast,out.forecast_states] = forecast_npi(model,repmat(train_u(:,end),1,10));
[out.optimal_controls,out.optimal_cases,out.optimal_states] = optimal_npi(model,10,[1;1.5],.3);
coder_folder = fullfile(root,'matlab','codegen'); addpath(coder_folder,'-begin'); cleanup = onCleanup(@() rmpath(coder_folder));
state = [.8;.1;.3;.05;.02;.03]; control = [1;2];
out.coder_state_margin = state_hard_margins(state,p); out.coder_obs_margin = obs_hard_margins(-.1,p);
[out.coder_u,out.coder_state] = nlin_state_update(control,state,zeros(6,1),p);
out.coder_obs = nlin_obs_update(control,state,0,p);
[out.coder_a,out.coder_b] = state_jacobians(control,state,zeros(6,1),p);
[out.coder_c,out.coder_d] = obs_jacobian(control,state,0,p);
[out.coder_fs,out.coder_cs,out.coder_fw,out.coder_cw] = state_hessian_terms(control,state,eye(6),zeros(6,1),eye(6),p);
[out.coder_gs,out.coder_gsp,out.coder_gv,out.coder_gvp] = obs_hessian_terms(control,state,eye(6),0,1,p);
save(output_file,'-struct','out','-v7'); fprintf('Saved %d MATLAB arrays to %s\n',numel(fieldnames(out)),output_file);
end
function out = add_filter_fields(out,prefix,names,values)
% ADD_FILTER_FIELDS Add named filter result arrays to the export struct.
for j = 1:numel(names), out.([prefix names{j}]) = values{j}; end
end
