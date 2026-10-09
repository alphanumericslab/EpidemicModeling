function tests = test_npi_spatial_layers()
% TEST_NPI_SPATIAL_LAYERS NNLS, training/control bounds, spatial invariants and layers.
setup_once([]);
tests = functiontests(localfunctions);
end
function setup_once(test_case)
% SETUP_ONCE Add function library.
addpath(fullfile(fileparts(fileparts(mfilename('fullpath'))),'functions'));
end
function test_nnls_solution(test_case)
% TEST_NNLS_SOLUTION Recover a nonnegative known regression coefficient vector.
x = [ones(40,1) (0:39)'/40 sin((0:39)').^2]; expected = [.1;.2;.3];
coef = nonnegative_least_squares(x,x*expected); verifyEqual(test_case,coef,expected,'AbsTol',1e-7);
end
function test_pipeline_bounds(test_case)
% TEST_PIPELINE_BOUNDS Fitted model forecasts remain valid and controls obey bounds.
day = 0:49; u = [2*(day>=25);2*(day>=35)];
truth.population = 1e6; truth.params = default_si_params([3;3]); truth.params.a = [.06;.04]; truth.params.b = .02;
truth.state = [.9999;.0001;.3]; truth.covariance = diag([1e-8 1e-8 .01]);
[cases,~] = forecast_npi(truth,u); model = fit_npi_model(cumsum(cases),u,1e6,[3;3]);
[control,predicted,states] = optimal_npi(model,15,[1;1.5],.3);
verifySize(test_case,control,[2 15]); verifySize(test_case,states,[3 15]);
verifyGreaterThanOrEqual(test_case,control,0); verifyLessThanOrEqual(test_case,control,3);
verifyEqual(test_case,predicted(1),1e6*prod(model.state),'AbsTol',1e-9);
end
function test_cost_definition(test_case)
% TEST_COST_DEFINITION Cost averages over both interventions and days.
[human,cost] = npi_cost([10 20],[1 2;2 3],[1;2]); verifyEqual(test_case,human,15); verifyEqual(test_case,cost,3.25);
end
function test_diffusion_and_reflection(test_case)
% TEST_DIFFUSION_AND_REFLECTION Periodic diffusion conserves mass and reflections stay bounded.
x = zeros(21); x(11,11) = 1; result = diffusion_2d(x,1,.2,1,40);
verifyEqual(test_case,squeeze(sum(sum(result,1),2)),ones(41,1),'AbsTol',1e-12); verifyGreaterThanOrEqual(test_case,result,0);
motion = population_motion_2d([.2 .3],[10 -20],1,3); verifyGreaterThanOrEqual(test_case,motion,0); verifyLessThanOrEqual(test_case,motion,1);
end
function test_neural_forward_pass(test_case)
% TEST_NEURAL_FORWARD_PASS Optional MATLAB layers reproduce supplied-weight formulas.
assumeTrue(test_case,exist('nnet.layer.Layer','class')==8,'Deep Learning Toolbox required.');
x = reshape(linspace(-2,2,21),1,1,1,[]); exp_activation = exp_layer(1,'exp'); exp_activation.Alpha = .7;
tanh_activation = my_tanh_layer(1,'tanh',1); tanh_activation.Alpha = .7;
verifyEqual(test_case,predict(exp_activation,x),exp(.7*x),'AbsTol',1e-12);
verifyEqual(test_case,predict(tanh_activation,x),.7*tanh(x/.7),'AbsTol',1e-12);
end
