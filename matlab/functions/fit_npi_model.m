function model = fit_npi_model(cumulative, inputs, population, u_max, regression_start)
% FIT_NPI_MODEL Fit a two-pass EKF/NNLS SI-alpha model to one geography.
% model = fit_npi_model(cumulative,inputs,population,u_max,regression_start)
% cumulative: N counts; inputs: interventions-by-N; population: positive count;
% u_max: intervention bounds; regression_start: one-based start (default 1).
% Output struct contains params, coefficients, final filtered state/covariance,
% population and training_days. Python uses a zero-based regression_start.
% The second regression uses the refined alpha series; see docs/migration.md.
% Author: Reza Sameni | Emory University
if nargin < 5, regression_start = 1; end
cumulative = cumulative(:)'; u_max = u_max(:); count = numel(cumulative);
assert(count >= 14 && population > 0 && isequal(size(inputs),[numel(u_max) count]), 'Invalid training dimensions.');
assert(regression_start >= 1 && regression_start <= count, 'Invalid regression_start.');
assert(all(isfinite(inputs(:))) && all(inputs >= 0,'all') && all(inputs <= u_max,'all'), 'Inputs outside bounds.');
[daily, smoothed] = prepare_cases(cumulative);
positives = smoothed(smoothed > 0); positives = positives(1:min(7,numel(positives)));
if isempty(positives), initial = 10; else, initial = max(10,mean(positives)); end
assert(initial < population, 'Population must exceed initial cases.');
fraction = initial/population; params = default_si_params(u_max);
state = [1-fraction; fraction; params.beta+log(2.5)];
covariance = 100*diag([fraction fraction .01].^2);
Q = diag([10*fraction 30*fraction .01].^2);
variance = max(var((daily-smoothed)/population),1e-14);
[~,~,~,~,first] = si_alpha_model_ekf(zeros(size(inputs)),smoothed/population,params,state,covariance,nan(3,1),nan(3),zeros(3,1),0,Q,variance,1,1,21,1);
design = (u_max-inputs)';
coef = nonnegative_least_squares(design(regression_start:end,:),first(3,regression_start:end)');
params.a = coef;
[~,~,~,second_plus,second,~,second_cov] = si_alpha_model_ekf(inputs,smoothed/population,params,state,covariance,nan(3,1),nan(3),zeros(3,1),0,Q,variance,1,1,21,1);
refined = nonnegative_least_squares(design(regression_start:end,:),second(3,regression_start:end)');
params.a = refined;
model.population = population; model.params = params;
model.coefficients = coef; model.refined_coefficients = refined;
model.state = second_plus(:,end); model.covariance = second_cov(:,:,end); model.training_days = count;
end
