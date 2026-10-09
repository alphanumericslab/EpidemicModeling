function results = train_predict_prescribe_npi(npi_weights, human_npi_cost_factor, start_train_date_str, end_train_date_str, start_regression_date_str, end_predict_prescribe_date_str, data_file, geo_file, populations_file, included_ip, npi_mins, npi_maxes, trained_model_params_file)
% TRAIN_PREDICT_PRESCRIBE_NPI Train and compare fixed versus optimized NPI forecasts.
% Inputs follow train_npi_prescriptor, with nonnegative intervention weights,
% human_npi_cost_factor in [0,1], regression dates, forecast endpoint, and
% explicit NPI bounds. Returns a cell array of per-geography result structs:
% fixed_inputs/cases/states, optimal_inputs/cases/states, and two cost pairs.
% Forecasting begins one day after the inclusive training endpoint.
% Author: Reza Sameni | Emory University
% Reference: Sameni (2022), doi:10.1109/JSTSP.2021.3129118.
bundle = train_npi_prescriptor(start_train_date_str,end_train_date_str,data_file,geo_file,populations_file,included_ip,npi_maxes,trained_model_params_file,start_regression_date_str);
horizon = days(datetime(end_predict_prescribe_date_str)-datetime(end_train_date_str));
assert(horizon >= 1, 'Forecast endpoint must follow training.'); results = cell(size(bundle.models));
for k = 1:numel(bundle.models)
    model = bundle.models{k}; model.params.u_min = npi_mins(:);
    fixed = repmat(model.last_input(:),1,horizon); [cases,states] = forecast_npi(model,fixed);
    [control,opt_cases,opt_states] = optimal_npi(model,horizon,npi_weights,human_npi_cost_factor);
    result.country_name = model.country_name; result.region_name = model.region_name;
    result.fixed_inputs = fixed; result.fixed_cases = cases; result.fixed_states = states;
    result.optimal_inputs = control; result.optimal_cases = opt_cases; result.optimal_states = opt_states;
    [j0,j1] = npi_cost(cases,fixed,npi_weights(:)); result.fixed_cost = [j0 j1];
    [j0,j1] = npi_cost(opt_cases,control,npi_weights(:)); result.optimal_cost = [j0 j1];
    results{k} = result;
end
end
