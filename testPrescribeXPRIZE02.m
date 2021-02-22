% A script for testing the prescription algorithm developed for the XPRIZE
% Pandemic Response Challenge
% (C) Reza Sameni, 2021
% reza.sameni@gmail.com
% https://github.com/alphanumericslab/EpidemicModeling

clear;
close all;
clc

START_DATE_TRAIN = "2020-01-01"; % start time
START_DATE_REGRESSION = "2020-01-15"; % start NPI model regression date
END_DATE_TRAIN = "2021-02-07"; % end time
END_DATE_PREDICT_PRESCRIBE = "2021-06-07"; % end time
% LATEST_DATA_FILE = './../covid-policy-tracker/data/OxCGRT_latest.csv'; % The historic data file cloned from: https://github.com/OxCGRT/covid-policy-tracker/tree/master/data
LATEST_DATA_FILE = 'xprize-sample-data/OxCGRT_latest.csv'; % The historic data file cloned from: https://github.com/OxCGRT/covid-policy-tracker/tree/master/data
GEO_FILE = "xprize-sample-data/countries_regions.csv"; % countries and regions to include
% GEO_FILE = "xprize-sample-data/countries_regions_short_list.csv"; % countries and regions to include
POPULATION_FILE = "xprize-sample-data/populations.csv"; % country and regional populations
% TRAINED_MODEL_PARAMS_FILE = "xprize-sample-data/prescription_trained_params_lasso.mat"; % file to log the trained model parameters
% TRAINED_MODEL_PARAMS_FILE = "xprize-sample-data/prescription_trained_params_nonnegls.mat"; % file to log the trained model parameters
TRAINED_MODEL_PARAMS_FILE = "xprize-sample-data/prescription_trained_params_train_and_prescribe.mat"; % file to log the trained model parameters
INCLUDED_IP = {'C1_School closing',... % see full descriptions at: https://github.com/OxCGRT/covid-policy-tracker/blob/master/documentation/codebook.md%'C1_Flag',...
    'C2_Workplace closing',... %'C2_Flag',...
    'C3_Cancel public events',... %'C3_Flag',...
    'C4_Restrictions on gatherings',... %'C4_Flag',...
    'C5_Close public transport',... %'C5_Flag',...
    'C6_Stay at home requirements',... %'C6_Flag',...
    'C7_Restrictions on internal movement',... 'C7_Flag',...
    'C8_International travel controls',... %'E1_Income support', 'E1_Flag', 'E2_Debt/contract relief',... % 'E3_Fiscal measures', 'E4_International support',...
    'H1_Public information campaigns',... %'H1_Flag',...
    'H2_Testing policy',...
    'H3_Contact tracing',... % 'H4_Emergency investment in healthcare', 'H5_Investment in vaccines',...
    'H6_Facial Coverings',... %'H6_Flag',... % {'H7_Vaccination policy', 'H7_Flag',...%'Holidays'...
    };
IP_MINS = zeros(12, 1); % Minimum values for the 12 phase NPI values defined in the oxford data set
IP_MAXES = [3, 3, 2, 4, 2, 3, 2, 4, 2, 3, 2, 4]'; % Maximum values for the 12 phase NPI values defined in the oxford data set
NumNPI = length(IP_MAXES);
num_days_before_opt_control = 30;
num_days_during_opt_control = 120;

IP_FILE = "prescriptions/robojudge_test_scenario.csv";
costs_file = "xprize-sample-data/uniform_random_costs.csv";
output_file = "xprize-sample-data/prescriptor2020-08-01_2020-08-04.csv";
ip_file_historic = "xprize-sample-data/2020-09-30_historical_ip.csv";

% human_npi_cost_factor = 1e-6;  % [0, 1]: 0 neglects NPI cost and 1 neglects human factor!
num_pareto_front_points = 250;
% human_npi_cost_factor = logspace(-12.0, -1e-8, num_pareto_front_points/2);
% human_npi_cost_factor = cat(2, human_npi_cost_factor, linspace(1e-8, 1-1e-8, num_pareto_front_points/2));
human_npi_cost_factor = logspace(-12.0, -eps, num_pareto_front_points/2);
human_npi_cost_factor = cat(2, human_npi_cost_factor, linspace(eps, 1-eps, num_pareto_front_points/2));

% Equal weights constant over time:
npi_weights = ones(1, NumNPI);
npi_weights = NumNPI*npi_weights/sum(npi_weights);
npi_weights_day_wise = diag(npi_weights) * ones(NumNPI, num_days_before_opt_control + num_days_during_opt_control);

% Random weights constant over time:
% npi_weights = rand(1, NumNPI);
% npi_weights = NumNPI*npi_weights/sum(npi_weights);
% npi_weights_day_wise = diag(npi_weights) * ones(NumNPI, num_days_before_opt_control + num_days_during_opt_control);

% Random weights over time
% npi_weights = rand(NumNPI, num_days_before_opt_control + num_days_during_opt_control);
% sm = sum(npi_weights, 1);
% npi_weights_day_wise = NumNPI*npi_weights./sm(ones(1, NumNPI), :);

% Train/Predict/Prescribe
TrainPredictPrescribeNPI(npi_weights, human_npi_cost_factor, START_DATE_TRAIN, END_DATE_TRAIN, START_DATE_REGRESSION, END_DATE_PREDICT_PRESCRIBE, LATEST_DATA_FILE, GEO_FILE, POPULATION_FILE, INCLUDED_IP, IP_MINS, IP_MAXES, TRAINED_MODEL_PARAMS_FILE);
