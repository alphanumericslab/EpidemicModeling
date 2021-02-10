% a script for testing the prescription algorithm developed for the XPRIZE
% Pandemic Response Challenge

clear;
close all;
clc

START_DATE = "2020-01-01"; % start time
END_DATE = "2021-02-07"; % end time
LATEST_DATA_FILE = './../covid-policy-tracker/data/OxCGRT_latest.csv'; % The historic data file cloned from: https://github.com/OxCGRT/covid-policy-tracker/tree/master/data
GEO_FILE = "xprize-sample-data/countries_regions.csv"; % countries and regions to include
POPULATION_FILE = "xprize-sample-data/populations.csv"; % country and regional populations
% TRAINED_MODEL_PARAMS_FILE = "xprize-sample-data/prescription_trained_params_lasso.mat"; % file to log the trained model parameters
TRAINED_MODEL_PARAMS_FILE = "xprize-sample-data/prescription_trained_params_nonnegls.mat"; % file to log the trained model parameters
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
IP_MAXES = [3, 3, 2, 4, 2, 3, 2, 4, 2, 3, 2, 4]'; % Maximum values for the 12 phase NPI values defined in the oxford data set

IP_FILE = "prescriptions/robojudge_test_scenario.csv";
costs_file = "xprize-sample-data/uniform_random_costs.csv";
output_file = "xprize-sample-data/prescriptor2020-08-01_2020-08-04.csv";
ip_file_historic = "xprize-sample-data/2020-09-30_historical_ip.csv";

% Train the model parameters
train_model = 0;
if(train_model)
    TrainNPIPrescriptor(START_DATE, END_DATE, LATEST_DATA_FILE, GEO_FILE, POPULATION_FILE, INCLUDED_IP, IP_MAXES, TRAINED_MODEL_PARAMS_FILE);
end

% generate scenarios
dt = 1;
I0 = 10;
alpha_min = 0;
alpha_max = 1;
gamma = 1/(dt*7.0); % input to contact influence rate (inverse time)
% alpha0 = 0.03;
prob_contagion_after_Tdays = 0.01;
Tdays = 21 * dt;
beta = -log(prob_contagion_after_Tdays)/Tdays; % recovery rate from being contagious (inverse time)
R0 = 2.5; % An assumption during outbreak
alpha0 = beta + log(R0)/dt; % The logic is that in the SIalpha model, during the outbreak R0 = exp(dt*(alpha - beta)) and alpha = beta is the metastable threshold (R0 = 1) %1.0/N_population; % the per-population normalization is needed

NumNPI = length(IP_MAXES);
num_days_before_opt_control = 30;
num_days_during_opt_control = 120;
load(TRAINED_MODEL_PARAMS_FILE);
ComumnNames = TrainedModelParams(1, :);
TrainedModelParams = TrainedModelParams(2 : end, :);

% Random weights constant over time:
% npi_weights = rand(1, NumNPI);
% npi_weights = NumNPI*npi_weights/sum(npi_weights);
% npi_weights_day_wise = diag(npi_weights) * ones(NumNPI, num_days_before_opt_control + num_days_during_opt_control);

% Random weights over time
npi_weights = rand(NumNPI, num_days_before_opt_control + num_days_during_opt_control);
sm = sum(npi_weights, 1);
npi_weights_day_wise = NumNPI*npi_weights./sm(ones(1, NumNPI), :);

s_noise_std = 1e-8;
i_noise_std = 1e-8;
alpha_noise_std = 1e-9;
for k = 84 : 84% length(TrainedModelParams) % the first row has the names
    k
    CountryName = TrainedModelParams{k, 1};
    RegionName = TrainedModelParams{k, 2};
    N_population = TrainedModelParams{k, 3};
    b1 = TrainedModelParams{k, 4};
    a1 = TrainedModelParams{k, 5};
    b2 = TrainedModelParams{k, 6};
    a2 = TrainedModelParams{k, 7};
    
    % select parameter set from training phase
    %     a = a1;
    %     b = b1;
    a = a2;
    b = b2;
    
    % ENFORCE MONOTONICITY OF THE LASSO MODEL
    %     a(a < 0 ) = 0;
    
    i0 = I0/N_population;
    s0 = (N_population - I0)/N_population;
    
    % generate historic data
    u_historic = zeros(NumNPI, num_days_before_opt_control);
    [s_historic, i_historic, alpha_historic] = SIalpha_Controlled(u_historic, s0, i0, alpha0, IP_MAXES, alpha_min, alpha_max, gamma, a, b, beta, s_noise_std, i_noise_std, alpha_noise_std, num_days_before_opt_control, dt);
    
    % generate zero control scenario
    u_zero_control = zeros(NumNPI, num_days_during_opt_control);
    [s_zero_control, i_zero_control, alpha_zero_control] = SIalpha_Controlled(u_zero_control, s_historic(end), i_historic(end), alpha_historic(end), IP_MAXES, alpha_min, alpha_max, gamma, a, b, beta, s_noise_std, i_noise_std, alpha_noise_std, num_days_during_opt_control, dt);
    s_zero_control = cat(2, s_historic, s_zero_control);
    i_zero_control = cat(2, i_historic, i_zero_control);
    alpha_zero_control = cat(2, alpha_historic, alpha_zero_control);
    u_zero_control = cat(2, u_historic, u_zero_control);
    [J0_zero_control, J1_zero_control] = NPICost(s_zero_control.*i_zero_control.*alpha_zero_control, u_zero_control, npi_weights_day_wise);
    
    % generate full (max) control scenario
    u_full_control = repmat(IP_MAXES, 1, num_days_during_opt_control);
    [s_full_control, i_full_control, alpha_full_control] = SIalpha_Controlled(u_full_control, s_historic(end), i_historic(end), alpha_historic(end), IP_MAXES, alpha_min, alpha_max, gamma, a, b, beta, s_noise_std, i_noise_std, alpha_noise_std, num_days_during_opt_control, dt);
    s_full_control = cat(2, s_historic, s_full_control);
    i_full_control = cat(2, i_historic, i_full_control);
    alpha_full_control = cat(2, alpha_historic, alpha_full_control);
    u_full_control = cat(2, u_historic, u_full_control);
    [J0_full_control, J1_full_control] = NPICost(s_full_control.*i_full_control.*alpha_full_control, u_full_control, npi_weights_day_wise);
    
    % generate optimal control scenario
    % % %     I0 = 100;%N_population*i_historic(end); % initial number of cases
    min_cases = 10; % the minimum cases start date for processing each region/country
    params.dt = 1.0; % temporal time scale
    %     control_input = nan(NumNPI, num_days_during_opt_control); % The input is to be estimated
    %     control_input = u_zero_control; % The input is to be estimated
    %     u_zero_control(:, end - num_days_during_opt_control + 1 : end) = nan;
    control_input = [u_historic, nan(NumNPI, num_days_during_opt_control)];
    params.w = npi_weights; % NPI costs
    params.a = a; % input influence weight vector
    params.b = b; % input influence bias constant
    params.u_min = zeros(NumNPI, 1); % minimum input values
    params.u_max = IP_MAXES(:); % maximum input values according to the OXFORD dataset
    params.alpha_min = 0.0; % minimum alpha
    params.alpha_max = inf; % maximum alpha
    params.gamma = gamma;%1/(params.dt*7.0); % input to contact influence rate (inverse time)
    params.obs_type = 'NEWCASES'; % TOTALCASES or NEWCASES
    % see the following for the background on beta: https://www.cdc.gov/coronavirus/2019-ncov/hcp/duration-isolation.html
    %prob_contagion_after_Tdays = 0.01;
    %Tdays = 21*params.dt;
    params.beta = beta;%-log(prob_contagion_after_Tdays)/Tdays; % recovery rate from being contagious (inverse time)
    %R0 = 1.1; % An assumption during outbreak
    % % %     alpha0 = alpha_historic(end);%params.beta + log(R0)/params.dt; % The logic is that in the SIalpha model, during the outbreak R0 = exp(dt*(alpha - beta)) and alpha = beta is the metastable threshold (R0 = 1) %1.0/N_population; % the per-population normalization is needed
    params.sigma = 10000; % sigmoid function slope
    beta_ekf = .9;%0.9; % Observation noise update factor (set to 1 for no update)
    gamma_ekf = 0.995; % Kalman gain stability factor (set very close to 1, or equal to 1 to disable the feature)
    inv_monitor_len_ekf = 21; % Window length for innovations process whiteness monitoring
    order = 1; % 1 for standard EKF; 2 for second-order EKF
    q_alpha = 1e-2;
    lambda0 = 1;
    q_lambda = 10;
    %     Q_w = (params.dt)^2*diag([10.0*s_historic(end), 10.0*i_historic(end), q_alpha, q_lambda, q_lambda, q_lambda].^2); % Process noise covariance matrix
    Q_w = (params.dt)^2*diag([10.0*i0, 30.0*i0, q_alpha, q_lambda, q_lambda, q_lambda].^2); % Process noise covariance matrix
    w_bar = zeros(6, 1); % mean value of process noises
    v_bar = 0; % mean value of observation noise
    R_v = var(5.0e3/N_population); % Observation noise variance
    
    %     num_pareto_front_points = 1;
    %     human_npi_cost_factor = 0.5;
    num_pareto_front_points = 1000;
    human_npi_cost_factor = logspace(-9.0, 0, num_pareto_front_points);
    human_npi_cost_factor = cat(2, human_npi_cost_factor, linspace(0, 1, num_pareto_front_points));
    num_pareto_front_points = num_pareto_front_points * 2;
    
    J0_opt_control = zeros(1, num_pareto_front_points);
    J1_opt_control = zeros(1, num_pareto_front_points);
    %     controlled_indexes = length(s_zero_control) - num_days_during_opt_control + 1 : length(s_zero_control);
    controlled_indexes = 1 : length(s_zero_control);
    NewCases = s_zero_control(controlled_indexes) .* i_zero_control(controlled_indexes) .* alpha_zero_control(controlled_indexes);
    ConfirmedCases = cumsum(s_zero_control(controlled_indexes) .* i_zero_control(controlled_indexes) .* alpha_zero_control(controlled_indexes));
    S_MINUS = zeros(6, length(NewCases), length(human_npi_cost_factor));
    S_PLUS = zeros(6, length(NewCases), length(human_npi_cost_factor));
    S_SMOOTH = zeros(6, length(NewCases), length(human_npi_cost_factor));
    S_MINUS_bkw = zeros(6, length(NewCases), length(human_npi_cost_factor));
    S_PLUS_bkw = zeros(6, length(NewCases), length(human_npi_cost_factor));
    S_SMOOTH_bkw = zeros(6, length(NewCases), length(human_npi_cost_factor));
    for ll = 1 : length(human_npi_cost_factor)
        params.epsilon = human_npi_cost_factor(ll); % [0, 1]: 0 neglects NPI cost and 1 neglects human factor!
        %     s_init = [s_historic(end) ; i_historic(end) ; alpha0 ; lambda0 ; lambda0 ; lambda0]; % initial state vector
        s_init = [s0 ; i0 ; alpha0 ; lambda0 ; lambda0 ; lambda0]; % initial state vector
        Ps_init = 100.0*(params.dt)^2*diag([i0, i0, q_alpha, q_lambda, q_lambda, q_lambda].^2); % Covariance matrix of initial states
        % Set the finite horizon end points
        s_final = nan(6, 1);
        % % %         s_final(2) = 0; % zero infection rates
        s_final(4) = 0; % zero costates
        s_final(5) = 0; % zero costates
        s_final(6) = 0; % zero costates
        Ps_final = zeros(6);
        Ps_final(1:3, 1:3) = nan; % don't touch these end-point states during smoothing
        % % %         Ps_final(2, 2) = 1e-4; % zero infection rates
        %     Ps_final(4 : 6, 4 : 6) = 1e-6;
        % % %     Ps_final(1, 1) = nan; % don't touch these end-point states during smoothing
        % % %     Ps_final(2, 2) = nan; % don't touch these end-point states during smoothing
        % % %     Ps_final(3, 3) = nan; % don't touch these end-point states during smoothing
        Ps_final(4, 4) = 1e-3; % zero costates
        Ps_final(5, 5) = 1e-3; % zero costates
        Ps_final(6, 6) = 1e-3; % zero costates
        if(isequal(params.obs_type, 'TOTALCASES'))
            [~, u_opt_control, S_MINUS(:, :, ll), S_PLUS(:, :, ll), S_SMOOTH(:, :, ll), P_MINUS, P_PLUS, P_SMOOTH, ~, ~, rho] = SIAlphaModelEKFOptControlled(control_input, ConfirmedCases(:)', params, s_init, Ps_init, s_final, Ps_final, w_bar, v_bar, Q_w, R_v, beta_ekf, gamma_ekf, inv_monitor_len_ekf, order);
            % % %             % replace final nans with highly uncertain arbitrary values
            % % %             % before backward EKF
            % % %             s_final(1) = S_PLUS(1, end, ll);
            % % %             s_final(2) = S_PLUS(2, end, ll);
            % % %             s_final(3) = S_PLUS(3, end, ll);
            % % %             nan_indexes = isnan(Ps_final);
            % % %             Ps_final(nan_indexes) = 1e10;
            % % %             Ps_init(4, 4) = 1e5;
            % % %             Ps_init(5, 5) = 1e5;
            % % %             Ps_init(6, 6) = 1e5;
            % % %             [u_opt_control_bkw, S_MINUS_bkw(:, :, ll), S_PLUS_bkw(:, :, ll), S_SMOOTH_bkw(:, :, ll), P_MINUS_bkw, P_PLUS_bkw, P_SMOOTH_bkw, ~, ~, rho_bkw] = SIAlphaModelBackwardEKFOptControlled(control_input, ConfirmedCases(:)', params, s_init, Ps_init, s_final, Ps_final, w_bar, v_bar, Q_w, R_v, beta_ekf, gamma_ekf, inv_monitor_len_ekf, order);
        elseif(isequal(params.obs_type, 'NEWCASES'))
            [~, u_opt_control, S_MINUS(:, :, ll), S_PLUS(:, :, ll), S_SMOOTH(:, :, ll), P_MINUS, P_PLUS, P_SMOOTH, ~, ~, rho] = SIAlphaModelEKFOptControlled(control_input, NewCases(:)', params, s_init, Ps_init, s_final, Ps_final, w_bar, v_bar, Q_w, R_v, beta_ekf, gamma_ekf, inv_monitor_len_ekf, order);
            % % %             % replace final nans with highly uncertain arbitrary values
            % % %             % before backward EKF
            % % %             s_final(1) = S_PLUS(1, end, ll);
            % % %             s_final(2) = S_PLUS(2, end, ll);
            % % %             s_final(3) = S_PLUS(3, end, ll);
            % % %             nan_indexes = isnan(Ps_final);
            % % %             Ps_final(nan_indexes) = 1e10;
            % % %             Ps_init(4, 4) = 1e5;
            % % %             Ps_init(5, 5) = 1e5;
            % % %             Ps_init(6, 6) = 1e5;
            % % %             [u_opt_control_bkw, S_MINUS_bkw(:, :, ll), S_PLUS_bkw(:, :, ll), S_SMOOTH_bkw(:, :, ll), P_MINUS_bkw, P_PLUS_bkw, P_SMOOTH_bkw, ~, ~, rho_bkw] = SIAlphaModelBackwardEKFOptControlled(control_input, NewCases(:)', params, s_init, Ps_init, s_final, Ps_final, w_bar, v_bar, Q_w, R_v, beta_ekf, gamma_ekf, inv_monitor_len_ekf, order);
        end
        %         [J0_opt_control(ll), J1_opt_control(ll)] = NPICost(S_SMOOTH(1, num_days_before_opt_control + 1: end, ll).*S_SMOOTH(2, num_days_before_opt_control + 1: end, ll).*S_SMOOTH(3, num_days_before_opt_control + 1: end, ll), opt_control_input(:, num_days_before_opt_control + 1: end), npi_weights_day_wise);
        
        %         figure;
        %         plot(u_opt_control');
        %         grid
        %         close all
        
        [s_opt_control, i_opt_control, alpha_opt_control] = SIalpha_Controlled(u_opt_control(:, end - num_days_during_opt_control + 1 : end), s_historic(end), i_historic(end), alpha_historic(end), IP_MAXES, alpha_min, alpha_max, gamma, a, b, beta, s_noise_std, i_noise_std, alpha_noise_std, num_days_during_opt_control, dt);
        s_opt_control = cat(2, s_historic, s_opt_control);
        i_opt_control = cat(2, i_historic, i_opt_control);
        alpha_opt_control = cat(2, alpha_historic, alpha_opt_control);
        % % %         u_opt_control = cat(2, u_historic, u_opt_control);
        
        
        
        
        % Use the EKF estimates
        %         [J0_opt_control(ll), J1_opt_control(ll)] = NPICost(squeeze(S_SMOOTH(1, :, ll).*S_SMOOTH(2, :, ll).*S_SMOOTH(3, :, ll)), u_opt_control, npi_weights_day_wise);
        [J0_opt_control(ll), J1_opt_control(ll)] = NPICost(s_opt_control.*i_opt_control.*alpha_opt_control, u_opt_control, npi_weights_day_wise);
        
        
        
        
        
        
        if(0)
            figure
            subplot(211);
            plot(NewCases, 'k', 'linewidth', 2);
            hold on
            plot(S_MINUS(1, :, ll).*S_MINUS(2, :, ll).*S_MINUS(3, :, ll));
            plot(S_PLUS(1, :, ll).*S_PLUS(2, :, ll).*S_PLUS(3, :, ll));
            plot(S_SMOOTH(1, :, ll).*S_SMOOTH(2, :, ll).*S_SMOOTH(3, :, ll));
            grid
            subplot(212);
            plot(opt_control_input');
            grid
        end
        if(0)
            figure
            title([CountryName , ' ', RegionName]);
            for mm = 1 : 6
                subplot(6, 1, mm);
                plot(S_PLUS(mm, :, ll));
                hold on
                plot(S_SMOOTH(mm, :, ll));
                legend('S_PLUS', 'S_SMOOTH');
                % % %                 plot(S_PLUS_bkw(mm, :, ll));
                % % %                 plot(S_SMOOTH_bkw(mm, :, ll));
                % % %                 legend('S_PLUS', 'S_SMOOTH', 'S_PLUS_bkw', 'S_SMOOTH_bkw');
                grid
            end
        end
    end
    
    
    % generate random control scenario
    num_random_input_monte_carlo_runs = 500;
    J0 = zeros(1, num_random_input_monte_carlo_runs);
    J1 = zeros(1, num_random_input_monte_carlo_runs);
    for scenario = 1 : num_random_input_monte_carlo_runs
        u = zeros(NumNPI, num_days_during_opt_control);
        for jj = 1 : NumNPI
            if(scenario < num_random_input_monte_carlo_runs/2) % random over NPI, constant over time
                u(jj, :) = randi([0, IP_MAXES(jj)]);
            else
                for t = 1 : num_days_during_opt_control % random over NPI and time
                    u(jj, t) = randi([0, IP_MAXES(jj)]);
                end
            end
        end
        %         u = u + 0.1*rand(size(u)); % add noise to input
        [s_controlled, i_controlled, alpha_controlled] = SIalpha_Controlled(u, s_historic(end), i_historic(end), alpha_historic(end), IP_MAXES, alpha_min, alpha_max, gamma, a, b, beta, s_noise_std, i_noise_std, alpha_noise_std, num_days_during_opt_control, dt);
        
        s = [s_historic, s_controlled];
        i = [i_historic, i_controlled];
        alpha = [alpha_historic, alpha_controlled];
        u = cat(2, u_historic, u);
        
        [J0(scenario), J1(scenario)] = NPICost(s.*i.*alpha, u, npi_weights_day_wise);
        
        if(scenario == 1)
            figure
        end
        subplot(311);
        plot(N_population * s .* i .* alpha);
        if(scenario == 1)
            hold on
            grid on
            title([CountryName '-' RegionName]);
        end
        subplot(312);
        plot(N_population * i);
        if(scenario == 1)
            hold on
            grid on
        end
    end
    subplot(311);
    %     plot(N_population * squeeze(S_MINUS(1, :, :).*S_MINUS(2, :, :).*S_MINUS(3, :, :)), 'k');
    %     plot(N_population * squeeze(S_PLUS(1, :, :).*S_PLUS(2, :, :).*S_PLUS(3, :, :)), 'k');
    %     plot(N_population * squeeze(S_SMOOTH(1, :, :).*S_SMOOTH(2, :, :).*S_SMOOTH(3, :, :)), 'k');
    plot(N_population * s_opt_control.*i_opt_control.*alpha_opt_control, 'k', 'linewidth', 2);
    ylabel('New cases');
    subplot(312);
    %     plot(N_population * squeeze(S_MINUS(2, :, :)), 'k');
    %     plot(N_population * squeeze(S_PLUS(2, :, :)), 'k');
    %     plot(N_population * squeeze(S_SMOOTH(2, :, :)), 'k');
    plot(N_population * i_opt_control, 'k', 'linewidth', 2);
    ylabel('Active cases: I(t)');
    
    subplot(313);
    
    figure
    %     plot(log(J0), log(J1), 'bo');
    plot(J0, J1, 'bo');
    hold on
    grid
    %     plot(log(J0_zero_control), log(J1_zero_control), 'kx', 'MarkerSize',12);
    %     plot(log(J0_full_control), log(J1_full_control), 'gx', 'MarkerSize',12);
    %     plot(log(J0_opt_control), log(J1_opt_control), 'ro');
    plot(J0_zero_control, J1_zero_control, 'kx', 'MarkerSize',12);
    plot(J0_full_control, J1_full_control, 'gx', 'MarkerSize',12);
    plot(J0_opt_control, J1_opt_control, 'ro');
    %     axis square
    xlabel('Human factor');
    ylabel('NPI cost');
    axis square
    % Find optimal control
    %     [u_opt, S_MINUS, S_PLUS, S_SMOOTH, P_MINUS, P_PLUS, P_SMOOTH, K_GAIN, innovations, rho] = SIAlphaModelEKFOptControlled(u, x, params, s_init, Ps_init, s_final, Ps_final, w_bar, v_bar, Q_w, R_v, beta, gamma, inv_monitor_len, order)
    %
    
    %     close all
end


% PrescribeNPI(START_DATE, END_DATE, ip_file, costs_file, output_file);

% license('inuse') % list the used licenses (useful for building standalone packages)