function TrainPredictPrescribeNPI(npi_weights, human_npi_cost_factor, start_train_date_str, end_train_date_str, end_predict_presscribe_date_str, data_file, geo_file, populations_file, included_IP, NPI_MAXES, trained_model_params_file)

plot_results = true; % plot the results or not

% Convert start date to number
start_train_date_chars = char(start_train_date_str);
dash_indexes = start_train_date_chars == '-';
start_train_date_chars(dash_indexes) = [];
start_train_date = str2double(string(start_train_date_chars));

% Convert end date to number
end_train_date_chars = char(end_train_date_str);
dash_indexes = end_train_date_chars == '-';
end_train_date_chars(dash_indexes) = [];
end_train_date = str2double(string(end_train_date_chars));

% Convert end of prediction/prescription date to number
end_predict_presscribe_date_chars = char(end_predict_presscribe_date_str);
dash_indexes = end_predict_presscribe_date_chars == '-';
end_predict_presscribe_date_chars(dash_indexes) = [];
end_predict_presscribe_date = str2double(string(end_predict_presscribe_date_chars));

% Assess: start_train_date <= end_train_date <= end_predict_presscribe_date
if(~(start_train_date <= end_train_date && end_train_date <= end_predict_presscribe_date))
    error('Invalid input times order.');
end

% parameters
% plot_figures = false; % plot per-region/country plots or not
min_cases = 10; % the minimum cases start date for processing each region/country
% % % ar_order = 24; % Autoregressive model
% % % ar_learninghistory = 120; % Autoregressive model estimation history look back
% % % predict_ahead_num_days = 0; % number of days to predict ahead
SmoothingWinLen = 7; % window length used for smoothing the new cases

% Read the list of country/regions file
import_opt = detectImportOptions(geo_file);
import_opt.VariableNamingRule = 'preserve';
geolist_to_study = readtable(geo_file, import_opt);
SelectedGeoIDs = strcat(string(geolist_to_study.CountryName), " ", string(geolist_to_study.RegionName));

% Read the regional populations file
import_opt = detectImportOptions(populations_file);
import_opt.VariableNamingRule = 'preserve';
all_populations_data = readtable(populations_file, import_opt);

% Read Oxford time-series data
import_opt = detectImportOptions(data_file);
import_opt.VariableNamingRule = 'preserve';
import_opt.VariableTypes(6:end) = {'double'}; % All except the first five columns are 'double' (correct if not correctly identified)
all_data = readtable(data_file, import_opt);
AllCountryNames = all_data.CountryName; % All country names
AllRegionNames = all_data.RegionName; % All region names

% Make AllGeoIDs code (A unique combination of region and country names)
AllGeoIDs = strcat(string(AllCountryNames), " ", string(AllRegionNames));

% Make a dictionary of country-regions
[CountryAndRegionList, ~, ~] = unique(AllGeoIDs, 'stable');

NumGeoLocations = length(CountryAndRegionList); % Number of country-region pairs

NumNPI = length(included_IP);

TrainedModelParams = [{'CountryName'}, {'RegionName'}, {'N_population'}, {'coef0'}, {'coef'}, {'coef0_2'}, {'coef_2'}];


% Construct a table for the results
% TableHeaders = [{'CountryName'}; {'RegionName'}; {'Bias'}; included_IP(:)];
% CountryName = [];
% RegionName = [];
% b = [];
% a0 = [];
% a1 = [];
% a2 = [];
% a3 = [];
% a4 = [];
% a5 = [];
% a6 = [];
% a7 = [];
% a8 = [];
% a9 = [];
% a10 = [];
% a11 = [];
% a12 = [];

% NPI weights (not used during training. just needed as inputs to the function)
for k = 1 : NumGeoLocations%index_ : index_%240 %122 : 125%225 %1 : NumGeoLocations
    if(find(SelectedGeoIDs == CountryAndRegionList(k)))
        disp([num2str(k) '- ' char(CountryAndRegionList(k))]);
        geoid_all_row_indexes = AllGeoIDs == CountryAndRegionList(k) & all_data.Date >= start_train_date & all_data.Date <= end_train_date;
        ConfirmedCases = all_data.ConfirmedCases(geoid_all_row_indexes); % Confirmed cases
        ConfirmedDeaths = all_data.ConfirmedDeaths(geoid_all_row_indexes); % Confirmed deaths
        population_index = find(strcat(all_populations_data.CountryName, " ", all_populations_data.RegionName) == CountryAndRegionList(k));
        N_population = all_populations_data.Population2020(population_index);
        %         CaseFatalityJHDB = all_populations_data.CaseFatalityJHDBFeb2021(population_index);
        CountryName = all_data.CountryName(find(geoid_all_row_indexes == 1, 1, 'first'));
        RegionName = all_data.RegionName(find(geoid_all_row_indexes == 1, 1, 'first'));
        
        % Intervention plans
        InterventionPlans = all_data{geoid_all_row_indexes, included_IP}; % Region/Country intervention plans
        % Replace all N/A IP with previous IP
        for j = 1 : size(InterventionPlans, 2)
            for i = 2 : size(InterventionPlans, 1)
                if(isnan(InterventionPlans(i, j)) && not(isnan(InterventionPlans(i - 1, j))))
                    InterventionPlans(i, j) = InterventionPlans(i - 1, j);
                end
            end
        end
        InterventionPlans(isnan(InterventionPlans)) = 0; % Replace any remaining N/A IP with no IP
        
        % Calculate the number of new cases
        NewCases = [nan; diff(ConfirmedCases)]; % calculate the new daily cases
        NewCases(NewCases < 0) = 0;
        if(length(NewCases) < 2)
            warning(strcat("Insufficient data for ", CountryAndRegionList(k), ". Skipping the country/region."));
            continue;
        end
        % Replace nans with 0 and the last nan with the previous number
        NewCasesRefined = NewCases;
        NewCasesRefined(isnan(NewCasesRefined)) = 0;
        if(isnan(NewCases(end))) % Just fill the last missing date (if any) with the previous day
            NewCasesRefined(end) = NewCasesRefined(end-1);
        end
        % Smooth the newcases time series
        NewCasesSmoothed = filter(ones(1, SmoothingWinLen), SmoothingWinLen, NewCasesRefined); % causal
        NewCasesSmoothedZeroLag = filtfilt(ones(1, round(SmoothingWinLen/2)), round(SmoothingWinLen/2), NewCasesRefined); % causal
        NewCasesSmoothedNormalized = NewCasesSmoothed / N_population;
        
        % A smooth version of the confirmed cases time series
        ConfirmedCasesSmoothed = cumsum(NewCasesSmoothed);
        ConfirmedCasesSmoothedNormalized = ConfirmedCasesSmoothed / N_population;
        
        % Calculate the number of new deaths
        NewDeaths = [nan; diff(ConfirmedDeaths)]; % calculate the new daily deaths
        NewDeaths(NewDeaths < 0) = 0;
        % Replace nans with 0 and the last nan with the previous number
        NewDeathsRefined = NewDeaths;
        NewDeathsRefined(isnan(NewDeathsRefined)) = 0;
        if(isnan(NewDeaths(end))) % Just fill the last missing date (if any) with the previous day
            NewDeathsRefined(end) = 0;
        end
        % Smooth the deathcases time series
        NewDeathsSmoothed = filter(ones(1, SmoothingWinLen), SmoothingWinLen, NewDeathsRefined); % causal
        ConfirmedDeathsSmoothed = cumsum(NewDeathsSmoothed);
        
        % 1- Apply an EKF to find first estimates of alpha (with zero control assumption)
        %     I0 = min_cases; % initial number of cases
        first_case_indexes = find(NewCasesSmoothed > 0, 7, 'first'); % find the first non-zero indexes (average over the first week)
        I0 = max(min_cases, mean(NewCasesSmoothed(first_case_indexes))); % initial number of cases
        params.dt = 1.0; % temporal time scale
        control_input = zeros(size(InterventionPlans')); % The first round assumes zero inputs
        params.w = nan; % No NPI costs when input is zero
        params.a = zeros(NumNPI, 1); % input influence weight vector (zero in the first round)
        params.b = 0; % input influence bias constant (zero in the first round)
        params.u_min = zeros(NumNPI, 1); % minimum input values
        params.u_max = NPI_MAXES(:); % maximum input values according to the OXFORD dataset
        params.alpha_min = 0.0; % minimum alpha
        params.alpha_max = inf; % maximum alpha
        params.epsilon = human_npi_cost_factor; % [0, 1]: 0 neglects NPI cost and 1 neglects human factor!
        params.gamma = 1/(params.dt*7.0); % input to contact influence rate (inverse time)
        params.obs_type = 'NEWCASES'; % TOTALCASES or NEWCASES
        
        % see the following for the background on beta: https://www.cdc.gov/coronavirus/2019-ncov/hcp/duration-isolation.html
        prob_contagion_after_Tdays = 0.01;
        Tdays = 21*params.dt;
        params.beta = -log(prob_contagion_after_Tdays)/Tdays; % recovery rate from being contagious (inverse time)
        
        R0 = 2.5; % An assumption during outbreak
        alpha0 = params.beta + log(R0)/params.dt; % The logic is that in the SIalpha model, during the outbreak R0 = exp(dt*(alpha - beta)) and alpha = beta is the metastable threshold (R0 = 1) %1.0/N_population; % the per-population normalization is needed
        
        params.sigma = 10000; % sigmoid function slope
        beta = .9;%0.9; % Observation noise update factor (set to 1 for no update)
        gamma = 0.995; % Kalman gain stability factor (set very close to 1, or equal to 1 to disable the feature)
        inv_monitor_len = 21; % Window length for innovations process whiteness monitoring
        order = 1; % 1 for standard EKF; 2 for second-order EKF
        q_alpha = 1e-2;
        %         lambda_0 = 1.0;
        %         q_lambda = 1e-3;
        %         Q_w = 1*(params.dt)^2*diag([10.0*I0/N_population, 10.0*I0/N_population, q_alpha, q_lambda, q_lambda, q_lambda].^2); % Process noise covariance matrix
        Q_w = (params.dt)^2*diag([10.0*I0/N_population, 30.0*I0/N_population, q_alpha].^2); % Process noise covariance matrix
        %         w_bar = zeros(6, 1); % mean value of process noises
        w_bar = zeros(3, 1); % mean value of process noises
        v_bar = 0; % mean value of observation noise
        %         s_init = [(N_population - I0)/N_population ; I0/N_population ; alpha0 ; lambda_0 ; lambda_0 ; lambda_0]; % initial state vector
        s_init = [(N_population - I0)/N_population ; I0/N_population ; alpha0]; % initial state vector
        %         Ps_init = 1.0*(params.dt)^2*diag([10.0*I0/N_population, 10.0*I0/N_population, 10*q_alpha, 10*q_lambda, 10*q_lambda, 10*q_lambda].^2); % Covariance matrix of initial states
        Ps_init = 100.0*(params.dt)^2*diag([I0/N_population, I0/N_population, q_alpha].^2); % Covariance matrix of initial states
        %         s_final = nan(6, 1); % Set the finite horizon end points (not required during first round)
        s_final = nan(3, 1); % Set the finite horizon end points (not required during first round)
        %         Ps_final = nan(6); % Set the finite horizon end points (not required during first round)
        Ps_final = nan(3); % Set the finite horizon end points (not required during first round)
        if(isequal(params.obs_type, 'TOTALCASES'))
            R_v = var((NewCasesSmoothedZeroLag - NewCasesRefined)/N_population); % Observation noise variance
            %             [~, S_MINUS, S_PLUS, S_SMOOTH, P_MINUS, P_PLUS, P_SMOOTH, ~, ~, rho] = NewCaseEKFEstimatorWithOptimalNPI(control_input, ConfirmedCasesSmoothedNormalized(:)', params, s_init, Ps_init, s_final, Ps_final, w_bar, v_bar, Q_w, R_v, beta, gamma, inv_monitor_len, order);
            [~, ~, S_MINUS, S_PLUS, S_SMOOTH, P_MINUS, P_PLUS, P_SMOOTH, ~, ~, rho] = SIAlphaModelEKF(control_input, ConfirmedCasesSmoothedNormalized(:)', params, s_init, Ps_init, s_final, Ps_final, w_bar, v_bar, Q_w, R_v, beta, gamma, inv_monitor_len, order);
        elseif(isequal(params.obs_type, 'NEWCASES'))
            R_v = var((NewCasesSmoothedZeroLag - NewCasesRefined)/N_population); % Observation noise variance
            %             [~, S_MINUS, S_PLUS, S_SMOOTH, P_MINUS, P_PLUS, P_SMOOTH, ~, ~, rho] = NewCaseEKFEstimatorWithOptimalNPI(control_input, NewCasesSmoothedNormalized(:)', params, s_init, Ps_init, s_final, Ps_final, w_bar, v_bar, Q_w, R_v, beta, gamma, inv_monitor_len, order);
            [~, ~, S_MINUS, S_PLUS, S_SMOOTH, P_MINUS, P_PLUS, P_SMOOTH, ~, ~, rho] = SIAlphaModelEKF(control_input, NewCasesSmoothedNormalized(:)', params, s_init, Ps_init, s_final, Ps_final, w_bar, v_bar, Q_w, R_v, beta, gamma, inv_monitor_len, order);
        end
        
        % 2- Apply LASSO over alpha
        x_data_train = NPI_MAXES(:, ones(size(InterventionPlans, 1), 1))' - InterventionPlans;
        y_data_train = S_SMOOTH(3, :)'; % alpha
        %         y_data_train = S_PLUS(3, :)'; % alpha
        
        REGRESSION_TYPE = 'NONNEGATIVELS'; % 'LASSO' or 'NONNEGATIVELS'
        if(isequal(REGRESSION_TYPE, 'LASSO'))
            [B_LASSO, FitInfo] = lasso(x_data_train, y_data_train, 'CV',10);
            idxLambda1SE = FitInfo.Index1SE;
            coef = B_LASSO(:, idxLambda1SE);
            coef0 = FitInfo.Intercept(idxLambda1SE);
        elseif(isequal(REGRESSION_TYPE, 'NONNEGATIVELS'))
            coef = lsqnonneg(x_data_train, y_data_train);
            coef0 = 0;
        end
        y_pred_lasso = x_data_train * coef + coef0;
        %     axTrace = lassoPlot(B_LASSO, FitInfo);
        
        %         TWO_ROUND_TRAINING = true;
        %         if(TWO_ROUND_TRAINING)
        % 3- Apply the SECOND round EKF to refine alpha based on real inputs
        %         I0 = max(1, round(N_population * S_SMOOTH(2, 1)));
        %         I0 = max(1, round(N_population * S_SMOOTH(2, 1)));
        params.a = coef; % input influence weight vector
        params.b = coef0; % input influence bias constant
        control_input = InterventionPlans'; % The second round works with real inputs
        %s_init = [(N_population - I0)/N_population ; I0/N_population ; S_SMOOTH(3, 1) ; S_SMOOTH(4, 1) ; S_SMOOTH(5, 1) ; S_SMOOTH(6, 1)]; % initial state vector
        if(isequal(params.obs_type, 'TOTALCASES'))
            R_v = 0.1 * var((NewCasesSmoothed - NewCasesRefined)/N_population); % Observation noise variance
            %             [~, S_MINUS, S_PLUS, S_SMOOTH, P_MINUS, P_PLUS, P_SMOOTH, ~, ~, rho] = NewCaseEKFEstimatorWithOptimalNPI(control_input, ConfirmedCasesSmoothedNormalized(:)', params, s_init, Ps_init, s_final, Ps_final, w_bar, v_bar, Q_w, R_v, beta, gamma, inv_monitor_len, order);
            [~, ~, S_MINUS2, S_PLUS2, S_SMOOTH2, P_MINUS2, P_PLUS2, P_SMOOTH2, ~, ~, rho2] = SIAlphaModelEKF(control_input, ConfirmedCasesSmoothedNormalized(:)', params, s_init, Ps_init, s_final, Ps_final, w_bar, v_bar, Q_w, R_v, beta, gamma, inv_monitor_len, 1);
        elseif(isequal(params.obs_type, 'NEWCASES'))
            R_v = 0.1 * var((NewCasesSmoothed - NewCasesRefined)/N_population); % Observation noise variance
            %         [~, S_MINUS2, S_PLUS2, S_SMOOTH2, P_MINUS2, P_PLUS2, P_SMOOTH2, ~, ~, rho2] = NewCaseEKFEstimatorWithOptimalNPI(control_input, NewCasesSmoothedNormalized(:)', params, s_init, Ps_init, s_final, Ps_final, w_bar, v_bar, Q_w, R_v, beta, gamma, inv_monitor_len, 1);
            [~, ~, S_MINUS2, S_PLUS2, S_SMOOTH2, P_MINUS2, P_PLUS2, P_SMOOTH2, ~, ~, rho2] = SIAlphaModelEKF(control_input, NewCasesSmoothedNormalized(:)', params, s_init, Ps_init, s_final, Ps_final, w_bar, v_bar, Q_w, R_v, beta, gamma, inv_monitor_len, 1);
            % % %             [~, S_MINUS_backward, S_PLUS_backward, S_SMOOTH_backward, P_MINUS_backward, P_PLUS_backward, P_SMOOTH_backward, ~, ~, rho_backward] = SIAlphaModelBackwardEKF(control_input, NewCasesSmoothedNormalized(:)', params, s_init, Ps_init, s_final, Ps_final, w_bar, v_bar, Q_w, R_v, beta, gamma, inv_monitor_len, order);
        end
        
        % 4- Apply second LASSO over refined alpha
        x_data_train2 = NPI_MAXES(:, ones(size(InterventionPlans, 1), 1))' - InterventionPlans;
        y_data_train2 = S_SMOOTH2(3, :)';
        %         y_data_train2 = S_PLUS2(3, :)';
        
        if(isequal(REGRESSION_TYPE, 'LASSO'))
            [B_LASSO2, FitInfo2] = lasso(x_data_train2, y_data_train2, 'CV',10);
            idxLambda1SE2 = FitInfo2.Index1SE;
            coef_2 = B_LASSO2(:, idxLambda1SE2);
            coef0_2 = FitInfo2.Intercept(idxLambda1SE2);
        elseif(isequal(REGRESSION_TYPE, 'NONNEGATIVELS'))
            coef_2 = lsqnonneg(x_data_train, y_data_train);
            coef0_2 = 0;
        end
        y_pred_lasso2 = x_data_train2 * coef_2 + coef0_2;
        %         axTrace2 = lassoPlot(B_LASSO2, FitInfo2);
        %         end
        %
        
        forecast_time = end_predict_presscribe_date - end_train_date;
        if(forecast_time > 0)
            params.a = coef_2; % input influence weight vector
            params.b = coef0_2; % input influence bias constant
            
            % 5- Forecast/prescription phase with fixed input
            IP = InterventionPlans';
            IP_last = IP(:, end);
            control_input = [IP, 0*IP_last(:, ones(1, forecast_time))]; % The second round works with real inputs
            %s_init = [(N_population - I0)/N_population ; I0/N_population ; S_SMOOTH(3, 1) ; S_SMOOTH(4, 1) ; S_SMOOTH(5, 1) ; S_SMOOTH(6, 1)]; % initial state vector
            if(isequal(params.obs_type, 'TOTALCASES'))
                R_v = 0.1 * var((NewCasesSmoothed - NewCasesRefined)/N_population); % Observation noise variance
                observations = [ConfirmedCasesSmoothedNormalized(:)', nan(1, forecast_time)];
                %             [~, S_MINUS, S_PLUS, S_SMOOTH, P_MINUS, P_PLUS, P_SMOOTH, ~, ~, rho] = NewCaseEKFEstimatorWithOptimalNPI(control_input, ConfirmedCasesSmoothedNormalized(:)', params, s_init, Ps_init, s_final, Ps_final, w_bar, v_bar, Q_w, R_v, beta, gamma, inv_monitor_len, order);
            elseif(isequal(params.obs_type, 'NEWCASES'))
                R_v = 0.1 * var((NewCasesSmoothed - NewCasesRefined)/N_population); % Observation noise variance
                observations = [NewCasesSmoothedNormalized(:)', nan(1, forecast_time)];
                %         [~, S_MINUS2, S_PLUS2, S_SMOOTH2, P_MINUS2, P_PLUS2, P_SMOOTH2, ~, ~, rho2] = NewCaseEKFEstimatorWithOptimalNPI(control_input, NewCasesSmoothedNormalized(:)', params, s_init, Ps_init, s_final, Ps_final, w_bar, v_bar, Q_w, R_v, beta, gamma, inv_monitor_len, 1);
                % % %             [~, S_MINUS_backward, S_PLUS_backward, S_SMOOTH_backward, P_MINUS_backward, P_PLUS_backward, P_SMOOTH_backward, ~, ~, rho_backward] = SIAlphaModelBackwardEKF(control_input, NewCasesSmoothedNormalized(:)', params, s_init, Ps_init, s_final, Ps_final, w_bar, v_bar, Q_w, R_v, beta, gamma, inv_monitor_len, order);
            end
            [zero_control_input, ~, S_MINUS5, S_PLUS5, S_SMOOTH5, P_MINUS5, P_PLUS5, P_SMOOTH5, ~, ~, rho5] = SIAlphaModelEKF(control_input, observations, params, s_init, Ps_init, s_final, Ps_final, w_bar, v_bar, Q_w, R_v, beta, gamma, inv_monitor_len, 1);

            % 6- Forecast/prescription phase with fixed input
            IP = InterventionPlans';
            IP_last = IP(:, end);
            control_input = [IP, IP_last(:, ones(1, forecast_time))]; % The second round works with real inputs
            %s_init = [(N_population - I0)/N_population ; I0/N_population ; S_SMOOTH(3, 1) ; S_SMOOTH(4, 1) ; S_SMOOTH(5, 1) ; S_SMOOTH(6, 1)]; % initial state vector
            if(isequal(params.obs_type, 'TOTALCASES'))
                R_v = 0.1 * var((NewCasesSmoothed - NewCasesRefined)/N_population); % Observation noise variance
                observations = [ConfirmedCasesSmoothedNormalized(:)', nan(1, forecast_time)];
                %             [~, S_MINUS, S_PLUS, S_SMOOTH, P_MINUS, P_PLUS, P_SMOOTH, ~, ~, rho] = NewCaseEKFEstimatorWithOptimalNPI(control_input, ConfirmedCasesSmoothedNormalized(:)', params, s_init, Ps_init, s_final, Ps_final, w_bar, v_bar, Q_w, R_v, beta, gamma, inv_monitor_len, order);
            elseif(isequal(params.obs_type, 'NEWCASES'))
                R_v = 0.1 * var((NewCasesSmoothed - NewCasesRefined)/N_population); % Observation noise variance
                observations = [NewCasesSmoothedNormalized(:)', nan(1, forecast_time)];
                %         [~, S_MINUS2, S_PLUS2, S_SMOOTH2, P_MINUS2, P_PLUS2, P_SMOOTH2, ~, ~, rho2] = NewCaseEKFEstimatorWithOptimalNPI(control_input, NewCasesSmoothedNormalized(:)', params, s_init, Ps_init, s_final, Ps_final, w_bar, v_bar, Q_w, R_v, beta, gamma, inv_monitor_len, 1);
                % % %             [~, S_MINUS_backward, S_PLUS_backward, S_SMOOTH_backward, P_MINUS_backward, P_PLUS_backward, P_SMOOTH_backward, ~, ~, rho_backward] = SIAlphaModelBackwardEKF(control_input, NewCasesSmoothedNormalized(:)', params, s_init, Ps_init, s_final, Ps_final, w_bar, v_bar, Q_w, R_v, beta, gamma, inv_monitor_len, order);
            end
            [fixed_control_input, ~, S_MINUS3, S_PLUS3, S_SMOOTH3, P_MINUS3, P_PLUS3, P_SMOOTH3, ~, ~, rho3] = SIAlphaModelEKF(control_input, observations, params, s_init, Ps_init, s_final, Ps_final, w_bar, v_bar, Q_w, R_v, beta, gamma, inv_monitor_len, 1);
            
            % 6- Forecast/prescription phase with optimal input
            params.w = npi_weights;
            control_input = [InterventionPlans', nan(size(InterventionPlans, 2), forecast_time)]; % The second round works with real inputs
            %s_init = [(N_population - I0)/N_population ; I0/N_population ; S_SMOOTH(3, 1) ; S_SMOOTH(4, 1) ; S_SMOOTH(5, 1) ; S_SMOOTH(6, 1)]; % initial state vector
            if(isequal(params.obs_type, 'TOTALCASES'))
                R_v = 0.1 * var((NewCasesSmoothed - NewCasesRefined)/N_population); % Observation noise variance
                observations = [ConfirmedCasesSmoothedNormalized(:)', nan(1, forecast_time)];
                %             [~, S_MINUS, S_PLUS, S_SMOOTH, P_MINUS, P_PLUS, P_SMOOTH, ~, ~, rho] = NewCaseEKFEstimatorWithOptimalNPI(control_input, ConfirmedCasesSmoothedNormalized(:)', params, s_init, Ps_init, s_final, Ps_final, w_bar, v_bar, Q_w, R_v, beta, gamma, inv_monitor_len, order);
            elseif(isequal(params.obs_type, 'NEWCASES'))
                R_v = 0.1 * var((NewCasesSmoothed - NewCasesRefined)/N_population); % Observation noise variance
                observations = [NewCasesSmoothedNormalized(:)', nan(1, forecast_time)];
                %         [~, S_MINUS2, S_PLUS2, S_SMOOTH2, P_MINUS2, P_PLUS2, P_SMOOTH2, ~, ~, rho2] = NewCaseEKFEstimatorWithOptimalNPI(control_input, NewCasesSmoothedNormalized(:)', params, s_init, Ps_init, s_final, Ps_final, w_bar, v_bar, Q_w, R_v, beta, gamma, inv_monitor_len, 1);
                % % %             [~, S_MINUS_backward, S_PLUS_backward, S_SMOOTH_backward, P_MINUS_backward, P_PLUS_backward, P_SMOOTH_backward, ~, ~, rho_backward] = SIAlphaModelBackwardEKF(control_input, NewCasesSmoothedNormalized(:)', params, s_init, Ps_init, s_final, Ps_final, w_bar, v_bar, Q_w, R_v, beta, gamma, inv_monitor_len, order);
            end
            
            lambda0 = 1;
            q_lambda = 10;
            %     s_init = [s_historic(end) ; i_historic(end) ; alpha0 ; lambda0 ; lambda0 ; lambda0]; % initial state vector
            s_init = cat(1, s_init , [lambda0 ; lambda0 ; lambda0]); % initial state vector
            Q_w = blkdiag(Q_w, 10.0*(params.dt)^2*diag([q_lambda, q_lambda, q_lambda].^2)); % Covariance matrix of initial states
            Ps_init = blkdiag(Ps_init, 100.0*(params.dt)^2*diag([q_lambda, q_lambda, q_lambda].^2)); % Covariance matrix of initial states
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
            
            %             [opt_control_input, ~, S_MINUS3, S_PLUS3, S_SMOOTH3, P_MINUS3, P_PLUS3, P_SMOOTH3, ~, ~, rho3] = SIAlphaModelEKF(control_input, observations, params, s_init, Ps_init, s_final, Ps_final, w_bar, v_bar, Q_w, R_v, beta, gamma, inv_monitor_len, 1);
            [opt_control_input, opt_control_input_smooth, S_MINUS4, S_PLUS4, S_SMOOTH4, P_MINUS4, P_PLUS4, P_SMOOTH4, K_GAIN4, innovations4, rho4] = SIAlphaModelEKFOptControlled(control_input, observations, params, s_init, Ps_init, s_final, Ps_final, w_bar, v_bar, Q_w, R_v, beta, gamma, inv_monitor_len, 1);
        end
        % % % %         % Forecast future alpha values
        % % % %         if(TWO_ROUND_TRAINING)
        % % %             ar_train_segment = y_data_train2(end - ar_learninghistory + 1 : end);
        % % % %         else
        % % % %             ar_train_segment = y_data_train(end - ar_learninghistory + 1 : end);
        % % % %         end
        % % %         ar_sys = ar(ar_train_segment, ar_order);
        % % %         A_ar = get(ar_sys, 'A');
        % % %         noisevar_ar = get(ar_sys, 'NoiseVariance');
        % % %         zi = filtic(sqrt(noisevar_ar), A_ar, ar_train_segment(end:-1:1));
        % % %         y_pred_ar = filter(sqrt(noisevar_ar), A_ar, randn(1, predict_ahead_num_days), zi)';
        % % %         AlphaHatARX = [ar_train_segment ; y_pred_ar]';
        % % %         AlphaHatARX(AlphaHatARX < 0 ) = 0;
        
        MortalityRate = ConfirmedDeathsSmoothed ./ ConfirmedCasesSmoothed;
        MortalityRate(isnan(MortalityRate)) = 0;
%         MedMortalityRate = median(MortalityRate);
%         MedRecentMortalityRate = median(MortalityRate(round(0.75*end) : end));
        
%         CumInfections = cumsum(N_population * S_SMOOTH2(2, :));
%         DeathToCumInfectionsRatio = ConfirmedDeaths ./ CumInfections(:);
%         
%         BetaEstimate = DeathToCumInfectionsRatio/MedRecentMortalityRate;
%         MedRecentBetaEstimate = median(BetaEstimate(round(0.75*end) : end));
        
        %         TableData = cat(2, TableData, [ThisCountryName; ThisRegionName; num2cell(coef0); num2cell(coef(:))]);
        %         Age = [38;43;38;40;49];
        %         Smoker = logical([1;0;1;0;1]);
        %         Height = [71;69;64;67;64];
        %         Weight = [176;163;131;133;119];
        %         BloodPressure = [124 93; 109 77; 125 83; 117 75; 122 80];
        %
        %
        %         OutTable = table
        
        figure
        subplot(211);
        plot(N_population * S_PLUS5(1, :).*S_PLUS5(2, :).*S_PLUS5(3, :));
        hold on
        plot(N_population * S_PLUS3(1, :).*S_PLUS3(2, :).*S_PLUS3(3, :));
        plot(N_population * S_PLUS4(1, :).*S_PLUS4(2, :).*S_PLUS4(3, :));
        grid
        legend('new cases - no NPI', 'new cases - fixed NPI', 'new cases - optimal NPI');
        title(CountryAndRegionList(k), 'interpreter', 'none');
        subplot(212);
        plot(S_PLUS5(3, :));
        hold on
        plot(S_PLUS3(3, :));
        plot(S_PLUS4(3, :));
        legend('alpha - no NPI', 'alpha - fixed NPI', 'alpha - optimal NPI');
        grid
        
        pause(2)
        
        if(plot_results)
            figure
            subplot(411)
            plot(N_population * NewCasesSmoothedNormalized, 'linewidth', 2);
            hold on
            plot(N_population * S_PLUS(1, :).*S_PLUS(2, :).*S_PLUS(3, :));
            plot(N_population * S_SMOOTH(1, :).*S_SMOOTH(2, :).*S_SMOOTH(3, :));
            plot(N_population * S_PLUS2(1, :).*S_PLUS2(2, :).*S_PLUS2(3, :));
            plot(N_population * S_SMOOTH2(1, :).*S_SMOOTH2(2, :).*S_SMOOTH2(3, :));
            legend('NewCases', 'PLUS', 'SMOOTH', 'PLUS2', 'SMOOTH2');
            title(CountryAndRegionList(k), 'interpreter', 'none');
            grid;
            
            subplot(412);
            plot(rho)
            hold on
            plot(rho2, 'r')
            legend('rho', 'rho2');
            grid
            
            subplot(413);
            plot(y_data_train);
            hold on;
            plot(y_data_train2);
            plot(y_pred_lasso);
            plot(y_pred_lasso2);
            grid
            legend('alpha', 'alpha 2', 'alpha lasso', 'alpha lasso 2');
            
            subplot(414);
            plot(InterventionPlans);
            grid
            % % %         plot(AlphaHatARX)
            % % %         hold on
            % % %         plot(ar_train_segment, 'r')
            % % %         grid
            
            % % % %         figure
            % % % % %         plot(NewCasesSmoothed);
            % % % %         plot(ConfirmedCases);
            % % % %         hold on
            % % % %         plot(ConfirmedDeaths);
            % % % %         plot(N_population * S_PLUS(1, :) .* S_PLUS(2, :) .* S_PLUS(3, :));
            % % % %         plot(cumsum(N_population * S_PLUS(2, :)));
            % % % %         grid
            % % % %         legend('ConfirmedCases', 'ConfirmedDeaths', 'Estimated new cases', 'Estimated infections');
            % % % % %         legend('NewCasesSmoothed', 'ConfirmedDeaths', 'Estimated new cases', 'Estimated infections');
            
            % % %         figure
            % % %         subplot(211);
            % % %         plot(100.0 * MortalityRate);
            % % %         hold on
            % % %         plot(100.0 * MedMortalityRate(ones(1, length(MortalityRate))), 'linewidth', 2);
            % % %         plot(100.0 * MedRecentMortalityRate(ones(1, length(MortalityRate))), 'linewidth', 2);
            % % %         plot(CaseFatalityJHDB(ones(1, length(MortalityRate))), 'linewidth', 2);
            % % %         title([CountryAndRegionList(k) ' mortality rate (%)'], 'interpreter', 'none');
            % % %         grid
            % % %         legend('MortalityRate', 'MedMortalityRate', 'MedRecentMortalityRate', 'CaseFatalityJHDB');
            % % %
            % % %         subplot(212);
            % % %         plot(BetaEstimate);
            % % %         hold on
            % % %         plot(MedRecentBetaEstimate(ones(1, length(BetaEstimate))), 'linewidth', 2);
            % % %         legend('\beta estimate', 'MedRecentBetaEstimate');
            % % %         grid
            
            % % %         figure
            % % %         hold on
            % % %         plot(diff(ConfirmedCases));
            % % %         plot(NewCasesSmoothed);
            % % %         plot(NewDeathsSmoothed/MortalityRate(end));
            % % %         %     plot(ConfirmedCases - ConfirmedDeaths);
            % % %         %     legend('ConfirmedCases', 'ConfirmedDeaths/MortalityRate', 'RecoveredCases');
            % % %         legend('diff(ConfirmedCases)', 'NewCasesSmoothed', 'NewDeathsSmoothedNormalized');
            % % %         title([CountryAndRegionList(k) ' mortality rate (%)'], 'interpreter', 'none');
            % % %         grid
            
            
            figure
            subplot(311);
            errorbar(N_population*S_MINUS(1, :), N_population*sqrt(squeeze(P_MINUS(1, 1, :))));
            hold on
            errorbar(N_population*S_PLUS(1, :), N_population*sqrt(squeeze(P_PLUS(1, 1, :))));
            errorbar(N_population*S_SMOOTH(1, :), N_population*sqrt(squeeze(P_SMOOTH(1, 1, :))));
            errorbar(N_population*S_MINUS2(1, :), N_population*sqrt(squeeze(P_MINUS2(1, 1, :))));
            errorbar(N_population*S_PLUS2(1, :), N_population*sqrt(squeeze(P_PLUS2(1, 1, :))));
            errorbar(N_population*S_SMOOTH2(1, :), N_population*sqrt(squeeze(P_SMOOTH2(1, 1, :))));
            %         errorbar(N_population*S_MINUS_backward(1, :), N_population*sqrt(squeeze(P_MINUS_backward(1, 1, :))));
            %         errorbar(N_population*S_PLUS_backward(1, :), N_population*sqrt(squeeze(P_PLUS_backward(1, 1, :))));
            %         errorbar(N_population*S_SMOOTH_backward(1, :), N_population*sqrt(squeeze(P_SMOOTH_backward(1, 1, :))));
            legend('MINUS', 'PLUS', 'SMOOTH', 'MINUS2', 'PLUS2', 'SMOOTH2');%, 'MINUS_backward', 'PLUS_backward', 'SMOOTH_backward');
            title(CountryAndRegionList(k), 'interpreter', 'none');
            grid
            subplot(312);
            errorbar(N_population*S_MINUS(2, :), N_population*sqrt(squeeze(P_MINUS(2, 2, :))));
            hold on
            errorbar(N_population*S_PLUS(2, :), N_population*sqrt(squeeze(P_PLUS(2, 2, :))));
            errorbar(N_population*S_SMOOTH(2, :), N_population*sqrt(squeeze(P_SMOOTH(2, 2, :))));
            %         plot(diff(ConfirmedDeathsSmoothed)/(MedRecentMortalityRate*params.beta), 'linewidth', 2);
            errorbar(N_population*S_MINUS2(2, :), N_population*sqrt(squeeze(P_MINUS2(2, 2, :))));
            errorbar(N_population*S_PLUS2(2, :), N_population*sqrt(squeeze(P_PLUS2(2, 2, :))));
            errorbar(N_population*S_SMOOTH2(2, :), N_population*sqrt(squeeze(P_SMOOTH2(2, 2, :))));
            %         errorbar(N_population*S_MINUS_backward(2, :), N_population*sqrt(squeeze(P_MINUS_backward(2, 2, :))));
            %         errorbar(N_population*S_PLUS_backward(2, :), N_population*sqrt(squeeze(P_PLUS_backward(2, 2, :))));
            %         errorbar(N_population*S_SMOOTH_backward(2, :), N_population*sqrt(squeeze(P_SMOOTH_backward(2, 2, :))));
            plot(diff(ConfirmedDeathsSmoothed)./(MortalityRate(2:end)*params.beta), 'linewidth', 2);
            legend('MINUS', 'PLUS', 'SMOOTH', 'MINUS2', 'PLUS2', 'SMOOTH2', 'diff(ConfirmedDeathsSmoothed)/(MedRecentMortalityRate*params.beta)');
            title(CountryAndRegionList(k), 'interpreter', 'none');
            grid
            subplot(313);
            errorbar(S_MINUS(3, :), sqrt(squeeze(P_MINUS(3, 3, :))));
            hold on
            errorbar(S_PLUS(3, :), sqrt(squeeze(P_PLUS(3, 3, :))));
            errorbar(S_SMOOTH(3, :), sqrt(squeeze(P_SMOOTH(3, 3, :))));
            errorbar(S_MINUS2(3, :), sqrt(squeeze(P_MINUS2(3, 3, :))));
            errorbar(S_PLUS2(3, :), sqrt(squeeze(P_PLUS2(3, 3, :))));
            errorbar(S_SMOOTH2(3, :), sqrt(squeeze(P_SMOOTH2(3, 3, :))));
            %         errorbar(S_MINUS_backward(3, :), sqrt(squeeze(P_MINUS_backward(3, 3, :))));
            %         errorbar(S_PLUS_backward(3, :), sqrt(squeeze(P_PLUS_backward(3, 3, :))));
            %         errorbar(S_SMOOTH_backward(3, :), sqrt(squeeze(P_SMOOTH_backward(3, 3, :))));
            legend('MINUS', 'PLUS', 'SMOOTH', 'MINUS2', 'PLUS2', 'SMOOTH2');%, 'MINUS_backward', 'PLUS_backward', 'SMOOTH_backward');
            title(CountryAndRegionList(k), 'interpreter', 'none');
            grid
            
            % % %         figure
            % % %         subplot(311);
            % % %         errorbar(S_MINUS(4, :), sqrt(squeeze(P_MINUS(4, 4, :))));
            % % %         hold on
            % % %         errorbar(S_PLUS(4, :), sqrt(squeeze(P_PLUS(4, 4, :))));
            % % %         errorbar(S_SMOOTH(4, :), sqrt(squeeze(P_SMOOTH(4, 4, :))));
            % % %         errorbar(S_MINUS2(4, :), sqrt(squeeze(P_MINUS2(1, 1, :))));
            % % %         errorbar(S_PLUS2(4, :), sqrt(squeeze(P_PLUS2(4, 4, :))));
            % % %         errorbar(S_SMOOTH2(4, :), sqrt(squeeze(P_SMOOTH2(4, 4, :))));
            % % %         legend('MINUS', 'PLUS', 'SMOOTH', 'MINUS2', 'PLUS2', 'SMOOTH2');
            % % %         title(CountryAndRegionList(k), 'interpreter', 'none');
            % % %         grid
            % % %         subplot(312);
            % % %         errorbar(S_MINUS(5, :), sqrt(squeeze(P_MINUS(5, 5, :))));
            % % %         hold on
            % % %         errorbar(S_PLUS(5, :), sqrt(squeeze(P_PLUS(5, 5, :))));
            % % %         errorbar(S_SMOOTH(5, :),sqrt(squeeze(P_SMOOTH(5, 5, :))));
            % % %         errorbar(S_MINUS2(5, :), sqrt(squeeze(P_MINUS2(5, 5, :))));
            % % %         errorbar(S_PLUS2(5, :), sqrt(squeeze(P_PLUS2(5, 5, :))));
            % % %         errorbar(S_SMOOTH2(5, :), sqrt(squeeze(P_SMOOTH2(5, 5, :))));
            % % %         legend('MINUS', 'PLUS', 'SMOOTH', 'MINUS2', 'PLUS2', 'SMOOTH2');
            % % %         title(CountryAndRegionList(k), 'interpreter', 'none');
            % % %         grid
            % % %         subplot(313);
            % % %         errorbar(S_MINUS(6, :), sqrt(squeeze(P_MINUS(6, 6, :))));
            % % %         hold on
            % % %         errorbar(S_PLUS(6, :), sqrt(squeeze(P_PLUS(6, 6, :))));
            % % %         errorbar(S_SMOOTH(6, :), sqrt(squeeze(P_SMOOTH(6, 6, :))));
            % % %         errorbar(S_MINUS2(6, :), sqrt(squeeze(P_MINUS2(6, 6, :))));
            % % %         errorbar(S_PLUS2(6, :), sqrt(squeeze(P_PLUS2(6, 6, :))));
            % % %         errorbar(S_SMOOTH2(6, :), sqrt(squeeze(P_SMOOTH2(6, 6, :))));
            % % %         legend('MINUS', 'PLUS', 'SMOOTH', 'MINUS2', 'PLUS2', 'SMOOTH2');
            % % %         title(CountryAndRegionList(k), 'interpreter', 'none');
            % % %         grid
            
            % % %             params.beta
            % % %             N_population
            
            close all
        end
        
        TrainedModelParams = cat(1, TrainedModelParams, [CountryName, RegionName, {N_population}, {coef0}, {coef}, {coef0_2}, {coef_2}]);
    end
end
save(char(trained_model_params_file), 'TrainedModelParams');

