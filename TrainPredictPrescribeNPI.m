function TrainPredictPrescribeNPI(npi_weights, human_npi_cost_factor, start_train_date_str, end_train_date_str, start_regression_date_str, end_predict_presscribe_date_str, data_file, geo_file, populations_file, included_IP, NPI_MINS, NPI_MAXES, trained_model_params_file)

% A function for the prediction and prescription algorithm developed for the XPRIZE
% Pandemic Response Challenge
% (C) Reza Sameni, 2021
% reza.sameni@gmail.com
% https://github.com/alphanumericslab/EpidemicModeling

% npi_weights: NPI weights; not used during training; only used during prescription phase.

% parameters
plot_results = true; % plot per-region/country plots or not
SmoothingWinLen = 7; % window length used for smoothing the new cases
min_cases = 1; % the absolute minimum number of cases at any time
first_num_days_for_case_estimation = 7; % the first few days used to find an initial estimate of the first cases (only for EKF initialization)
model_gamma_param = 7; % input to contact influence time constant
observation_type = 'NEWCASES'; % TOTALCASES or NEWCASES
num_days_for_beta_calculation = 21;
prob_contagion_after_Tdays = 0.01;
R0 = 2.5; % An assumption during outbreak (https://doi.org/10.1016/S1473-3099(20)30484-9)
REGRESSION_TYPE = 'NONNEGATIVELS'; % 'LASSO'/'NONNEGATIVELS'/'NONNEGATIVELS-ELEMENT-WISE'
NONNEGATIVELS_IRERATIONS = 100;

% Convert training start date string to number
start_train_date_chars = char(start_train_date_str);
start_train_date = datetime(str2double(start_train_date_chars(1:4)), str2double(start_train_date_chars(6:7)), str2double(start_train_date_chars(9:10)));
dash_indexes = start_train_date_chars == '-';
start_train_date_chars(dash_indexes) = [];
start_train_date_number = str2double(string(start_train_date_chars));

% Convert training end date string to number
end_train_date_chars = char(end_train_date_str);
end_train_date = datetime(str2double(end_train_date_chars(1:4)), str2double(end_train_date_chars(6:7)), str2double(end_train_date_chars(9:10)));
dash_indexes = end_train_date_chars == '-';
end_train_date_chars(dash_indexes) = [];
end_train_date_number = str2double(string(end_train_date_chars));

% Convert LASSO/LSTM/... regression start date string to number
start_regression_date_chars = char(start_regression_date_str);
start_regression_date = datetime(str2double(start_regression_date_chars(1:4)), str2double(start_regression_date_chars(6:7)), str2double(start_regression_date_chars(9:10)));
% dash_indexes = start_regression_date_chars == '-';
% start_regression_date_chars(dash_indexes) = [];
% start_regression_date_number = str2double(string(start_regression_date_chars));

% Convert end of prediction/prescription date string to number
end_predict_presscribe_date_chars = char(end_predict_presscribe_date_str);
end_predict_date = datetime(str2double(end_predict_presscribe_date_chars(1:4)), str2double(end_predict_presscribe_date_chars(6:7)), str2double(end_predict_presscribe_date_chars(9:10)));
dash_indexes = end_predict_presscribe_date_chars == '-';
end_predict_presscribe_date_chars(dash_indexes) = [];
end_predict_presscribe_date_number = str2double(string(end_predict_presscribe_date_chars));

% find the number of days for forecasting
dates = [start_train_date, start_regression_date, end_train_date, end_predict_date];
num_forecast_days = split(between(dates(3), dates(4), 'days'), 'd');
num_regression_days = split(between(dates(2), dates(3), 'days'), 'd');

% Assess: start_train_date_number <= end_train_date_number <= end_predict_presscribe_date_number
if(~(start_train_date_number <= end_train_date_number && end_train_date_number <= end_predict_presscribe_date_number))
    error('Invalid input times order.');
end

% Read the list of country/regions file to be processed
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
all_data = readtable(data_file, import_opt); % The big data matrix (all historic cases and NPIs)
AllCountryNames = all_data.CountryName; % All country names
AllRegionNames = all_data.RegionName; % All region names

% Make AllGeoIDs code (A unique combination of region and country names)
AllGeoIDs = strcat(string(AllCountryNames), " ", string(AllRegionNames));

% Make a dictionary of country-regions
[CountryAndRegionList, ~, ~] = unique(AllGeoIDs, 'stable');

NumGeoLocations = length(CountryAndRegionList); % Number of country-region pairs

NumNPI = length(included_IP); % Number of NPI

TrainedModelParams = [{'CountryName'}, {'RegionName'}, {'N_population'}, {'reg_coef_b'}, {'reg_coef_a'}, {'reg_coef_b2'}, {'reg_coef_a2'}];

for k = 1 : NumGeoLocations
    if(~isempty(find(SelectedGeoIDs == CountryAndRegionList(k), 1)) && isequal(CountryAndRegionList(k), "Argentina ")) % make sure the current country/region is among the ones to be processed
        % 1) READ AND CLEAN THE DATA
        disp([num2str(k) '- ' char(CountryAndRegionList(k))]);
        geoid_all_row_indexes = AllGeoIDs == CountryAndRegionList(k) & all_data.Date >= start_train_date_number & all_data.Date <= end_train_date_number; % fetch all rows corresponding to the country/region from start date to end date
        CountryName = all_data.CountryName(find(geoid_all_row_indexes == 1, 1, 'first')); % The current country name
        RegionName = all_data.RegionName(find(geoid_all_row_indexes == 1, 1, 'first')); % The current region name
        ConfirmedCases = all_data.ConfirmedCases(geoid_all_row_indexes); % Confirmed cases
        ConfirmedDeaths = all_data.ConfirmedDeaths(geoid_all_row_indexes); % Confirmed deaths
        
        % Get coutry/region population
        population_index = find(strcat(all_populations_data.CountryName, " ", all_populations_data.RegionName) == CountryAndRegionList(k));
        N_population = all_populations_data.Population2020(population_index);
        
        % Intervention plans
        InterventionPlans = all_data{geoid_all_row_indexes, included_IP}; % Region/Country intervention plans
        NumNPIdays = size(InterventionPlans, 1); % Number of days with NPI
        
        % Preprocess the N/A IPs
        for j = 1 : NumNPI % Replace all N/A IP with previous IP
            for i = 2 : NumNPIdays
                if(isnan(InterventionPlans(i, j)) && not(isnan(InterventionPlans(i - 1, j))))
                    InterventionPlans(i, j) = InterventionPlans(i - 1, j);
                end
            end
        end
        InterventionPlans(isnan(InterventionPlans)) = 0; % Replace any remaining N/A IP with no IP
        
        % Calculate the number of new cases
        NewCases = diff([ConfirmedCases(1) ; ConfirmedCases]); % calculate the new daily cases
        
        % Preprocess the new cases for trivial issues
        NewCases(NewCases < 0) = 0;
        if(length(NewCases) < 2)
            warning(strcat("Insufficient data for ", CountryAndRegionList(k), ". Skipping the country/region."));
            continue;
        end
        % Replace the last nan with the last valid number
        NewCasesRefined = NewCases;
        if(isnan(NewCases(end))) % Just fill the last missing date (if any) with the previous day
            last_number = find(not(isnan(NewCases)), 1, 'last');
            NewCasesRefined(end) = NewCasesRefined(last_number);
        end
        NewCasesRefined(isnan(NewCasesRefined)) = 0; % replace any remaining nans with 0
        
        % Smooth the newcases time series
        NewCasesSmoothed = filter(ones(1, SmoothingWinLen), SmoothingWinLen, NewCasesRefined); % causal smoothed; used for processing
        NewCasesSmoothedZeroLag = filtfilt(ones(1, round(SmoothingWinLen/2)), round(SmoothingWinLen/2), NewCasesRefined); % non-causal; used for noise variance calculation only
        NewCasesSmoothedNormalized = NewCasesSmoothed / N_population; % the fraction of new cases
        
        % A smooth version of the confirmed cases time series
        ConfirmedCasesSmoothed = cumsum(NewCasesSmoothed);
        ConfirmedCasesSmoothedZeroLag = cumsum(NewCasesSmoothedZeroLag);
        ConfirmedCasesSmoothedNormalized = ConfirmedCasesSmoothed / N_population; % the fraction of total cases
        
        % Calculate the number of new deaths
        NewDeaths = diff([ConfirmedDeaths(1) ; ConfirmedDeaths]); % calculate the new daily deaths
        NewDeaths(NewDeaths < 0) = 0; % people are not bord due to covid19!
        NewDeathsRefined = NewDeaths;
        if(isnan(NewDeaths(end))) % Just fill the last missing date (if any) with the previous day
            last_number = find(not(isnan(NewDeaths)), 1, 'last');
            NewDeathsRefined(end) = NewDeathsRefined(last_number);
        end
        NewDeathsRefined(isnan(NewDeathsRefined)) = 0; % replace any remaining nans with 0
        
        % Smooth the deathcases time series
        NewDeathsSmoothed = filter(ones(1, SmoothingWinLen), SmoothingWinLen, NewDeathsRefined); % causal smoothed; used for processing
        ConfirmedDeathsSmoothed = cumsum(NewDeathsSmoothed); % A smooth version of the confirmed deaths time series
        
        FatalityRate = ConfirmedDeathsSmoothed ./ ConfirmedCasesSmoothed;
        FatalityRate(isnan(FatalityRate)) = 0;
        
        % 2) APPLY A TRI-STATE EKF TO FIND A FIRST ESTIMATE OF ALPHA (ZERO CONTROL ASSUMPTION)
        first_case_indexes = find(NewCasesSmoothed > 0, first_num_days_for_case_estimation, 'first'); % find the first non-zero indexes over the first_num_days_for_case_estimation
        I0 = max(min_cases, mean(NewCasesSmoothed(first_case_indexes))); % initial number of cases
        params.dt = 1.0; % temporal time scale
        control_input = zeros(NumNPI, NumNPIdays); % The first round assumes zero inputs
        params.w = nan; % No NPI costs when input is zero
        params.a = zeros(NumNPI, 1); % input influence weight vector (zero in the first round)
        params.b = 0; % input influence bias constant (zero in the first round)
        params.u_min = NPI_MINS(:); % minimum input values
        params.u_max = NPI_MAXES(:); % maximum input values according to the OXFORD dataset
        params.s_min = min_cases/N_population; % minimum s
        params.i_min = min_cases/N_population; % minimum i
        
        params.alpha_min = 1e-8; % minimum alpha
        params.alpha_max = 100; % maximum alpha
        
        params.epsilon = nan; % No NPI costs during the first round
        params.gamma = 1/(params.dt * model_gamma_param); % input to contact influence rate (inverse time)
        params.obs_type = observation_type;
        
        % see the following for the background on beta: https://www.cdc.gov/coronavirus/2019-ncov/hcp/duration-isolation.html
        Tdays = num_days_for_beta_calculation * params.dt;
        params.beta = -log(prob_contagion_after_Tdays)/Tdays; % recovery rate from being contagious (inverse time)
        alpha0 = params.beta + log(R0)/params.dt; % The logic is that in the SIalpha model, during the outbreak R0 = exp(dt*(alpha - beta)) and alpha = beta is the metastable threshold (R0 = 1) %1.0/N_population; % the per-population normalization is needed
        
        params.sigma = 1000000; % sigmoid function slope
        beta_ekf = 1.0;%.9; % Observation noise update factor (set to 1 for no update)
        gamma_ekf = 0.995; % Kalman gain stability factor (set very close to 1, or equal to 1 to disable the feature)
        inv_monitor_len_ekf = 21; % Window length for innovations process whiteness monitoring
        order = 1; % 1 for standard EKF; 2 for second-order EKF
        s_noise_std = 10.0*I0/N_population;
        i_noise_std = 30.0*I0/N_population;
        alpha_noise_std = 1e-2;
        Q_w = (params.dt)^2*diag([s_noise_std, i_noise_std, alpha_noise_std].^2); % Process noise covariance matrix
        w_bar = zeros(3, 1); % mean value of process noises
        v_bar = 0; % mean value of observation noise
        s_init = [(N_population - I0)/N_population ; I0/N_population ; alpha0]; % initial state vector
        %         Ps_init = 100.0*(params.dt)^2*diag([I0/N_population, I0/N_population, alpha_noise_std].^2); % Covariance matrix of initial states
        Ps_init = (params.dt)^2*diag([10*s_noise_std, 10*i_noise_std, 10*alpha_noise_std].^2); % Covariance matrix of initial states
        s_final = nan(3, 1); % Set the finite horizon end points (not required during first round)
        Ps_final = nan(3); % Set the finite horizon end points (not required during first round)
        R_v = 0.1 * ((NewCasesSmoothedZeroLag' - NewCasesRefined')/N_population).^2; % Observation noise variance estimate
        if(isequal(params.obs_type, 'TOTALCASES'))
            observations = ConfirmedCasesSmoothedNormalized(:)';
        elseif(isequal(params.obs_type, 'NEWCASES'))
            observations = NewCasesSmoothedNormalized(:)';
            %             R_v = var((NewCasesSmoothedZeroLag - NewCasesRefined)/N_population); % Observation noise variance estimate
            %             R_v = BaseLine1(((NewCasesSmoothedZeroLag' - NewCasesRefined')/N_population).^2, SmoothingWinLen, 'mn'); % Observation noise variance estimate
        end
        [~, ~, S_MINUS_round1_noinput, S_PLUS_round1_noinput, S_SMOOTH_round1_noinput, P_MINUS_round1_noinput, P_PLUS_round1_noinput, P_SMOOTH_round1_noinput, ~, ~, rho_round1_noinput] = SIAlphaModelEKF(control_input, observations, params, s_init, Ps_init, s_final, Ps_final, w_bar, v_bar, Q_w, R_v, beta_ekf, gamma_ekf, inv_monitor_len_ekf, order);
        
        % 3) APPLY LASSO OVER ALPHA ESTIMATE TO FIND INITIAL LINEAR REGRESSION MODEL PARAMS
        x_data_train = NPI_MAXES(:, ones(NumNPIdays, 1))' - InterventionPlans;
        y_data_train = S_SMOOTH_round1_noinput(3, :)'; % alpha
        % y_data_train = S_PLUS_round1_noinput(3, :)'; % alpha
        if(isequal(REGRESSION_TYPE, 'LASSO'))
            [B_LASSO, FitInfo] = lasso(x_data_train(end-num_regression_days+1:end, :), y_data_train(end-num_regression_days+1:end), 'CV',50);
            lassoPlot(B_LASSO,FitInfo,'PlotType','CV');
            legend('show') % Show legend
            % idxLambda1SE = FitInfo.Index1SE;
            idxLambda1 = FitInfo.IndexMinMSE;
            reg_coef_a = B_LASSO(:, idxLambda1);
            %             reg_coef_a(reg_coef_a < 0) = 0;
            reg_coef_b = FitInfo.Intercept(idxLambda1);
        elseif(isequal(REGRESSION_TYPE, 'NONNEGATIVELS'))
            reg_coef_a = lsqnonneg(x_data_train(end-num_regression_days+1:end, :), y_data_train(end-num_regression_days+1:end));
            reg_coef_b = 0;
            min_err = sum((y_data_train(end-num_regression_days+1:end) - x_data_train(end-num_regression_days+1:end, :) * reg_coef_a).^2);
            for jj = 1 : NONNEGATIVELS_IRERATIONS
                coef_temp = lsqnonneg(x_data_train(end-num_regression_days+1:end, :), y_data_train(end-num_regression_days+1:end) - reg_coef_b);
                coef0_temp = mean(y_data_train(end-num_regression_days+1:end) - x_data_train(end-num_regression_days+1:end, :) * reg_coef_a);
                error_temp = sum((y_data_train(end-num_regression_days+1:end) - x_data_train(end-num_regression_days+1:end, :) * reg_coef_a - coef0_temp).^2);
                if(error_temp < min_err)
                    reg_coef_a = coef_temp;
                    reg_coef_b = coef0_temp;
                    min_err = error_temp;
                else
                    break;
                end
            end
        elseif(isequal(REGRESSION_TYPE, 'NONNEGATIVELS-ELEMENT-WISE'))
            % Affine
            S = fitoptions('Method', 'NonlinearLeastSquares', 'Robust', 'on', 'Lower', [0, -inf], 'Upper', [inf, inf], 'MaxIter', 10000, 'Startpoint', [0, 0]);
            f = fittype('a * x + b', 'options', S);
            %             % Quadratic
            %             S = fitoptions('Method', 'NonlinearLeastSquares', 'Robust', 'on', 'Lower', [0, 0, -inf], 'Upper', [inf, inf, inf], 'MaxIter', 10000, 'Startpoint', [0, 0, 0]);
            %             f = fittype('a*x^2 + b * x + c', 'options', S);
            reg_coef_a = zeros(NumNPI, 1);
            for kk = 1 : length(reg_coef_a)
                ffit = fit(x_data_train(end-num_regression_days+1:end, kk), y_data_train(end-num_regression_days+1:end), f);
                reg_coef_a(kk) = ffit.a;
            end
            reg_coef_b = mean(y_data_train(end-num_regression_days+1:end) - x_data_train(end-num_regression_days+1:end, :) * reg_coef_a);
        end
        y_pred_lasso = x_data_train * reg_coef_a + reg_coef_b;
        
        % 4) APPLY THE SECOND ROUND TRI-STATE EKF TO REFINE ALPHA BASED ON REAL HISTORIC NPI INPUTS
        params.a = reg_coef_a; % input influence weight vector
        params.b = reg_coef_b; % input influence bias constant
        control_input = InterventionPlans'; % The second round works with real inputs
        %         I0 = max(1, round(N_population * S_SMOOTH_round1_noinput(2, 1)));
        %         I0 = max(1, round(N_population * S_SMOOTH_round1_noinput(2, 1)));
        %s_init = [(N_population - I0)/N_population ; I0/N_population ; S_SMOOTH_round1_noinput(3, 1) ; S_SMOOTH_round1_noinput(4, 1) ; S_SMOOTH_round1_noinput(5, 1) ; S_SMOOTH_round1_noinput(6, 1)]; % initial state vector
        if(isequal(params.obs_type, 'TOTALCASES'))
            observations = ConfirmedCasesSmoothedNormalized(:)';
        elseif(isequal(params.obs_type, 'NEWCASES'))
            observations = NewCasesSmoothedNormalized(:)';
        end
        [~, ~, S_MINUS_round2_withinput, S_PLUS_round2_withinput, S_SMOOTH_round2_withinput, P_MINUS_round2_withinput, P_PLUS_round2_withinput, P_SMOOTH_round2_withinput, ~, ~, rho_round2_withinput] = SIAlphaModelEKF(control_input, observations, params, s_init, Ps_init, s_final, Ps_final, w_bar, v_bar, Q_w, R_v, beta_ekf, gamma_ekf, inv_monitor_len_ekf, 1);
        
        % 5) APPLY SECOND LASSO OVER REFINED ALPHA
        x_data_train2 = NPI_MAXES(:, ones(NumNPIdays, 1))' - InterventionPlans;
        y_data_train2 = S_SMOOTH_round2_withinput(3, :)';
        %         y_data_train2 = S_PLUS_round2_withinput(3, :)';
        if(isequal(REGRESSION_TYPE, 'LASSO'))
            [B_LASSO2, FitInfo2] = lasso(x_data_train2(end-num_regression_days+1:end, :), y_data_train2(end-num_regression_days+1:end), 'CV',50);
            %             idxLambda = FitInfo2.Index1SE;
            lassoPlot(B_LASSO2,FitInfo2,'PlotType','CV');
            legend('show') % Show legend
            idxLambda2 = FitInfo2.IndexMinMSE;
            reg_coef_a2 = B_LASSO2(:, idxLambda2);
            
            %             reg_coef_a2(reg_coef_a2 < 0) = 0;
            
            reg_coef_b2 = FitInfo2.Intercept(idxLambda2);
        elseif(isequal(REGRESSION_TYPE, 'NONNEGATIVELS'))
            reg_coef_a2 = lsqnonneg(x_data_train2(end-num_regression_days+1:end, :), y_data_train2(end-num_regression_days+1:end));
            reg_coef_b2 = 0;
            min_err = sum((y_data_train2(end-num_regression_days+1:end) - x_data_train2(end-num_regression_days+1:end, :) * reg_coef_a2).^2);
            for jj = 1 : NONNEGATIVELS_IRERATIONS
                coef_2_temp = lsqnonneg(x_data_train2(end-num_regression_days+1:end, :), y_data_train2(end-num_regression_days+1:end) - reg_coef_b2);
                coef0_2_temp = mean(y_data_train2(end-num_regression_days+1:end) - x_data_train2(end-num_regression_days+1:end, :) * reg_coef_a2);
                error_temp = sum((y_data_train2(end-num_regression_days+1:end) - x_data_train2(end-num_regression_days+1:end, :) * reg_coef_a2 - coef0_2_temp).^2);
                if(error_temp < min_err)
                    reg_coef_a2 = coef_2_temp;
                    reg_coef_b2 = coef0_2_temp;
                    min_err = error_temp;
                else
                    break;
                end
            end
        elseif(isequal(REGRESSION_TYPE, 'NONNEGATIVELS-ELEMENT-WISE'))
            % Affine
            S = fitoptions('Method', 'NonlinearLeastSquares', 'Robust', 'on', 'Lower', [0, -inf], 'Upper', [inf, inf], 'MaxIter', 10000, 'Startpoint', [0, 0]);
            f = fittype('a * x + b', 'options', S);
            %             % Quadratic
            %             S = fitoptions('Method', 'NonlinearLeastSquares', 'Robust', 'on', 'Lower', [0, 0, -inf], 'Upper', [inf, inf, inf], 'MaxIter', 10000, 'Startpoint', [0, 0, 0]);
            %             f = fittype('c*x^2 + a * x + b', 'options', S);
            reg_coef_a2 = zeros(NumNPI, 1);
            for kk = 1 : length(reg_coef_a2)
                ffit = fit(x_data_train2(end-num_regression_days+1:end, kk), y_data_train2(end-num_regression_days+1:end), f);
                reg_coef_a2(kk) = ffit.a;
            end
            reg_coef_b2 = mean(y_data_train2(end-num_regression_days+1:end) - x_data_train2(end-num_regression_days+1:end, :) * reg_coef_a2);
        end
        y_pred_lasso2 = x_data_train2 * reg_coef_a2 + reg_coef_b2;
        
        % 6) FORECAST/PRESCRIPTION PHASE
        if(num_forecast_days > 0)
            params.a = reg_coef_a2; % input influence weight vector
            params.b = reg_coef_b2; % input influence bias constant
            
            % A) Forecast/prescription phase with zero input (3-state EKF)
            IP = InterventionPlans';
            control_input = [IP, NPI_MINS(:, ones(1, num_forecast_days))]; % The second round works with real inputs
            R_v = cat(2, R_v, ones(1, num_forecast_days) * mean(R_v)); % Fill in the R_v matrix with the average variance
            %s_init = [(N_population - I0)/N_population ; I0/N_population ; S_SMOOTH_round1_noinput(3, 1) ; S_SMOOTH_round1_noinput(4, 1) ; S_SMOOTH_round1_noinput(5, 1) ; S_SMOOTH_round1_noinput(6, 1)]; % initial state vector
            if(isequal(params.obs_type, 'TOTALCASES'))
                observations = [ConfirmedCasesSmoothedNormalized(:)', nan(1, num_forecast_days)];
            elseif(isequal(params.obs_type, 'NEWCASES'))
                observations = [NewCasesSmoothedNormalized(:)', nan(1, num_forecast_days)];
            end
            [zero_control_input, ~, S_MINUS_zeroinput, S_PLUS_zeroinput, S_SMOOTH_zeroinput, P_MINUS_zeroinput, P_PLUS_zeroinput, P_SMOOTH_zeroinput, ~, ~, rho_zeroinput] = SIAlphaModelEKF(control_input, observations, params, s_init, Ps_init, s_final, Ps_final, w_bar, v_bar, Q_w, R_v, beta_ekf, gamma_ekf, inv_monitor_len_ekf, 1);
            
            % B) Forecast/prescription phase with fixed input (3-state EKF)
            IP = InterventionPlans';
            IP_last = IP(:, end);
            control_input = [IP, IP_last(:, ones(1, num_forecast_days))]; % The second round works with real inputs
            %s_init = [(N_population - I0)/N_population ; I0/N_population ; S_SMOOTH_round1_noinput(3, 1) ; S_SMOOTH_round1_noinput(4, 1) ; S_SMOOTH_round1_noinput(5, 1) ; S_SMOOTH_round1_noinput(6, 1)]; % initial state vector
            if(isequal(params.obs_type, 'TOTALCASES'))
                observations = [ConfirmedCasesSmoothedNormalized(:)', nan(1, num_forecast_days)];
            elseif(isequal(params.obs_type, 'NEWCASES'))
                observations = [NewCasesSmoothedNormalized(:)', nan(1, num_forecast_days)];
            end
            [fixed_control_input, ~, S_MINUS_fixedinput, S_PLUS_fixedinput, S_SMOOTH_fixedinput, P_MINUS_fixedinput, P_PLUS_fixedinput, P_SMOOTH_fixedinput, ~, ~, rho_fixedinput] = SIAlphaModelEKF(control_input, observations, params, s_init, Ps_init, s_final, Ps_final, w_bar, v_bar, Q_w, R_v, beta_ekf, gamma_ekf, inv_monitor_len_ekf, 1);
            
            % C) Forecast/prescription phase with optimal input (6-state EKF)
            num_pareto_front_points = length(human_npi_cost_factor);
            J0_opt_control = zeros(1, num_pareto_front_points);
            J1_opt_control = zeros(1, num_pareto_front_points);
            for ll = 1 : num_pareto_front_points
                disp(['Pareto front #' num2str(ll)]);
                params.w = npi_weights;
                params.epsilon = human_npi_cost_factor(ll); % [0, 1]: 0 neglects NPI cost and 1 neglects human factor!
                lambda0 = 1.0;
                q_lambda = 0.01;
                %     s_init = [s_historic(end) ; i_historic(end) ; alpha0 ; lambda0 ; lambda0 ; lambda0]; % initial state vector
                ss_init = cat(1, s_init , [lambda0 ; lambda0 ; lambda0]); % initial state vector
                QQ_w = blkdiag(Q_w, (params.dt)^2*diag([q_lambda, q_lambda, q_lambda].^2)); % Covariance matrix of initial states
                PPs_init = blkdiag(Ps_init, 10.0*(params.dt)^2*diag([q_lambda, q_lambda, q_lambda].^2)); % Covariance matrix of initial states
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
                Ps_final(4, 4) = 1e-8; % zero costates
                Ps_final(5, 5) = 1e-8; % zero costates
                Ps_final(6, 6) = 1e-8; % zero costates
                control_input = [InterventionPlans', nan(NumNPI, num_forecast_days)]; % The second round works with real inputs
                %s_init = [(N_population - I0)/N_population ; I0/N_population ; S_SMOOTH_round1_noinput(3, 1) ; S_SMOOTH_round1_noinput(4, 1) ; S_SMOOTH_round1_noinput(5, 1) ; S_SMOOTH_round1_noinput(6, 1)]; % initial state vector
                if(isequal(params.obs_type, 'TOTALCASES'))
                    observations = [ConfirmedCasesSmoothedNormalized(:)', nan(1, num_forecast_days)];
                elseif(isequal(params.obs_type, 'NEWCASES'))
                    observations = [NewCasesSmoothedNormalized(:)', nan(1, num_forecast_days)];
                end
                [opt_control_input, opt_control_input_smooth, S_MINUS_optinput, S_PLUS_optinput, S_SMOOTH_optinput, P_MINUS_optinput, P_PLUS_optinput, P_SMOOTH_optinput, ~, ~, rho_optinput] = SIAlphaModelEKFOptControlled(control_input, observations, params, ss_init, PPs_init, s_final, Ps_final, w_bar, v_bar, QQ_w, R_v, beta_ekf, gamma_ekf, inv_monitor_len_ekf, 1);
                
                % Backward filtering (under test)
                s_final = [1, 0, 0, 0, 0, 0];
                Ps_final = diag([.1, .1, 1e-2, 1e-8, 1e-8, 1e-8]);
                [~, ~, S_MINUS_backward, S_PLUS_backward, S_SMOOTH_backward, P_MINUS_backward, P_PLUS_backward, P_SMOOTH_backward, ~, ~, rho_backward] = SIAlphaModelBackwardEKFOptControlled(control_input, observations, params, ss_init, PPs_init, s_final, Ps_final, w_bar, v_bar, QQ_w, R_v, beta_ekf, gamma_ekf, inv_monitor_len_ekf, order);
                
                %                 S_FRW_BCK = zeros(6, length(observations));
                %                 P_FRW_BCK = zeros(6, 6, length(observations));
                %                 for hh = 1 : length(observations)
                %                     P_FRW_BCK(:, :, hh) = (P_PLUS_optinput(:, :, hh) + P_PLUS_backward(:, :, hh))\(P_PLUS_optinput(:, :, hh) * P_PLUS_backward(:, :, hh));
                %                     S_FRW_BCK(:, hh) = pinv(P_PLUS_optinput(:, :, hh) + P_PLUS_backward(:, :, hh))*(P_PLUS_backward(:, :, hh) * S_PLUS_optinput(:, hh) + P_PLUS_optinput(:, :, hh) * S_PLUS_backward(:, hh));
                %                 end
                
                % Generate test scenario
                s_historic = S_SMOOTH_optinput(1, 1 : NumNPIdays);
                i_historic = S_SMOOTH_optinput(2, 1 : NumNPIdays);
                alpha_historic = S_SMOOTH_optinput(3, 1 : NumNPIdays);
                [s_opt_control, i_opt_control, alpha_opt_control] = SIalpha_Controlled(opt_control_input_smooth(:, NumNPIdays + 1 : end), s_historic(end), i_historic(end), alpha_historic(end), NPI_MAXES, params.alpha_min, params.alpha_max, params.gamma, params.a, params.b, params.beta, s_noise_std, i_noise_std, alpha_noise_std, num_forecast_days, params.dt);
                s_opt_control = cat(2, s_historic, s_opt_control);
                i_opt_control = cat(2, i_historic, i_opt_control);
                alpha_opt_control = cat(2, alpha_historic, alpha_opt_control);
                npi_weights_day_wise = diag(npi_weights) * ones(NumNPI, size(opt_control_input_smooth, 2));
                % % %         u_opt_control = cat(2, u_historic, u_opt_control);
                % Use the EKF estimates
                %         [J0_opt_control(ll), J1_opt_control(ll)] = NPICost(squeeze(S_SMOOTH(1, :, ll).*S_SMOOTH(2, :, ll).*S_SMOOTH(3, :, ll)), u_opt_control, npi_weights_day_wise);
                [J0_opt_control(ll), J1_opt_control(ll)] = NPICost(s_opt_control.*i_opt_control.*alpha_opt_control, opt_control_input_smooth, npi_weights_day_wise);
            end
            
            % D) Generate random control scenario
            num_random_input_monte_carlo_runs = 500;
            J0 = zeros(1, num_random_input_monte_carlo_runs);
            J1 = zeros(1, num_random_input_monte_carlo_runs);
            for scenario = 1 : num_random_input_monte_carlo_runs
                u = zeros(NumNPI, num_forecast_days);
                for jj = 1 : NumNPI
                    if(scenario < num_random_input_monte_carlo_runs/2) % random over NPI, constant over time
                        u(jj, :) = randi([NPI_MINS(jj), NPI_MAXES(jj)]);
                    else
                        for t = 1 : num_forecast_days % random over NPI and time
                            u(jj, t) = randi([NPI_MINS(jj), NPI_MAXES(jj)]);
                        end
                    end
                end
                %         u = u + 0.1*rand(size(u)); % add noise to input
                [s_controlled, i_controlled, alpha_controlled] = SIalpha_Controlled(u, s_historic(end), i_historic(end), alpha_historic(end), NPI_MAXES, params.alpha_min, params.alpha_max, params.gamma, params.a, params.b, params.beta, s_noise_std, i_noise_std, alpha_noise_std, num_forecast_days, params.dt);
                
                s = [s_historic, s_controlled];
                i = [i_historic, i_controlled];
                alpha = [alpha_historic, alpha_controlled];
                u = cat(2, IP, u);
                
                [J0(scenario), J1(scenario)] = NPICost(s.*i.*alpha, u, npi_weights_day_wise);
                
            end
            
        end
        
        %         MedFatalityRate = median(FatalityRate);
        %         MedRecentFatalityRate = median(FatalityRate(round(0.75*end) : end));
        
        %         CumInfections = cumsum(N_population * S_SMOOTH_round2_withinput(2, :));
        %         DeathToCumInfectionsRatio = ConfirmedDeaths ./ CumInfections(:);
        %
        %         BetaEstimate = DeathToCumInfectionsRatio/MedRecentFatalityRate;
        %         MedRecentBetaEstimate = median(BetaEstimate(round(0.75*end) : end));
        
        % PLOT RESULTS
        if(plot_results)
            varnames = {'s', 'i', '\alpha'};
            figure
            for mm = 1 : 3
                subplot(3, 1, mm);
                errorbar(S_MINUS_round1_noinput(mm, :), sqrt(squeeze(P_MINUS_round1_noinput(mm, mm, :))));
                hold on
                errorbar(S_PLUS_round1_noinput(mm, :), sqrt(squeeze(P_PLUS_round1_noinput(mm, mm, :))));
                errorbar(S_SMOOTH_round1_noinput(mm, :), sqrt(squeeze(P_SMOOTH_round1_noinput(mm, mm, :))));
                legend('S_MINUS', 'S_PLUS', 'S_SMOOTH');
                grid
                if(mm == 1)
                    title(strjoin([CountryName , ' ', RegionName ' state estimates assuming no input NPI']));
                end
                ylabel(varnames{mm});
            end
            
            figure
            for mm = 1 : 3
                subplot(3, 1, mm);
                errorbar(S_MINUS_round2_withinput(mm, :), sqrt(squeeze(P_MINUS_round2_withinput(mm, mm, :))));
                hold on
                errorbar(S_PLUS_round2_withinput(mm, :), sqrt(squeeze(P_PLUS_round2_withinput(mm, mm, :))));
                errorbar(S_SMOOTH_round2_withinput(mm, :), sqrt(squeeze(P_SMOOTH_round2_withinput(mm, mm, :))));
                legend('S_MINUS', 'S_PLUS', 'S_SMOOTH');
                grid
                if(mm == 1)
                    title(strjoin([CountryName , ' ', RegionName ' state estimates assuming with historic input NPI']));
                end
                ylabel(varnames{mm});
            end
            
            varnames = {'s', 'i', '\alpha', '\lambda_1', '\lambda_2', '\lambda_3'};
            figure
            for mm = 1 : 6
                subplot(6, 1, mm);
                errorbar(S_MINUS_optinput(mm, :), sqrt(squeeze(P_MINUS_optinput(mm, mm, :))));
                hold on
                errorbar(S_PLUS_optinput(mm, :), sqrt(squeeze(P_PLUS_optinput(mm, mm, :))));
                errorbar(S_SMOOTH_optinput(mm, :), sqrt(squeeze(P_SMOOTH_optinput(mm, mm, :))));
                legend('S_MINUS', 'S_PLUS', 'S_SMOOTH');
                grid
                if(mm == 1)
                    title(strjoin([CountryName , ' ', RegionName ' state estimates assuming with optimal input NPI']));
                end
                ylabel(varnames{mm});
            end
            
            if(0)
                figure
                for mm = 1 : 6
                    subplot(6, 1, mm);
                    errorbar(S_PLUS_optinput(mm, :), sqrt(squeeze(P_PLUS_optinput(mm, mm, :))));
                    hold on
                    errorbar(S_SMOOTH_optinput(mm, :), sqrt(squeeze(P_SMOOTH_optinput(mm, mm, :))));
                    errorbar(S_PLUS_backward(mm, :), sqrt(squeeze(P_PLUS_backward(mm, mm, :))));
                    errorbar(S_SMOOTH_backward(mm, :), sqrt(squeeze(P_SMOOTH_backward(mm, mm, :))));
                    %                 errorbar(S_FRW_BCK(mm, :), sqrt(squeeze(P_FRW_BCK(mm, mm, :))));
                    legend('S_PLUS', 'S_SMOOTH', 'S_PLUS_bck', 'S_SMOOTH_bck');%, 'S_FRW_BCK');
                    grid
                    if(mm == 1)
                        title(strjoin([CountryName , ' ', RegionName ' forward vs backwards state estimates']));
                    end
                    ylabel(varnames{mm});
                end
            end
            
            figure
            subplot(211);
            plot(N_population * S_PLUS_zeroinput(1, :).*S_PLUS_zeroinput(2, :).*S_PLUS_zeroinput(3, :));
            hold on
            plot(N_population * S_PLUS_fixedinput(1, :).*S_PLUS_fixedinput(2, :).*S_PLUS_fixedinput(3, :));
            plot(N_population * S_PLUS_optinput(1, :).*S_PLUS_optinput(2, :).*S_PLUS_optinput(3, :));
            grid
            legend('new cases - no NPI', 'new cases - fixed NPI', 'new cases - optimal NPI');
            title(CountryAndRegionList(k), 'interpreter', 'none');
            subplot(212);
            plot(S_PLUS_zeroinput(3, :));
            hold on
            plot(S_PLUS_fixedinput(3, :));
            plot(S_PLUS_optinput(3, :));
            legend('alpha - no NPI', 'alpha - fixed NPI', 'alpha - optimal NPI');
            grid
            
            figure
            subplot(411)
            plot(N_population * NewCasesSmoothedNormalized, 'linewidth', 2);
            hold on
            plot(N_population * S_PLUS_round1_noinput(1, :).*S_PLUS_round1_noinput(2, :).*S_PLUS_round1_noinput(3, :));
            plot(N_population * S_SMOOTH_round1_noinput(1, :).*S_SMOOTH_round1_noinput(2, :).*S_SMOOTH_round1_noinput(3, :));
            plot(N_population * S_PLUS_round2_withinput(1, :).*S_PLUS_round2_withinput(2, :).*S_PLUS_round2_withinput(3, :));
            plot(N_population * S_SMOOTH_round2_withinput(1, :).*S_SMOOTH_round2_withinput(2, :).*S_SMOOTH_round2_withinput(3, :));
            legend('NewCases', 'PLUS', 'SMOOTH', 'PLUS2', 'SMOOTH2');
            title(CountryAndRegionList(k), 'interpreter', 'none');
            grid;
            
            subplot(412);
            plot(rho_round1_noinput)
            hold on
            plot(rho_round2_withinput, 'r')
            legend('rho_round1_noinput', 'rho_round2_withinput');
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
            
            if(0)
                figure
                %         plot(NewCasesSmoothed);
                plot(ConfirmedCases);
                hold on
                plot(ConfirmedDeaths);
                plot(N_population * S_PLUS_round1_noinput(1, :) .* S_PLUS_round1_noinput(2, :) .* S_PLUS_round1_noinput(3, :));
                plot(cumsum(N_population * S_PLUS_round1_noinput(2, :)));
                grid
                legend('ConfirmedCases', 'ConfirmedDeaths', 'Estimated new cases', 'Estimated infections');
                %         legend('NewCasesSmoothed', 'ConfirmedDeaths', 'Estimated new cases', 'Estimated infections');
                
                figure
                subplot(211);
                plot(100.0 * FatalityRate);
                hold on
                plot(100.0 * MedFatalityRate(ones(1, length(FatalityRate))), 'linewidth', 2);
                plot(100.0 * MedRecentFatalityRate(ones(1, length(FatalityRate))), 'linewidth', 2);
                plot(CaseFatalityJHDB(ones(1, length(FatalityRate))), 'linewidth', 2);
                title([CountryAndRegionList(k) ' mortality rate (%)'], 'interpreter', 'none');
                grid
                legend('FatalityRate', 'MedFatalityRate', 'MedRecentFatalityRate', 'CaseFatalityJHDB');
                
                subplot(212);
                plot(BetaEstimate);
                hold on
                plot(MedRecentBetaEstimate(ones(1, length(BetaEstimate))), 'linewidth', 2);
                legend('\beta estimate', 'MedRecentBetaEstimate');
                grid
                
                figure
                hold on
                plot(diff(ConfirmedCases));
                plot(NewCasesSmoothed);
                plot(NewDeathsSmoothed/FatalityRate(end));
                %     plot(ConfirmedCases - ConfirmedDeaths);
                %     legend('ConfirmedCases', 'ConfirmedDeaths/FatalityRate', 'RecoveredCases');
                legend('diff(ConfirmedCases)', 'NewCasesSmoothed', 'NewDeathsSmoothedNormalized');
                title([CountryAndRegionList(k) ' mortality rate (%)'], 'interpreter', 'none');
                grid
            end
            
            figure
            subplot(311);
            errorbar(N_population*S_MINUS_round1_noinput(1, :), N_population*sqrt(squeeze(P_MINUS_round1_noinput(1, 1, :))));
            hold on
            errorbar(N_population*S_PLUS_round1_noinput(1, :), N_population*sqrt(squeeze(P_PLUS_round1_noinput(1, 1, :))));
            errorbar(N_population*S_SMOOTH_round1_noinput(1, :), N_population*sqrt(squeeze(P_SMOOTH_round1_noinput(1, 1, :))));
            errorbar(N_population*S_MINUS_round2_withinput(1, :), N_population*sqrt(squeeze(P_MINUS_round2_withinput(1, 1, :))));
            errorbar(N_population*S_PLUS_round2_withinput(1, :), N_population*sqrt(squeeze(P_PLUS_round2_withinput(1, 1, :))));
            errorbar(N_population*S_SMOOTH_round2_withinput(1, :), N_population*sqrt(squeeze(P_SMOOTH_round2_withinput(1, 1, :))));
            %         errorbar(N_population*S_MINUS_backward(1, :), N_population*sqrt(squeeze(P_MINUS_backward(1, 1, :))));
            %         errorbar(N_population*S_PLUS_backward(1, :), N_population*sqrt(squeeze(P_PLUS_backward(1, 1, :))));
            %         errorbar(N_population*S_SMOOTH_backward(1, :), N_population*sqrt(squeeze(P_SMOOTH_backward(1, 1, :))));
            legend('MINUS', 'PLUS', 'SMOOTH', 'MINUS2', 'PLUS2', 'SMOOTH2');%, 'MINUS_backward', 'PLUS_backward', 'SMOOTH_backward');
            title(CountryAndRegionList(k), 'interpreter', 'none');
            grid
            subplot(312);
            errorbar(N_population*S_MINUS_round1_noinput(2, :), N_population*sqrt(squeeze(P_MINUS_round1_noinput(2, 2, :))));
            hold on
            errorbar(N_population*S_PLUS_round1_noinput(2, :), N_population*sqrt(squeeze(P_PLUS_round1_noinput(2, 2, :))));
            errorbar(N_population*S_SMOOTH_round1_noinput(2, :), N_population*sqrt(squeeze(P_SMOOTH_round1_noinput(2, 2, :))));
            %         plot(diff(ConfirmedDeathsSmoothed)/(MedRecentFatalityRate*params.beta), 'linewidth', 2);
            errorbar(N_population*S_MINUS_round2_withinput(2, :), N_population*sqrt(squeeze(P_MINUS_round2_withinput(2, 2, :))));
            errorbar(N_population*S_PLUS_round2_withinput(2, :), N_population*sqrt(squeeze(P_PLUS_round2_withinput(2, 2, :))));
            errorbar(N_population*S_SMOOTH_round2_withinput(2, :), N_population*sqrt(squeeze(P_SMOOTH_round2_withinput(2, 2, :))));
            %         errorbar(N_population*S_MINUS_backward(2, :), N_population*sqrt(squeeze(P_MINUS_backward(2, 2, :))));
            %         errorbar(N_population*S_PLUS_backward(2, :), N_population*sqrt(squeeze(P_PLUS_backward(2, 2, :))));
            %         errorbar(N_population*S_SMOOTH_backward(2, :), N_population*sqrt(squeeze(P_SMOOTH_backward(2, 2, :))));
            plot(diff(ConfirmedDeathsSmoothed)./(FatalityRate(2:end)*params.beta), 'linewidth', 2);
            legend('MINUS', 'PLUS', 'SMOOTH', 'MINUS2', 'PLUS2', 'SMOOTH2', 'diff(ConfirmedDeathsSmoothed)/(MedRecentFatalityRate*params.beta)');
            title(CountryAndRegionList(k), 'interpreter', 'none');
            grid
            subplot(313);
            errorbar(S_MINUS_round1_noinput(3, :), sqrt(squeeze(P_MINUS_round1_noinput(3, 3, :))));
            hold on
            errorbar(S_PLUS_round1_noinput(3, :), sqrt(squeeze(P_PLUS_round1_noinput(3, 3, :))));
            errorbar(S_SMOOTH_round1_noinput(3, :), sqrt(squeeze(P_SMOOTH_round1_noinput(3, 3, :))));
            errorbar(S_MINUS_round2_withinput(3, :), sqrt(squeeze(P_MINUS_round2_withinput(3, 3, :))));
            errorbar(S_PLUS_round2_withinput(3, :), sqrt(squeeze(P_PLUS_round2_withinput(3, 3, :))));
            errorbar(S_SMOOTH_round2_withinput(3, :), sqrt(squeeze(P_SMOOTH_round2_withinput(3, 3, :))));
            %         errorbar(S_MINUS_backward(3, :), sqrt(squeeze(P_MINUS_backward(3, 3, :))));
            %         errorbar(S_PLUS_backward(3, :), sqrt(squeeze(P_PLUS_backward(3, 3, :))));
            %         errorbar(S_SMOOTH_backward(3, :), sqrt(squeeze(P_SMOOTH_backward(3, 3, :))));
            legend('MINUS', 'PLUS', 'SMOOTH', 'MINUS2', 'PLUS2', 'SMOOTH2');%, 'MINUS_backward', 'PLUS_backward', 'SMOOTH_backward');
            title(CountryAndRegionList(k), 'interpreter', 'none');
            grid
            
            figure
            hold on
            %     plot(log(J0 * N_population), log(J1), 'bo');
            plot(J0 * N_population, J1, 'bo');
            plot(J0_opt_control * N_population, J1_opt_control, 'ro');
            grid
            %     plot(log(J0_zero_control), log(J1_zero_control), 'kx', 'MarkerSize',12);
            %     plot(log(J0_full_control), log(J1_full_control), 'gx', 'MarkerSize',12);
            %     plot(log(J0_opt_control), log(J1_opt_control), 'ro');
            %             plot(J0_zero_control, J1_zero_control, 'kx', 'MarkerSize',12);
            %             plot(J0_full_control, J1_full_control, 'gx', 'MarkerSize',12);
            %     axis square
            xlabel('Human factor');
            ylabel('NPI cost');
            axis tight
            axis square
            set(gca, 'fontsize', 16)
            set(gca, 'box', 'on')
            
            if(0)
                figure
                subplot(311);
                errorbar(S_MINUS_round1_noinput(4, :), sqrt(squeeze(P_MINUS_round1_noinput(4, 4, :))));
                hold on
                errorbar(S_PLUS_round1_noinput(4, :), sqrt(squeeze(P_PLUS_round1_noinput(4, 4, :))));
                errorbar(S_SMOOTH_round1_noinput(4, :), sqrt(squeeze(P_SMOOTH_round1_noinput(4, 4, :))));
                errorbar(S_MINUS_round2_withinput(4, :), sqrt(squeeze(P_MINUS_round2_withinput(1, 1, :))));
                errorbar(S_PLUS_round2_withinput(4, :), sqrt(squeeze(P_PLUS_round2_withinput(4, 4, :))));
                errorbar(S_SMOOTH_round2_withinput(4, :), sqrt(squeeze(P_SMOOTH_round2_withinput(4, 4, :))));
                legend('MINUS', 'PLUS', 'SMOOTH', 'MINUS2', 'PLUS2', 'SMOOTH2');
                title(CountryAndRegionList(k), 'interpreter', 'none');
                grid
                subplot(312);
                errorbar(S_MINUS_round1_noinput(5, :), sqrt(squeeze(P_MINUS_round1_noinput(5, 5, :))));
                hold on
                errorbar(S_PLUS_round1_noinput(5, :), sqrt(squeeze(P_PLUS_round1_noinput(5, 5, :))));
                errorbar(S_SMOOTH_round1_noinput(5, :),sqrt(squeeze(P_SMOOTH_round1_noinput(5, 5, :))));
                errorbar(S_MINUS_round2_withinput(5, :), sqrt(squeeze(P_MINUS_round2_withinput(5, 5, :))));
                errorbar(S_PLUS_round2_withinput(5, :), sqrt(squeeze(P_PLUS_round2_withinput(5, 5, :))));
                errorbar(S_SMOOTH_round2_withinput(5, :), sqrt(squeeze(P_SMOOTH_round2_withinput(5, 5, :))));
                legend('MINUS', 'PLUS', 'SMOOTH', 'MINUS2', 'PLUS2', 'SMOOTH2');
                title(CountryAndRegionList(k), 'interpreter', 'none');
                grid
                subplot(313);
                errorbar(S_MINUS_round1_noinput(6, :), sqrt(squeeze(P_MINUS_round1_noinput(6, 6, :))));
                hold on
                errorbar(S_PLUS_round1_noinput(6, :), sqrt(squeeze(P_PLUS_round1_noinput(6, 6, :))));
                errorbar(S_SMOOTH_round1_noinput(6, :), sqrt(squeeze(P_SMOOTH_round1_noinput(6, 6, :))));
                errorbar(S_MINUS_round2_withinput(6, :), sqrt(squeeze(P_MINUS_round2_withinput(6, 6, :))));
                errorbar(S_PLUS_round2_withinput(6, :), sqrt(squeeze(P_PLUS_round2_withinput(6, 6, :))));
                errorbar(S_SMOOTH_round2_withinput(6, :), sqrt(squeeze(P_SMOOTH_round2_withinput(6, 6, :))));
                legend('MINUS', 'PLUS', 'SMOOTH', 'MINUS2', 'PLUS2', 'SMOOTH2');
                title(CountryAndRegionList(k), 'interpreter', 'none');
                grid
            end
            
            %             pause(2)
            %             input('Press enter/return to process next country/region.');
            close all
        end
        
        TrainedModelParams = cat(1, TrainedModelParams, [CountryName, RegionName, {N_population}, {reg_coef_b}, {reg_coef_a}, {reg_coef_b2}, {reg_coef_a2}]);
    end
end
save(char(trained_model_params_file), 'TrainedModelParams');

