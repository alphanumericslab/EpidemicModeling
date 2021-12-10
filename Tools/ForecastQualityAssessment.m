function ForecastQualityAssessment(npi_weights, human_npi_cost_factor, start_train_date_str, end_train_date_str, start_regression_date_str, end_predict_presscribe_date_str, MaxLookAheadDays, data_file, geo_file, populations_file, included_IP, NPI_MINS, NPI_MAXES, trained_model_params_file)

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
    if(~isempty(find(SelectedGeoIDs == CountryAndRegionList(k), 1)) && isequal(CountryAndRegionList(k), "United States ")) % make sure the current country/region is among the ones to be processed
        disp([num2str(k) '- ' char(CountryAndRegionList(k))]);
        
        % Get coutry/region population
        population_index = find(strcat(all_populations_data.CountryName, " ", all_populations_data.RegionName) == CountryAndRegionList(k));
        N_population = all_populations_data.Population2020(population_index);
        
        % 0) READ UP TO THE END OF THE DATA (USED ONLY FOR VALIDATION, NOT PREDICTION)
        geoid_all_row_indexes_ENTIRE = AllGeoIDs == CountryAndRegionList(k) & all_data.Date >= start_train_date_number & all_data.Date <= end_predict_presscribe_date_number; % fetch all rows corresponding to the country/region from start date to end date
        ConfirmedCases_ENTIRE = all_data.ConfirmedCases(geoid_all_row_indexes_ENTIRE); % Confirmed cases
        %ConfirmedDeaths_ENTIRE = all_data.ConfirmedDeaths(geoid_all_row_indexes_ENTIRE); % Confirmed deaths
        NewCases_ENTIRE = diff([ConfirmedCases_ENTIRE(1) ; ConfirmedCases_ENTIRE]); % calculate the new daily cases
        NewCases_ENTIRE(NewCases_ENTIRE < 0) = 0;
        if(length(NewCases_ENTIRE) < 2)
            warning(strcat("Insufficient data for ", CountryAndRegionList(k), ". Skipping the country/region."));
            continue;
        end
        % Replace the last nan with the last valid number
        NewCasesRefined_ENTIRE = NewCases_ENTIRE;
        if(isnan(NewCases_ENTIRE(end))) % Just fill the last missing date (if any) with the previous day
            last_number = find(not(isnan(NewCases_ENTIRE)), 1, 'last');
            NewCasesRefined_ENTIRE(end) = NewCasesRefined_ENTIRE(last_number);
        end
        NewCasesRefined_ENTIRE(isnan(NewCasesRefined_ENTIRE)) = 0; % replace any remaining nans with 0
        NewCasesSmoothed_ENTIRE = filter(ones(1, SmoothingWinLen), SmoothingWinLen, NewCasesRefined_ENTIRE); % causal smoothed; used for processing
        ConfirmedCasesSmoothed_ENTIRE = cumsum(NewCasesSmoothed_ENTIRE);
        NewCasesSmoothedNormalized_ENTIRE = NewCasesSmoothed_ENTIRE / N_population; % the fraction of new cases
        ConfirmedCasesSmoothedNormalized_ENTIRE = ConfirmedCasesSmoothed_ENTIRE / N_population; % the fraction of total cases
        
        % Intervention plans
        InterventionPlans_ENTIRE = all_data{geoid_all_row_indexes_ENTIRE, included_IP}; % Region/Country intervention plans
        NumNPIdays_ENTIRE = size(InterventionPlans_ENTIRE, 1); % Number of days with NPI
        % Preprocess the N/A IPs
        for j = 1 : NumNPI % Replace all N/A IP with previous IP
            for i = 2 : NumNPIdays_ENTIRE
                if(isnan(InterventionPlans_ENTIRE(i, j)) && not(isnan(InterventionPlans_ENTIRE(i - 1, j))))
                    InterventionPlans_ENTIRE(i, j) = InterventionPlans_ENTIRE(i - 1, j);
                end
            end
        end
        InterventionPlans_ENTIRE(isnan(InterventionPlans_ENTIRE)) = 0; % Replace any remaining N/A IP with no IP
        
        % 1) READ AND CLEAN THE DATA
        geoid_all_row_indexes = AllGeoIDs == CountryAndRegionList(k) & all_data.Date >= start_train_date_number & all_data.Date <= end_train_date_number; % fetch all rows corresponding to the country/region from start date to end date
        CountryName = all_data.CountryName(find(geoid_all_row_indexes == 1, 1, 'first')); % The current country name
        RegionName = all_data.RegionName(find(geoid_all_row_indexes == 1, 1, 'first')); % The current region name
        ConfirmedCases = all_data.ConfirmedCases(geoid_all_row_indexes); % Confirmed cases
        ConfirmedDeaths = all_data.ConfirmedDeaths(geoid_all_row_indexes); % Confirmed deaths
        
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
        %         ConfirmedCasesSmoothedZeroLag = cumsum(NewCasesSmoothedZeroLag);
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
            R_v = cat(2, R_v, ones(1, num_forecast_days) * mean(R_v)); % Fill in the R_v matrix with the average variance
            if(isequal(params.obs_type, 'TOTALCASES'))
                observations = [ConfirmedCasesSmoothedNormalized(:)', nan(1, num_forecast_days)];
            elseif(isequal(params.obs_type, 'NEWCASES'))
                observations = [NewCasesSmoothedNormalized(:)', nan(1, num_forecast_days)];
            end
            
            % Z) Forecast with actual NPIs from the ground-truth data (used only to evaluate the estimator)
            control_input_ENTIRE = InterventionPlans_ENTIRE'; % The second round works with real inputs
            %s_init = [(N_population - I0)/N_population ; I0/N_population ; S_SMOOTH_round1_noinput(3, 1) ; S_SMOOTH_round1_noinput(4, 1) ; S_SMOOTH_round1_noinput(5, 1) ; S_SMOOTH_round1_noinput(6, 1)]; % initial state vector
            [~, actual_control_input, ~, S_PLUS_actualinput, S_SMOOTH_actualinput, ~, ~, ~, ~, ~, ~] = SIAlphaModelEKF(control_input_ENTIRE, observations, params, s_init, Ps_init, s_final, Ps_final, w_bar, v_bar, Q_w, R_v, beta_ekf, gamma_ekf, inv_monitor_len_ekf, 1);
            
            % Y) Forecast look ahead quality evaluation
            if(isequal(params.obs_type, 'TOTALCASES'))
                observations_ENTIRE = ConfirmedCasesSmoothedNormalized_ENTIRE(:)';
            elseif(isequal(params.obs_type, 'NEWCASES'))
                observations_ENTIRE = NewCasesSmoothedNormalized_ENTIRE(:)';
            end
            EstError_PLUS = zeros(num_forecast_days, MaxLookAheadDays);
            EstError_SMOOTH = zeros(num_forecast_days, MaxLookAheadDays);
            LL = length(observations_ENTIRE);
            for start = 1 : num_forecast_days
                observations_PARTIAL = observations_ENTIRE;
                observations_PARTIAL(LL - start + 1 : LL) = nan;
                [~, ~, ~, S_PLUS_partial, S_SMOOTH_partial, ~, ~, ~, ~, ~, ~] = SIAlphaModelEKF(control_input_ENTIRE, observations_PARTIAL, params, s_init, Ps_init, s_final, Ps_final, w_bar, v_bar, Q_w, R_v, beta_ekf, gamma_ekf, inv_monitor_len_ekf, 1);
                NewCasesSmoothed_EST_PLUS = N_population * S_PLUS_partial(1, :).*S_PLUS_partial(2, :).*S_PLUS_partial(3, :);
                NewCasesSmoothed_EST_SMOOTH = N_population * S_SMOOTH_partial(1, :).*S_SMOOTH_partial(2, :).*S_SMOOTH_partial(3, :);
                error_PLUS = 100 * abs(NewCasesSmoothed_ENTIRE' - NewCasesSmoothed_EST_PLUS) ./ NewCasesSmoothed_ENTIRE';
                error_SMOOTH = 100 * abs(NewCasesSmoothed_ENTIRE' - NewCasesSmoothed_EST_SMOOTH) ./ NewCasesSmoothed_ENTIRE';
                
                last_index = min(LL, LL - start + MaxLookAheadDays);
                EstError_PLUS(start, 1 : last_index - LL + start) = error_PLUS(LL - start + 1 : last_index);
                EstError_SMOOTH(start, 1 : last_index - LL + start) = error_SMOOTH(LL - start + 1 : last_index);
                
                start
                if 0
                    figure
                    hold on
                    plot(0 : length(NewCasesSmoothed_ENTIRE) - 1, NewCasesSmoothed_ENTIRE, 'linewidth', 4);
                    plot(0 : length(S_PLUS_actualinput(1, :)) - 1, N_population * S_PLUS_actualinput(1, :).*S_PLUS_actualinput(2, :).*S_PLUS_actualinput(3, :), 'linewidth', 2);
                    plot(0 : length(S_PLUS_partial(1, :)) - 1, NewCasesSmoothed_EST_PLUS, 'linewidth', 2);
                    plot(0 : length(S_SMOOTH_actualinput(1, :)) - 1, N_population * S_SMOOTH_actualinput(1, :).*S_SMOOTH_actualinput(2, :).*S_SMOOTH_actualinput(3, :), 'linewidth', 2);
                    plot(0 : length(S_SMOOTH_partial(1, :)) - 1, NewCasesSmoothed_EST_SMOOTH, 'linewidth', 2);
                    grid
                    legend('Ground truth', 'Ground truth estimate EKF', 'Ground truth estimate EKF (partial)', 'Ground truth estimate EKS', 'Ground truth estimate EKS (partial)');
                    title(strcat("Number of new cases for ", CountryAndRegionList(k)), 'interpreter', 'none')
                    xlabel(strcat("Days since ", start_train_date_str));
                    ylabel('Number of new cases')
                    set(gca, 'fontsize', 14)
                    %                 set(gca, 'YScale', 'log')
                    aa = axis;
                    aa(3) = 1;
                    axis(aa);
                end
            end
%             EstError_PLUS = EstError_PLUS(:, 1 : MaxLookAheadDays);
%             EstError_SMOOTH = EstError_SMOOTH(:, 1 : MaxLookAheadDays);
            
            figure
            hold on
            %subplot(211);
            %plot(1 : MaxLookAheadDays, EstError_PLUS', 'color', 0.6 * ones(1, 3))
            %hold on
            %plot(1 : MaxLookAheadDays, median(EstError_PLUS(MaxLookAheadDays:end, :), 1), 'k', 'linewidth', 3);
            %grid
            %subplot(212);
            plot(1 : MaxLookAheadDays, mean(EstError_SMOOTH(MaxLookAheadDays:end, :), 1), 'b', 'linewidth', 3);
            plot(1 : MaxLookAheadDays, median(EstError_SMOOTH(MaxLookAheadDays:end, :), 1), 'k', 'linewidth', 3);
            plot(1 : MaxLookAheadDays, EstError_SMOOTH', 'color', 0.6 * ones(1, 3));
            legend('Mean error', 'Median error', 'Error trace per start date');
            %Returns handles to the patch and line objects
            chi = get(gca, 'Children');
            %Reverse the stacking order so that the patch overlays the line
            set(gca, 'Children',flipud(chi))
            xlabel('Look ahead forecasting days')
            ylabel('New cases forecasting error (%)')
            set(gca, 'fontsize', 18)
            set(gca, 'box', 'on')
            grid
            aa = axis;
            aa(1) = 1;
            aa(2) = MaxLookAheadDays;
            axis(aa);

            figure
            hold on
            plot(1 : MaxLookAheadDays, mean(EstError_SMOOTH(MaxLookAheadDays:end, :), 1), 'b', 'linewidth', 3);
            errorbar(1 : MaxLookAheadDays, mean(EstError_SMOOTH(MaxLookAheadDays:end, :), 1), std(EstError_SMOOTH(MaxLookAheadDays:end, :), [], 1)/2, 'r');
            plot(1 : MaxLookAheadDays, EstError_SMOOTH', 'color', 0.6 * ones(1, 3));
            legend('Mean', '\pm STD', 'Error trace per start date');
            %Returns handles to the patch and line objects
            chi = get(gca, 'Children');
            %Reverse the stacking order so that the patch overlays the line
            set(gca, 'Children',flipud(chi))
            xlabel('Look ahead forecasting days')
            ylabel('New cases forecasting error (%)')
            set(gca, 'fontsize', 18)
            set(gca, 'box', 'on')
            grid
            aa = axis;
            aa(1) = 1;
            aa(2) = MaxLookAheadDays;
            axis(aa);
        end
        
        % PLOT RESULTS
        if(plot_results)
            if 1
                figure
                hold on
                plot(0 : length(NewCasesSmoothed_ENTIRE) - 1, NewCasesSmoothed_ENTIRE, 'linewidth', 2);
                plot(0 : length(S_PLUS_actualinput(1, :)) - 1, N_population * S_PLUS_actualinput(1, :).*S_PLUS_actualinput(2, :).*S_PLUS_actualinput(3, :), 'linewidth', 2);
                plot(0 : length(S_SMOOTH_actualinput(1, :)) - 1, N_population * S_SMOOTH_actualinput(1, :).*S_SMOOTH_actualinput(2, :).*S_SMOOTH_actualinput(3, :), 'linewidth', 2);
                grid
                legend('Ground truth', 'Ground truth estimate EKF', 'Ground truth estimate EKS');
                title(strcat("Number of new cases for ", CountryAndRegionList(k)), 'interpreter', 'none')
                xlabel(strcat("Days since ", start_train_date_str));
                ylabel('Number of new cases')
                set(gca, 'fontsize', 14)
                %                 set(gca, 'YScale', 'log')
                aa = axis;
                aa(3) = 1;
                axis(aa);
                %                 line([num_regression_days, num_regression_days], aa(3), aa(4), 'Color', 'k', 'LineStyle','--', 'linewidth', 3)
                
                figure
                hold on
                plot(0 : length(ConfirmedCasesSmoothed_ENTIRE) - 1, ConfirmedCasesSmoothed_ENTIRE, 'linewidth', 2);
                plot(0 : length(S_PLUS_actualinput(1, :)) - 1, N_population * (1.0 - S_PLUS_actualinput(1, :)), 'linewidth', 2);
                plot(0 : length(S_SMOOTH_actualinput(1, :)) - 1, N_population * (1.0 - S_SMOOTH_actualinput(1, :)), 'linewidth', 2);
                grid
                legend('Ground truth', 'Ground truth estimate EKF', 'Ground truth estimate EKS');
                title(strcat("Number of total cases for different NPI scenarios for ", CountryAndRegionList(k)), 'interpreter', 'none')
                xlabel(strcat("Days since ", start_train_date_str));
                ylabel('Number of total cases')
                set(gca, 'fontsize', 14)
                %                 set(gca, 'YScale', 'log')
                aa = axis;
                aa(3) = 1;
                axis(aa);
            end
            
            close all
        end
        
        TrainedModelParams = cat(1, TrainedModelParams, [CountryName, RegionName, {N_population}, {reg_coef_b}, {reg_coef_a}, {reg_coef_b2}, {reg_coef_a2}]);
    end
end
save(char(trained_model_params_file), 'TrainedModelParams');

