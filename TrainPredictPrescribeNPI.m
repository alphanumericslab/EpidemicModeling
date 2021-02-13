function TrainPredictPrescribeNPI(npi_weights, human_npi_cost_factor, start_train_date_str, end_train_date_str, end_predict_presscribe_date_str, data_file, geo_file, populations_file, included_IP, NPI_MINS, NPI_MAXES, trained_model_params_file)

% npi_weights: NPI weights; not used during training; only used during prescription phase.

% parameters
plot_results = true; % plot per-region/country plots or not
SmoothingWinLen = 7; % window length used for smoothing the new cases
min_cases = 10; % the absolute minimum number of cases at start date for processing each region/country
first_num_days_for_case_estimation = 7; % the first few days used to find an initial estimate of the first cases (only for EKF initialization)
model_gamma_param = 7; % input to contact influence time constant
observation_type = 'NEWCASES'; % TOTALCASES or NEWCASES
num_days_for_beta_calculation = 21;
prob_contagion_after_Tdays = 0.01;
R0 = 2.5; % An assumption during outbreak
REGRESSION_TYPE = 'NONNEGATIVELS'; % 'LASSO' or 'NONNEGATIVELS'

% Convert training start date string to number
start_train_date_chars = char(start_train_date_str);
dash_indexes = start_train_date_chars == '-';
start_train_date_chars(dash_indexes) = [];
start_train_date = str2double(string(start_train_date_chars));

% Convert training end date string to number
end_train_date_chars = char(end_train_date_str);
dash_indexes = end_train_date_chars == '-';
end_train_date_chars(dash_indexes) = [];
end_train_date = str2double(string(end_train_date_chars));

% Convert end of prediction/prescription date string to number
end_predict_presscribe_date_chars = char(end_predict_presscribe_date_str);
dash_indexes = end_predict_presscribe_date_chars == '-';
end_predict_presscribe_date_chars(dash_indexes) = [];
end_predict_presscribe_date = str2double(string(end_predict_presscribe_date_chars));

% Assess: start_train_date <= end_train_date <= end_predict_presscribe_date
if(~(start_train_date <= end_train_date && end_train_date <= end_predict_presscribe_date))
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

TrainedModelParams = [{'CountryName'}, {'RegionName'}, {'N_population'}, {'coef0'}, {'coef'}, {'coef0_2'}, {'coef_2'}];

for k = 1 : NumGeoLocations
    if(find(SelectedGeoIDs == CountryAndRegionList(k))) % make sure the current country/region is among the ones to be processed
        % 1) READ AND CLEAN THE DATA
        disp([num2str(k) '- ' char(CountryAndRegionList(k))]);
        geoid_all_row_indexes = AllGeoIDs == CountryAndRegionList(k) & all_data.Date >= start_train_date & all_data.Date <= end_train_date; % fetch all rows corresponding to the country/region from start date to end date
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
        params.alpha_min = 0.0; % minimum alpha
        params.alpha_max = inf; % maximum alpha
        params.epsilon = nan; % No NPI costs during the first round
        params.gamma = 1/(params.dt * model_gamma_param); % input to contact influence rate (inverse time)
        params.obs_type = observation_type;
        
        % see the following for the background on beta: https://www.cdc.gov/coronavirus/2019-ncov/hcp/duration-isolation.html
        Tdays = num_days_for_beta_calculation * params.dt;
        params.beta = -log(prob_contagion_after_Tdays)/Tdays; % recovery rate from being contagious (inverse time)
        alpha0 = params.beta + log(R0)/params.dt; % The logic is that in the SIalpha model, during the outbreak R0 = exp(dt*(alpha - beta)) and alpha = beta is the metastable threshold (R0 = 1) %1.0/N_population; % the per-population normalization is needed
        
        params.sigma = 1000000; % sigmoid function slope
        beta_ekf = .9; % Observation noise update factor (set to 1 for no update)
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
        % Ps_init = 100.0*(params.dt)^2*diag([I0/N_population, I0/N_population, alpha_noise_std].^2); % Covariance matrix of initial states
        Ps_init = 100.0*(params.dt)^2*diag([s_noise_std, i_noise_std, alpha_noise_std].^2); % Covariance matrix of initial states
        s_final = nan(3, 1); % Set the finite horizon end points (not required during first round)
        Ps_final = nan(3); % Set the finite horizon end points (not required during first round)
        R_v = var((NewCasesSmoothedZeroLag - NewCasesRefined)/N_population); % Observation noise variance
        if(isequal(params.obs_type, 'TOTALCASES'))
            observations = ConfirmedCasesSmoothedNormalized(:)';
        elseif(isequal(params.obs_type, 'NEWCASES'))
            observations = NewCasesSmoothedNormalized(:)';
        end
        [~, ~, S_MINUS_round1_noinput, S_PLUS_round1_noinput, S_SMOOTH_round1_noinput, P_MINUS, P_PLUS_round1_noinput, P_SMOOTH_round1_noinput, ~, ~, rho_round1_noinput] = SIAlphaModelEKF(control_input, observations, params, s_init, Ps_init, s_final, Ps_final, w_bar, v_bar, Q_w, R_v, beta_ekf, gamma_ekf, inv_monitor_len_ekf, order);
        
        % 3) APPLY LASSO OVER ALPHA ESTIMATE TO FIND INITIAL LINEAR REGRESSION MODEL PARAMS
        x_data_train = NPI_MAXES(:, ones(NumNPIdays, 1))' - InterventionPlans;
        y_data_train = S_SMOOTH_round1_noinput(3, :)'; % alpha
        % y_data_train = S_PLUS_round1_noinput(3, :)'; % alpha
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
        
        % 4) APPLY THE SECOND ROUND EKF TO REFINE ALPHA BASED ON REAL HISTORIC NPI INPUTS
        params.a = coef; % input influence weight vector
        params.b = coef0; % input influence bias constant
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
        
        % 6) FORECAST/PRESCRIPTION PHASE
        num_forecast_days = end_predict_presscribe_date - end_train_date;
        if(num_forecast_days > 0)
            params.a = coef_2; % input influence weight vector
            params.b = coef0_2; % input influence bias constant
            
            % A) Forecast/prescription phase with zero input
            IP = InterventionPlans';
            control_input = [IP, NPI_MINS(:, ones(1, num_forecast_days))]; % The second round works with real inputs
            %s_init = [(N_population - I0)/N_population ; I0/N_population ; S_SMOOTH_round1_noinput(3, 1) ; S_SMOOTH_round1_noinput(4, 1) ; S_SMOOTH_round1_noinput(5, 1) ; S_SMOOTH_round1_noinput(6, 1)]; % initial state vector
            if(isequal(params.obs_type, 'TOTALCASES'))
                observations = [ConfirmedCasesSmoothedNormalized(:)', nan(1, num_forecast_days)];
            elseif(isequal(params.obs_type, 'NEWCASES'))
                observations = [NewCasesSmoothedNormalized(:)', nan(1, num_forecast_days)];
            end
            [zero_control_input, ~, S_MINUS_zeroinput, S_PLUS_zeroinput, S_SMOOTH_zeroinput, P_MINUS_zeroinput, P_PLUS_zeroinput, P_SMOOTH_zeroinput, ~, ~, rho_zeroinput] = SIAlphaModelEKF(control_input, observations, params, s_init, Ps_init, s_final, Ps_final, w_bar, v_bar, Q_w, R_v, beta_ekf, gamma_ekf, inv_monitor_len_ekf, 1);
            
            % B) Forecast/prescription phase with fixed input
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
            
            % C) Forecast/prescription phase with optimal input
            num_pareto_front_points = length(human_npi_cost_factor);
            J0_opt_control = zeros(1, num_pareto_front_points);
            J1_opt_control = zeros(1, num_pareto_front_points);
            for ll = 1 : num_pareto_front_points
                disp(['Pareto front #' num2str(ll)]);
                params.w = npi_weights;
                params.epsilon = human_npi_cost_factor(ll); % [0, 1]: 0 neglects NPI cost and 1 neglects human factor!
                lambda0 = 1;
                q_lambda = 10;
                %     s_init = [s_historic(end) ; i_historic(end) ; alpha0 ; lambda0 ; lambda0 ; lambda0]; % initial state vector
                ss_init = cat(1, s_init , [lambda0 ; lambda0 ; lambda0]); % initial state vector
                QQ_w = blkdiag(Q_w, 10.0*(params.dt)^2*diag([q_lambda, q_lambda, q_lambda].^2)); % Covariance matrix of initial states
                PPs_init = blkdiag(Ps_init, 100.0*(params.dt)^2*diag([q_lambda, q_lambda, q_lambda].^2)); % Covariance matrix of initial states
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
                control_input = [InterventionPlans', nan(NumNPI, num_forecast_days)]; % The second round works with real inputs
                %s_init = [(N_population - I0)/N_population ; I0/N_population ; S_SMOOTH_round1_noinput(3, 1) ; S_SMOOTH_round1_noinput(4, 1) ; S_SMOOTH_round1_noinput(5, 1) ; S_SMOOTH_round1_noinput(6, 1)]; % initial state vector
                if(isequal(params.obs_type, 'TOTALCASES'))
                    observations = [ConfirmedCasesSmoothedNormalized(:)', nan(1, num_forecast_days)];
                elseif(isequal(params.obs_type, 'NEWCASES'))
                    observations = [NewCasesSmoothedNormalized(:)', nan(1, num_forecast_days)];
                end
                [opt_control_input, opt_control_input_smooth, S_MINUS_optinput, S_PLUS_optinput, S_SMOOTH_optinput, P_MINUS_optinput, P_PLUS_optinput, P_SMOOTH_optinput, ~, ~, rho_optinput] = SIAlphaModelEKFOptControlled(control_input, observations, params, ss_init, PPs_init, s_final, Ps_final, w_bar, v_bar, QQ_w, R_v, beta_ekf, gamma_ekf, inv_monitor_len_ekf, 1);
                % % %             [~, S_MINUS_backward, S_PLUS_backward, S_SMOOTH_backward, P_MINUS_backward, P_PLUS_backward, P_SMOOTH_backward, ~, ~, rho_backward] = SIAlphaModelBackwardEKF(control_input, NewCasesSmoothedNormalized(:)', params, s_init, Ps_init, s_final, Ps_final, w_bar, v_bar, Q_w, R_v, beta_ekf, gamma_ekf, inv_monitor_len_ekf, order);
                
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
        end
        
        %         MedFatalityRate = median(FatalityRate);
        %         MedRecentFatalityRate = median(FatalityRate(round(0.75*end) : end));
        
        %         CumInfections = cumsum(N_population * S_SMOOTH_round2_withinput(2, :));
        %         DeathToCumInfectionsRatio = ConfirmedDeaths ./ CumInfections(:);
        %
        %         BetaEstimate = DeathToCumInfectionsRatio/MedRecentFatalityRate;
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
        
        if(plot_results)
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
            % % %         plot(AlphaHatARX)
            % % %         hold on
            % % %         plot(ar_train_segment, 'r')
            % % %         grid
            
            % % % %         figure
            % % % % %         plot(NewCasesSmoothed);
            % % % %         plot(ConfirmedCases);
            % % % %         hold on
            % % % %         plot(ConfirmedDeaths);
            % % % %         plot(N_population * S_PLUS_round1_noinput(1, :) .* S_PLUS_round1_noinput(2, :) .* S_PLUS_round1_noinput(3, :));
            % % % %         plot(cumsum(N_population * S_PLUS_round1_noinput(2, :)));
            % % % %         grid
            % % % %         legend('ConfirmedCases', 'ConfirmedDeaths', 'Estimated new cases', 'Estimated infections');
            % % % % %         legend('NewCasesSmoothed', 'ConfirmedDeaths', 'Estimated new cases', 'Estimated infections');
            
            % % %         figure
            % % %         subplot(211);
            % % %         plot(100.0 * FatalityRate);
            % % %         hold on
            % % %         plot(100.0 * MedFatalityRate(ones(1, length(FatalityRate))), 'linewidth', 2);
            % % %         plot(100.0 * MedRecentFatalityRate(ones(1, length(FatalityRate))), 'linewidth', 2);
            % % %         plot(CaseFatalityJHDB(ones(1, length(FatalityRate))), 'linewidth', 2);
            % % %         title([CountryAndRegionList(k) ' mortality rate (%)'], 'interpreter', 'none');
            % % %         grid
            % % %         legend('FatalityRate', 'MedFatalityRate', 'MedRecentFatalityRate', 'CaseFatalityJHDB');
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
            % % %         plot(NewDeathsSmoothed/FatalityRate(end));
            % % %         %     plot(ConfirmedCases - ConfirmedDeaths);
            % % %         %     legend('ConfirmedCases', 'ConfirmedDeaths/FatalityRate', 'RecoveredCases');
            % % %         legend('diff(ConfirmedCases)', 'NewCasesSmoothed', 'NewDeathsSmoothedNormalized');
            % % %         title([CountryAndRegionList(k) ' mortality rate (%)'], 'interpreter', 'none');
            % % %         grid
            
            
            figure
            subplot(311);
            errorbar(N_population*S_MINUS_round1_noinput(1, :), N_population*sqrt(squeeze(P_MINUS(1, 1, :))));
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
            errorbar(N_population*S_MINUS_round1_noinput(2, :), N_population*sqrt(squeeze(P_MINUS(2, 2, :))));
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
            errorbar(S_MINUS_round1_noinput(3, :), sqrt(squeeze(P_MINUS(3, 3, :))));
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
            %     plot(log(J0), log(J1), 'bo');
            plot(J0_opt_control, J1_opt_control, 'ro');
            hold on
            %             plot(J0, J1, 'bo');
            grid
            %     plot(log(J0_zero_control), log(J1_zero_control), 'kx', 'MarkerSize',12);
            %     plot(log(J0_full_control), log(J1_full_control), 'gx', 'MarkerSize',12);
            %     plot(log(J0_opt_control), log(J1_opt_control), 'ro');
            %             plot(J0_zero_control, J1_zero_control, 'kx', 'MarkerSize',12);
            %             plot(J0_full_control, J1_full_control, 'gx', 'MarkerSize',12);
            %     axis square
            xlabel('Human factor');
            ylabel('NPI cost');
            axis square
            
            % % %         figure
            % % %         subplot(311);
            % % %         errorbar(S_MINUS_round1_noinput(4, :), sqrt(squeeze(P_MINUS(4, 4, :))));
            % % %         hold on
            % % %         errorbar(S_PLUS_round1_noinput(4, :), sqrt(squeeze(P_PLUS_round1_noinput(4, 4, :))));
            % % %         errorbar(S_SMOOTH_round1_noinput(4, :), sqrt(squeeze(P_SMOOTH_round1_noinput(4, 4, :))));
            % % %         errorbar(S_MINUS_round2_withinput(4, :), sqrt(squeeze(P_MINUS_round2_withinput(1, 1, :))));
            % % %         errorbar(S_PLUS_round2_withinput(4, :), sqrt(squeeze(P_PLUS_round2_withinput(4, 4, :))));
            % % %         errorbar(S_SMOOTH_round2_withinput(4, :), sqrt(squeeze(P_SMOOTH_round2_withinput(4, 4, :))));
            % % %         legend('MINUS', 'PLUS', 'SMOOTH', 'MINUS2', 'PLUS2', 'SMOOTH2');
            % % %         title(CountryAndRegionList(k), 'interpreter', 'none');
            % % %         grid
            % % %         subplot(312);
            % % %         errorbar(S_MINUS_round1_noinput(5, :), sqrt(squeeze(P_MINUS(5, 5, :))));
            % % %         hold on
            % % %         errorbar(S_PLUS_round1_noinput(5, :), sqrt(squeeze(P_PLUS_round1_noinput(5, 5, :))));
            % % %         errorbar(S_SMOOTH_round1_noinput(5, :),sqrt(squeeze(P_SMOOTH_round1_noinput(5, 5, :))));
            % % %         errorbar(S_MINUS_round2_withinput(5, :), sqrt(squeeze(P_MINUS_round2_withinput(5, 5, :))));
            % % %         errorbar(S_PLUS_round2_withinput(5, :), sqrt(squeeze(P_PLUS_round2_withinput(5, 5, :))));
            % % %         errorbar(S_SMOOTH_round2_withinput(5, :), sqrt(squeeze(P_SMOOTH_round2_withinput(5, 5, :))));
            % % %         legend('MINUS', 'PLUS', 'SMOOTH', 'MINUS2', 'PLUS2', 'SMOOTH2');
            % % %         title(CountryAndRegionList(k), 'interpreter', 'none');
            % % %         grid
            % % %         subplot(313);
            % % %         errorbar(S_MINUS_round1_noinput(6, :), sqrt(squeeze(P_MINUS(6, 6, :))));
            % % %         hold on
            % % %         errorbar(S_PLUS_round1_noinput(6, :), sqrt(squeeze(P_PLUS_round1_noinput(6, 6, :))));
            % % %         errorbar(S_SMOOTH_round1_noinput(6, :), sqrt(squeeze(P_SMOOTH_round1_noinput(6, 6, :))));
            % % %         errorbar(S_MINUS_round2_withinput(6, :), sqrt(squeeze(P_MINUS_round2_withinput(6, 6, :))));
            % % %         errorbar(S_PLUS_round2_withinput(6, :), sqrt(squeeze(P_PLUS_round2_withinput(6, 6, :))));
            % % %         errorbar(S_SMOOTH_round2_withinput(6, :), sqrt(squeeze(P_SMOOTH_round2_withinput(6, 6, :))));
            % % %         legend('MINUS', 'PLUS', 'SMOOTH', 'MINUS2', 'PLUS2', 'SMOOTH2');
            % % %         title(CountryAndRegionList(k), 'interpreter', 'none');
            % % %         grid
            
            % % %             params.beta
            % % %             N_population
            
            
            %             pause(2)
            input('Press enter/return to process next country/region.');
            close all
        end
        
        TrainedModelParams = cat(1, TrainedModelParams, [CountryName, RegionName, {N_population}, {coef0}, {coef}, {coef0_2}, {coef_2}]);
    end
end
save(char(trained_model_params_file), 'TrainedModelParams');

