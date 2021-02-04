function TrainNPIPrescriptor(start_date_str, end_date_str, data_file, geo_file, populations_file, included_IP, NPI_MAXES, trained_model_params_file)

% Convert start date to number
start_date_chars = char(start_date_str);
dash_indexes = start_date_chars == '-';
start_date_chars(dash_indexes) = [];
start_date = str2double(string(start_date_chars));

% Convert end date to number
end_date_chars = char(end_date_str);
dash_indexes = end_date_chars == '-';
end_date_chars(dash_indexes) = [];
end_date = str2double(string(end_date_chars));

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
SelectedGeoIDs = strcat(string(geolist_to_study.CountryName), string(geolist_to_study.RegionName));

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
AllGeoIDs = strcat(string(AllCountryNames), string(AllRegionNames));

% Make a dictionary of country-regions
[CountryAndRegionList, ~, ~] = unique(AllGeoIDs, 'stable');

NumGeoLocations = length(CountryAndRegionList); % Number of country-region pairs

NumNPI = length(included_IP);

% NPI weights (not used during training. just needed as inputs to the function)
human_npi_cost_factor = 0.5; % (not used during training. just needed as inputs to the function)
for k = 2 : NumGeoLocations%index_ : index_%240 %122 : 125%225 %1 : NumGeoLocations
    if(find(SelectedGeoIDs == CountryAndRegionList(k)))
        geoid_all_row_indexes = AllGeoIDs == CountryAndRegionList(k) & all_data.Date >= start_date & all_data.Date <= end_date;
        ConfirmedCases = all_data.ConfirmedCases(geoid_all_row_indexes); % Confirmed cases
        ConfirmedDeaths = all_data.ConfirmedDeaths(geoid_all_row_indexes); % Confirmed deaths
        population_index = find(all_populations_data.CountryName == CountryAndRegionList(k));
        N_population = all_populations_data.Population2020(population_index);
        CaseFatalityJHDB = all_populations_data.CaseFatalityJHDBFeb2021(population_index);
        
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
        first_case_indexes = find(NewCasesSmoothed > 0, 7, 'first'); % find the first non-zero indexes
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
        params.gamma = 1/(params.dt*50.0); % input to contact influence rate (inverse time)
        params.obs_type = 'NEWCASES'; % TOTALCASES or NEWCASES
        
        % see the following for the background on beta: https://www.cdc.gov/coronavirus/2019-ncov/hcp/duration-isolation.html
        prob_contagion_after_Tdays = 0.05;
        Tdays = 21*params.dt;
        params.beta = -log(prob_contagion_after_Tdays)/Tdays; % recovery rate from being contagious (inverse time)
        
        params.sigma = 100000; % sigmoid function slope
        beta = .9;%0.9; % Observation noise update factor (set to 1 for no update)
        gamma = 0.995; % Kalman gain stability factor (set very close to 1, or equal to 1 to disable the feature)
        inv_monitor_len = 21; % Window length for innovations process whiteness monitoring
        order = 1; % 1 for standard EKF; 2 for second-order EKF
        alpha0 = 1.0/N_population; % the per-population normalization is needed
        q_alpha = 2e-3;
        lambda_0 = 1.0;
        q_lambda = 1e-3;
        Q_w = 1*(params.dt)^2*diag([10.0*I0/N_population, 10.0*I0/N_population, q_alpha, q_lambda, q_lambda, q_lambda].^2); % Process noise covariance matrix
        w_bar = zeros(6, 1); % mean value of process noises
        v_bar = 0; % mean value of observation noise
        s_init = [(N_population - I0)/N_population ; I0/N_population ; alpha0 ; lambda_0 ; lambda_0 ; lambda_0]; % initial state vector
        Ps_init = 1.0*(params.dt)^2*diag([10.0*I0/N_population, 10.0*I0/N_population, 10*q_alpha, 10*q_lambda, 10*q_lambda, 10*q_lambda].^2); % Covariance matrix of initial states
        s_final = nan(6, 1); % Set the finite horizon end points (not required during first round)
        Ps_final = nan(6); % Set the finite horizon end points (not required during first round)
        if(isequal(params.obs_type, 'TOTALCASES'))
            R_v = var((NewCasesSmoothed - NewCasesRefined)/N_population); % Observation noise variance
            [~, S_MINUS, S_PLUS, S_SMOOTH, P_MINUS, P_PLUS, P_SMOOTH, ~, ~, rho] = NewCaseEKFEstimatorWithOptimalNPI(control_input, ConfirmedCasesSmoothedNormalized(:)', params, s_init, Ps_init, s_final, Ps_final, w_bar, v_bar, Q_w, R_v, beta, gamma, inv_monitor_len, order);
        elseif(isequal(params.obs_type, 'NEWCASES'))
            R_v = var((NewCasesSmoothed - NewCasesRefined)/N_population); % Observation noise variance
            [~, S_MINUS, S_PLUS, S_SMOOTH, P_MINUS, P_PLUS, P_SMOOTH, ~, ~, rho] = NewCaseEKFEstimatorWithOptimalNPI(control_input, NewCasesSmoothedNormalized(:)', params, s_init, Ps_init, s_final, Ps_final, w_bar, v_bar, Q_w, R_v, beta, gamma, inv_monitor_len, order);
        end
        
        % 2- Apply LASSO over alpha
        x_data_train = NPI_MAXES(:, ones(size(InterventionPlans, 1), 1))' - InterventionPlans;
        %     y_data_train = S_SMOOTH(3, :)'; % alpha
        y_data_train = S_PLUS(3, :)'; % alpha
        
        [B_LASSO, FitInfo] = lasso(x_data_train, y_data_train, 'CV',50);
        idxLambda1SE = FitInfo.Index1SE;
        coef = B_LASSO(:, idxLambda1SE);
        coef0 = FitInfo.Intercept(idxLambda1SE);
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
        s_init = [(N_population - I0)/N_population ; I0/N_population ; S_SMOOTH(3, 1) ; S_SMOOTH(4, 1) ; S_SMOOTH(5, 1) ; S_SMOOTH(6, 1)]; % initial state vector
        [~, S_MINUS2, S_PLUS2, S_SMOOTH2, P_MINUS2, P_PLUS2, P_SMOOTH2, ~, ~, rho2] = NewCaseEKFEstimatorWithOptimalNPI(control_input, NewCasesSmoothedNormalized(:)', params, s_init, Ps_init, s_final, Ps_final, w_bar, v_bar, Q_w, R_v, beta, gamma, inv_monitor_len, 1);
        
        % 4- Apply second LASSO over refined alpha
        x_data_train2 = NPI_MAXES(:, ones(size(InterventionPlans, 1), 1))' - InterventionPlans;
        %         y_data_train2 = S_SMOOTH2(3, :)';
        y_data_train2 = S_PLUS2(3, :)';
        
        [B_LASSO2, FitInfo2] = lasso(x_data_train2, y_data_train2, 'CV',50);
        idxLambda1SE2 = FitInfo2.Index1SE;
        coef_2 = B_LASSO2(:, idxLambda1SE2);
        coef0_2 = FitInfo2.Intercept(idxLambda1SE2);
        y_pred_lasso2 = x_data_train2 * coef_2 + coef0_2;
        %         axTrace2 = lassoPlot(B_LASSO2, FitInfo2);
        %         end
        %
        
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
        MortalityRate = ConfirmedDeathsSmoothed ./ ConfirmedCasesSmoothed;
        MortalityRate(isnan(MortalityRate)) = 0;
        MedMortalityRate = median(MortalityRate);
        MedRecentMortalityRate = median(MortalityRate(round(0.75*end) : end));
        
        CumInfections = cumsum(N_population * S_SMOOTH2(2, :));
        DeathToCumInfectionsRatio = ConfirmedDeaths ./ CumInfections(:);
        
        BetaEstimate = DeathToCumInfectionsRatio/MedRecentMortalityRate;
        MedRecentBetaEstimate = median(BetaEstimate(0.75*end : end));
        
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
        legend('MINUS', 'PLUS', 'SMOOTH', 'MINUS2', 'PLUS2', 'SMOOTH2');
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
        legend('MINUS', 'PLUS', 'SMOOTH', 'MINUS2', 'PLUS2', 'SMOOTH2');
        title(CountryAndRegionList(k), 'interpreter', 'none');
        grid
        
        params.beta
        N_population
        
        close all
    end
end

