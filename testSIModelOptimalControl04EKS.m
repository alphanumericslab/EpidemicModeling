
close all
clear
clc

% parameters
plot_figures = false; % plot per-region/country plots or not

start_date_criterion = 'DATE_BASED'; %'MIN_CASE_BASED'; % 'MIN_CASE_BASED'/'DATE_BASED'/'DATA_OR_MIN_CASE_BASED'
min_cases = 10; % the minimum cases start date for processing each region/country
start_date = 20200101; % start date
end_date = 20210128; % end date
ar_order = 24; % Autoregressive model
ar_learninghistory = 120; % Autoregressive model estimation history look back
predict_ahead_num_days = 90; % number of days to predict ahead
SmoothingWinLen = 7; % window length used for smoothing the new cases
filter_type = 'MOVINGAVERAGE-CAUSAL'; % 'BYPASS', 'MOVINGAVERAGE-NONCAUSAL' or 'MOVINGAVERAGE-CAUSAL' ' or 'MOVINGMEDIAN' or 'TIKHONOV'; % The last two call functions from the OSET package (oset.ir). Note: 'MOVINGAVERAGE-CAUSAL' is the contest standard and only evaluation algorithm

data_file = './../covid-policy-tracker/data/OxCGRT_latest.csv'; % The data-file cloned from: https://github.com/OxCGRT/covid-policy-tracker/tree/master/data
% data_file = './../covid-xprize/AlphanumericsTeam/data/files/OxCGRT_latest_aug.csv'; % The data-file cloned from: https://github.com/OxCGRT/covid-policy-tracker/tree/master/data
included_IP = {'C1_School closing',... % see full descriptions at: https://github.com/OxCGRT/covid-policy-tracker/blob/master/documentation/codebook.md%'C1_Flag',...
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

NPI_MAXES = [3, 3, 2, 4, 2, 3, 2, 4, 2, 3, 2, 4]'; % Maximum values for the 12 phase 2 NPIs

% Read Oxford time-series data
import_opt = detectImportOptions(data_file);
import_opt.VariableNamingRule = 'preserve';
% import_opt.VariableTypes(16) = {'double'}; % C5_Flag on column 16 is of type double (the format is incorrectly read as 'char')
import_opt.VariableTypes(6:end) = {'double'}; % All except the first five columns are 'double' (correct if not correctly identified)
all_data = readtable(data_file, import_opt);%'ReadVariableNames', true, 'VariableNamingRule', 'preserve', 'TextType', 'string', 'detectImportOptions', 'true');
ColumnHeaders = all_data.Properties.VariableNames; % Column titles
AllCountryNames = all_data.CountryName;%(:, ismember(ColumnHeaders,'CountryName')); % All country names
AllCountryCodes = all_data.CountryCode;%(:, ismember(ColumnHeaders,'CountryCode')); % All country codes
AllRegionNames = all_data.RegionName;%(:, ismember(ColumnHeaders,'RegionName')); % All region names
AllRegionCodes = all_data.RegionCode;%(:, ismember(ColumnHeaders,'RegionCode')); % All region codes

% Make GeoID code (A combination of region and country codes)
GeoID = strcat(string(AllCountryCodes), string(AllRegionCodes));

% Make a dictionary of country-regions
[CountryAndRegionList, IA, IC] = unique(GeoID, 'stable');

NumGeoLocations = length(CountryAndRegionList); % Number of country-region pairs

% Select regions/countries
index_ = find(CountryAndRegionList == 'USA') % 241 Maryland
N_population = 331e6; % US population
%     N_population = 6.046e6; % Maryland population
% index_ = find(CountryAndRegionList == 'FRA')
% index_ = find(CountryAndRegionList == 'SPA')
% N_population = 67.06e6; % France population
% N_population = 83.02e6; % Germany population
% N_population = 46.94e6; % Spain population
NumNPI = length(included_IP);
% NPI weights
npi_weights = ones(NumNPI, 1);
npi_weights = NumNPI*npi_weights/sum(npi_weights);

human_npi_cost_factor = 0.001; % [0, 1]: 0 neglects NPI cost and 1 neglects human factor!

for k = 3 : NumGeoLocations%index_ : index_%240 %122 : 125%225 %1 : NumGeoLocations
    k
    switch start_date_criterion
        case 'MIN_CASE_BASED'
            geoid_all_row_indexes = GeoID == CountryAndRegionList(k) & all_data.ConfirmedCases > min_cases & all_data.Date <= end_date;
        case 'DATE_BASED'
            geoid_all_row_indexes = GeoID == CountryAndRegionList(k) & all_data.Date >= start_date & all_data.Date <= end_date;
        case 'DATA_OR_MIN_CASE_BASED'
            geoid_all_row_indexes = GeoID == CountryAndRegionList(k) & all_data.ConfirmedCases > min_cases & all_data.Date >= start_date & all_data.Date <= end_date;
    end
    ConfirmedCases = all_data.ConfirmedCases(geoid_all_row_indexes); % Confirmed cases
    ConfirmedDeaths = all_data.ConfirmedDeaths(geoid_all_row_indexes); % Death cases
    
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
    InterventionPlansAvg= nanmean(InterventionPlans, 2);
    InterventionPlansAvgInt = cumsum(InterventionPlansAvg);
    InterventionPlansAvgIntNorm = (InterventionPlansAvgInt - nanmean(InterventionPlansAvgInt))/nanstd(InterventionPlansAvgInt);
    
    % Calculate the number of new cases
    NewCases = [nan; diff(ConfirmedCases)]; % calculate the new daily cases
    NewCases(NewCases < 0) = 0;
    
    if(length(NewCases) < 2)
        warning(strcat("Insufficient data for ", CountryAndRegionList(k), ". Skipping the country/region."));
        continue;
    end
    
    % Replace nans with 0 and the last one with the previous number
    NewCasesRefined = NewCases;
    NewCasesRefined(isnan(NewCasesRefined)) = 0;
    if(isnan(NewCases(end))) % Just fill the last missing date (if any) with the previous day
        NewCasesRefined(end) = NewCasesRefined(end-1);
    end
    
    % Smooth the newcases time series
    switch filter_type
        case 'BYPASS' % Bypass (no smoothing)
            NewCasesSmoothed = NewCasesRefined;
        case 'TIKHONOV' % Tikhonov regularization
            % DiffOrderOrFilterCoefs = [1 -2 1]; % Smoothness filter coefs
            DiffOrderOrFilterCoefs = 2; % Smoothness filter order
            TikhonovGamma = 25.0; % Tikhonov roughness penalty
            NewCasesSmoothed = TikhonovRegularization(NewCasesRefined', DiffOrderOrFilterCoefs, TikhonovGamma)';
        case 'MOVINGAVERAGE-CAUSAL'
            NewCasesSmoothed = filter(ones(1, SmoothingWinLen), SmoothingWinLen, NewCasesRefined); % causal
        case 'MOVINGAVERAGE-NONCAUSAL'
            NewCasesSmoothed = BaseLine1(NewCasesRefined', SmoothingWinLen, 'mn')'; % non-causal (zero-phase)
        case 'MOVINGMEDIAN'
            NewCasesSmoothed = BaseLine1(NewCasesRefined', floor(SmoothingWinLen/2), 'md')';
            NewCasesSmoothed = BaseLine1(NewCasesSmoothed', SmoothingWinLen, 'mn')';
        otherwise
            error('Unknown filter type');
    end
    % % % % %     NewCasesSmoothed(end - predict_ahead_num_days + 1 : end) = nan; % number of days to forecast
    
    % 1- Apply an EKF to find first estimates of alpha (with zero control assumption)
    %     I0 = min_cases; % initial number of cases
    I0 = max(1, NewCasesSmoothed(1)); % initial number of cases
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
    params.gamma = 1/(params.dt*100); % input to contact influence rate (inverse time)
    params.beta = 1/(params.dt*75); % recovery rate (inverse time)
    params.sigma = 100000; % sigmoid function slope
    beta = 0.9; % Observation noise update factor (set to 1 for no update)
    gamma = 0.995; % Kalman gain stability factor (set very close to 1, or equal to 1 to disable the feature)
    inv_monitor_len = 21; % Window length for innovations process whiteness monitoring
    order = 1; % 1 for standard EKF; 2 for second-order EKF
    Q_w = (params.dt)^2*diag([0.01, 0.01, 0.1, 10, 10, 10].^2); % Process noise covariance matrix
    R_v = 1e-6; % Observation noise variance
    w_bar = zeros(6, 1); % mean value of process noises
    v_bar = 0; % mean value of observation noise
    NewCasesSmoothedNormalized = NewCasesSmoothed / N_population;
    alpha0 = 0.01;
    s_init = [(N_population - I0)/N_population ; I0/N_population ; alpha0 ; 1.0 ; 1.0 ; 1.0]; % initial state vector
    Ps_init = 1000 * Q_w; % Covariance matrix of initial states
    s_final = nan(6, 1); % Set the finite horizon end points (not required during first round)
    Ps_final = nan(6); % Set the finite horizon end points (not required during first round)
    [~, S_MINUS, S_PLUS, P_MINUS, P_PLUS, K_GAIN, S_SMOOTH, P_SMOOTH, innovations, rho] = NewCaseEKFEstimatorWithOptimalNPI(control_input, NewCasesSmoothedNormalized(:)', params, s_init, Ps_init, s_final, Ps_final, w_bar, v_bar, Q_w, R_v, beta, gamma, inv_monitor_len, 1);
    
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
    
    TWO_ROUND_TRAINING = false;
    if(TWO_ROUND_TRAINING)
        % 3- Apply the SECOND round EKF to refine alpha based on real inputs
        I0 = max(1, round(N_population * S_SMOOTH(2, 1)));
        params.a = coef; % input influence weight vector
        params.b = coef0; % input influence bias constant
        control_input = InterventionPlans'; % The second round works with real inputs
        s_init = [(N_population - I0)/N_population ; I0/N_population ; S_SMOOTH(3, 1) ; S_SMOOTH(4, 1) ; S_SMOOTH(5, 1) ; S_SMOOTH(6, 1)]; % initial state vector
        [~, S_MINUS2, S_PLUS2, P_MINUS2, P_PLUS2, K_GAIN2, S_SMOOTH2, P_SMOOTH2, innovations2, rho2] = NewCaseEKFEstimatorWithOptimalNPI(control_input, NewCasesSmoothedNormalized(:)', params, s_init, Ps_init, s_final, Ps_final, w_bar, v_bar, Q_w, R_v, beta, gamma, inv_monitor_len, 1);
        
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
    end
    
    if(TWO_ROUND_TRAINING)
        %         I0 = N_population * S_SMOOTH2(2, 1);
        params.a = coef_2;% + eps; % input influence weight vector
        params.b = coef0_2; % input influence bias constant
        %         s_init = [(N_population - I0)/N_population ; I0/N_population ; S_SMOOTH2(3, 1) ; S_SMOOTH2(4, 1) ; S_SMOOTH2(5, 1) ; S_SMOOTH2(6, 1)]; % initial state vector
    else
        %         I0 = N_population * S_SMOOTH(2, 1);
        params.a = coef;% + eps; % input influence weight vector
        params.b = coef0; % input influence bias constant
        %         s_init = [(N_population - I0)/N_population ; I0/N_population ; S_SMOOTH(3, 1) ; S_SMOOTH(4, 1) ; S_SMOOTH(5, 1) ; S_SMOOTH(6, 1)]; % initial state vector
    end
    
    % Forecast future alpha values
    if(TWO_ROUND_TRAINING)
        ar_train_segment = y_data_train2(end - ar_learninghistory + 1 : end);
    else
        ar_train_segment = y_data_train(end - ar_learninghistory + 1 : end);
    end
    ar_sys = ar(ar_train_segment, ar_order);
    A_ar = get(ar_sys, 'A');
    noisevar_ar = get(ar_sys, 'NoiseVariance');
    zi = filtic(sqrt(noisevar_ar), A_ar, ar_train_segment(end:-1:1));
    y_pred_ar = filter(sqrt(noisevar_ar), A_ar, randn(1, predict_ahead_num_days), zi)';
    AlphaHatARX = [ar_train_segment ; y_pred_ar]';
    AlphaHatARX(AlphaHatARX < 0 ) = 0;
    
    figure
    plot(AlphaHatARX)
    hold on
    plot(ar_train_segment, 'r')
    grid
    
    % Generate a test scenario with the same set of weights
    num_random_input_monte_carlo_runs = 10;
    J0 = zeros(1, num_random_input_monte_carlo_runs);
    J1 = zeros(1, num_random_input_monte_carlo_runs);
    J = zeros(1, num_random_input_monte_carlo_runs);
    
    for kk = 1 : num_random_input_monte_carlo_runs
%         u = zeros(NumNPI, length(NewCasesSmoothedNormalized));
        u = zeros(NumNPI, predict_ahead_num_days);
        for t = 1 : predict_ahead_num_days
            for jj = 1 : NumNPI
                u(jj, t) = randi([params.u_min(jj), params.u_max(jj)]);
            end
        end
        
        AlphaHatARX = [ar_train_segment ; y_pred_ar + params.gamma*(u' * params.a + params.b)]';
        AlphaHatARX(AlphaHatARX < 0 ) = 0;
        
%         [s, i] = SI_Controlled(AlphaHatARX, params.beta, (N_population - I0)/N_population, I0/N_population, length(AlphaHatARX), params.dt);
        [s, i] = SI_Controlled(AlphaHatARX, params.beta, (N_population - I0)/N_population, I0/N_population, length(AlphaHatARX), params.dt);
        
    hold on
    plot(AlphaHatARX, 'g')
    grid

    figure
    hold on
    plot(s);
    plot(i);
    grid
        %alpha = S_SMOOTH_opt(3, :);
        %         [s, i, alpha] = SI_Controlled(u, params.u_max, params.alpha_min, params.alpha_max, params.gamma, S_SMOOTH_opt(3, 1), (N_population - I0)/N_population, I0/N_population, length(NewCasesSmoothedNormalized), params.dt);
        %[s, i] = SI_Controlled(alpha, params.beta, (N_population - I0)/N_population, I0/N_population, length(NewCasesSmoothedNormalized), params.dt);
        %[J0(kk), J1(kk), J(kk)] = NPICost(s.*i.*alpha, control_input_opt, diag(npi_weights) * ones(NumNPI, size(control_input_opt, 2)), human_npi_cost_factor);
    end

    figure
    hold on
    plot(J0, J1, 'bo');
    plot(J0_opt, J1_opt, 'ro');
    grid
    
    
    % 5- Find the optimal controls for the trained model
    % No inputs at this stage. The algorithm predicts the optimal inputs
    control_input = InterventionPlans';
    control_input(:, end - predict_ahead_num_days + 1 : end) = nan;
    %     params.w = diag(npi_weights) * ones(NumNPI, length(NewCasesSmoothed));
    params.w = npi_weights; % Get the NPI weights
    
    % % %     % Set the finite horizon end points
    % % %     s_final = nan(6, 1);
    % % %     % % %     s_final(2) = 0; % zero infection rates
    % % %     s_final(4) = 0; % zero costates
    % % %     s_final(5) = 0; % zero costates
    % % %     s_final(6) = 0; % zero costates
    % % %     Ps_final = zeros(6);
    % % %     % %     Ps_final(2, 2) = 1e-8; % zero infection rates
    % % %     %     Ps_final(4 : 6, 4 : 6) = 1e-6;
    % % %     Ps_final(1, 1) = nan; % don't touch these end-point states during smoothing
    % % %     Ps_final(2, 2) = nan; % don't touch these end-point states during smoothing
    % % %     Ps_final(3, 3) = nan; % don't touch these end-point states during smoothing
    % % %     Ps_final(4, 4) = 1e-2; % zero costates
    % % %     Ps_final(5, 5) = 1e-2; % zero costates
    % % %     Ps_final(6, 6) = 1e-2; % zero costates
    [control_input_opt, S_MINUS_opt, S_PLUS_opt, P_MINUS_opt, P_PLUS_opt, K_GAIN_opt, S_SMOOTH_opt, P_SMOOTH_opt, innovations_opt, rho_opt] = NewCaseEKFEstimatorWithOptimalNPI(control_input, NewCasesSmoothedNormalized(:)', params, s_init, Ps_init, s_final, Ps_final, w_bar, v_bar, Q_w, R_v, beta, gamma, inv_monitor_len, 1);
    
    
    
    
    
    
    
    
    % % %     figure
    % % %     subplot(211)
    % % %     plot(InterventionPlans)
    % % %     grid
    % % %     subplot(212)
    % % %     plot(control_input_opt')
    % % %     grid
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    figure
    subplot(211)
    plot(N_population * NewCasesSmoothedNormalized, 'linewidth', 2);
    hold on
    plot(N_population * S_PLUS(1, :).*S_PLUS(2, :).*S_PLUS(3, :));
    plot(N_population * S_SMOOTH(1, :).*S_SMOOTH(2, :).*S_SMOOTH(3, :));
    plot(N_population * S_PLUS2(1, :).*S_PLUS2(2, :).*S_PLUS2(3, :));
    plot(N_population * S_SMOOTH2(1, :).*S_SMOOTH2(2, :).*S_SMOOTH2(3, :));
    legend('NewCases', 'PLUS', 'SMOOTH', 'PLUS2', 'SMOOTH2');
    title(CountryAndRegionList(k), 'interpreter', 'none');
    grid;
    
    subplot(212);
    plot(y_data_train);
    hold on;
    plot(y_data_train2);
    plot(y_pred_lasso);
    plot(y_pred_lasso2);
    grid
    legend('alpha', 'alpha 2', 'alpha lasso', 'alpha lasso 2');
    
    % % %     reply = input('Do you want more regiond? Y/N [Y]:','s');
    % % %     if isempty(reply)
    % % %         reply = 'Y';
    % % %     end
    % % %     if(~isequal(upper(reply), 'Y'))
    % % %         break;
    % % %     end
    close all
    
    if(0)
        
        new_case_opt = S_SMOOTH_opt(1, :) .* S_SMOOTH_opt(2, :) .* S_SMOOTH_opt(3, :);
        [J0_opt, J1_opt, J_opt] = NPICost(new_case_opt, control_input_opt, diag(npi_weights) * ones(NumNPI, size(control_input_opt, 2)), human_npi_cost_factor);
        
        figure
        subplot(411);
        plot(N_population * NewCasesSmoothedNormalized);
        hold on
        plot(N_population * S_PLUS(1, :).*S_PLUS(2, :).*S_PLUS(3, :));
        plot(N_population * S_SMOOTH(1, :).*S_SMOOTH(2, :).*S_SMOOTH(3, :));
        plot(N_population * S_PLUS_opt(1, :).*S_PLUS_opt(2, :).*S_PLUS_opt(3, :));
        plot(N_population * S_SMOOTH_opt(1, :).*S_SMOOTH_opt(2, :).*S_SMOOTH_opt(3, :));
        grid;
        legend('NewCasesSmoothed', 'PLUS', 'SMOOTH', 'PLUS-opt', 'SMOOTH-opt');
        %     legend('NewCasesSmoothed', 'PLUS', 'SMOOTH');
        title(CountryAndRegionList(k), 'interpreter', 'none');
        
        subplot(412);
        plot(N_population * S_PLUS_opt(1, :));
        hold on
        plot(N_population * S_SMOOTH_opt(1, :));
        plot(N_population * S_PLUS_opt(1, :));
        plot(N_population * S_SMOOTH_opt(1, :));
        grid;
        legend('S-PLUS', 'S-SMOOTH', 'S-PLUS-opt', 'S-SMOOTH-opt');
        %     legend('S-PLUS', 'S-SMOOTH');
        %     title('S');
        
        subplot(413);
        plot(N_population * S_PLUS_opt(2, :));
        hold on
        plot(N_population * S_SMOOTH_opt(2, :));
        plot(N_population * S_PLUS_opt(2, :));
        plot(N_population * S_SMOOTH_opt(2, :));
        grid;
        legend('I-PLUS', 'I-SMOOTH', 'I-PLUS-opt', 'I-SMOOTH-opt');
        %     legend('I-PLUS', 'I-SMOOTH');
        %     title('I');
        
        subplot(414);
        plot(S_PLUS_opt(3, :));
        hold on
        plot(S_SMOOTH_opt(3, :));
        plot(S_PLUS_opt(3, :));
        plot(S_SMOOTH_opt(3, :));
        grid;
        legend('Alpha-PLUS', 'Alpha-SMOOTH', 'Alpha-PLUS-opt', 'Alpha-SMOOTH-opt');
        %     legend('Alpha-PLUS', 'Alpha-SMOOTH');
        %     title('Alpha');
        
        figure
        hold on
        plot(rho);
        plot(rho_opt);
        grid;
        legend('Rho', 'Rho-opt');
        %     legend('Rho', 'Rho2');
        %     legend('Rho');
        
        
        
        
        FEATURE_XCORR_ANALYSIS = false;
        if(FEATURE_XCORR_ANALYSIS)
            xycorrs = zeros(2*length(y_data_train) - 1, size(x_data_train, 2));
            for ii = 1 : size(x_data_train, 2)
                xycorrs(:, ii) = xcorr(y_data_train - mean(y_data_train), x_data_train(:, ii) - mean(x_data_train(:, ii)));
                xycorrs(:, ii) = xycorrs(:, ii) / max(xycorrs(:, ii));
            end
            yycorr = xcorr(y_data_train, y_data_train);
            yycorr = yycorr / max(yycorr);
            
            figure
            plot(xycorrs);
            hold on
            plot(yycorr, 'k', 'linewidth', 2);
            grid
        end
        
        % % % %     % 5- Apply the THIRD round EKF to find optimal control
        
        % input to alpha filter model:
        % % %     numerator = [0, params.dt * params.gamma];
        % % %     denominator = [1, -(1 - params.dt * params.gamma)];
        % % %     x_data_train_filtered = filter(numerator, denominator, x_data_train, [], 1);
        
        figure
        plot(x_data_train);
        %     plot(x_data_train_filtered);
        hold on
        plot(max(x_data_train(:))*y_data_train/max(y_data_train), 'color', 'k', 'linewidth', 3);
        %     plot(y_data_train, 'color', 'k', 'linewidth', 3);
        grid
        
        
        % % %     % SECOND ROUND
        % % %     control_input = InterventionPlans';
        % % %     control_input(:, end - 100 : end) = nan;
        % % %     s2_init = [NewCasesSmoothedNormalized(1) ; I0/N_population ; S_SMOOTH_opt(3, 1) ; S_SMOOTH_opt(4, 1) ; S_SMOOTH_opt(5, 1) ; S_SMOOTH_opt(6, 1)]; % initial state vector
        % % %     Ps2_init = 0.1 * Ps2_init; % Process noise covariance matrix
        % % %     [control_input_opt, S_MINUS_opt, S_PLUS_opt, P_MINUS_opt, P_PLUS_opt, K_GAIN_opt, S_SMOOTH_opt, P_SMOOTH_opt, innovations_opt, rho_opt] = NewCaseEKFEstimatorWithOptimalNPI(control_input, NewCasesSmoothedNormalized(:)', params, s2_init, Ps2_init, s_final, Ps_final, w2_bar, v2_bar, Q2_w, R2_v, beta, gamma, inv_monitor_len, 1);
        % % %     params.a = 0.01 * rand(NumNPI, 1); % input influence weight vector
        
        % % %     % THIRD ROUND
        % % %     s2_init = [NewCasesSmoothedNormalized(1) ; I0/N_population ; S_SMOOTH_opt(3, 1) ; S_SMOOTH_opt(4, 1) ; S_SMOOTH_opt(5, 1) ; S_SMOOTH_opt(6, 1)]; % initial state vector
        % % %     Ps2_init = 0.1 * Ps2_init; % Process noise covariance matrix
        % % %     [control_input_opt, S_MINUS_opt, S_PLUS_opt, P_MINUS_opt, P_PLUS_opt, K_GAIN_opt, S_SMOOTH_opt, P_SMOOTH_opt, innovations_opt, rho_opt] = NewCaseEKFEstimatorWithOptimalNPI(control_input, NewCasesSmoothedNormalized(:)', params, s2_init, Ps2_init, s_final, Ps_final, w2_bar, v2_bar, Q2_w, R2_v, beta, gamma, inv_monitor_len, 1);
        
        
        figure
        hold on
        %     plot(S_PLUS_opt2(4, :));
        %     plot(S_PLUS_opt2(5, :));
        %     plot(S_PLUS_opt2(6, :));
        plot(S_SMOOTH_opt(4, :));
        plot(S_SMOOTH_opt(5, :));
        plot(S_SMOOTH_opt(6, :));
        plot(S_SMOOTH_opt2(4, :));
        plot(S_SMOOTH_opt2(5, :));
        plot(S_SMOOTH_opt2(6, :));
        legend('\lambda_1 S', '\lambda_2 S', '\lambda_3 S', '\lambda_1 S2', '\lambda_2 S2', '\lambda_3 S2');
        grid;
        title('Co-states');
        
        
        figure
        hold on
        plot(control_input_opt2' - InterventionPlans);
        %     plot(control_input_opt', 'b');
        %     plot(control_input_opt2', 'r');
        title('difference in control inputs');
        grid
        
        % Find inclining/declining time instants
        %     rate_incline_decline_index = Lambda_GeoGenRatios;
        %     rate_incline_decline_index = Lambda_NonLinLS;
        rate_incline_decline_index = S_SMOOTH(2, :);
        
        rates_inclining = find(rate_incline_decline_index >= 0);
        rates_declining = find(rate_incline_decline_index < 0);
        
        if(plot_figures)
            figure
            plot(rho);
            grid
            title(CountryAndRegionList(k), 'interpreter', 'none');
            set(gca, 'fontsize', 18);
            set(gca, 'box', 'on');
            ylabel('Rho');
            xlabel('Days since 100th case');
            
            lgn = {};
            ksigma = 3.0; % k-sigma envelope plots
            %     options.handle = figure;
            figure
            hold on
            %     errorbar(S_MINUS(1, :), ksigma*sqrt(squeeze(P_MINUS(1, 1, :)))); lgn = cat(2, lgn, {'S_MINUS'});
            errorbar(S_PLUS(1, :), ksigma*sqrt(squeeze(P_PLUS(1, 1, :)))); lgn = cat(2, lgn, {'\pm 3\sigma EKF Envelopes'});
            %     errorbar(S_PLUS2(1, :), ksigma*sqrt(squeeze(P_PLUS2(1, 1, :)))); lgn = cat(2, lgn, {'\pm 3\sigma EKF2 Envelopes'});
            %     plot_areaerrorbar(S_PLUS(1, :), ksigma*sqrt(squeeze(P_PLUS(1, 1, :))), options); lgn = cat(2, lgn, {'EKF'});
            errorbar(S_SMOOTH(1, :), ksigma*sqrt(squeeze(P_SMOOTH(1, 1, :)))); lgn = cat(2, lgn, {'\pm 3\sigma EKS Envelopes'});
            %     errorbar(S_SMOOTH2(1, :), ksigma*sqrt(squeeze(P_SMOOTH2(1, 1, :)))); lgn = cat(2, lgn, {'\pm 3\sigma EKS2 Envelopes'});
            %     plot_areaerrorbar(S_SMOOTH(1, :), ksigma*sqrt(squeeze(P_SMOOTH(1, 1, :))), options); lgn = cat(2, lgn, {'EKS'});
            plot(NewCases, 'linewidth', 1); lgn = cat(2, lgn, {'New Cases'});
            plot(NewCasesSmoothed, 'linewidth', 3); lgn = cat(2, lgn, {'7-Day MA'});
            plot(S_PLUS(1, :), 'linewidth', 3); lgn = cat(2, lgn, {'EKF'});
            %     plot(S_PLUS2(1, :), 'linewidth', 3); lgn = cat(2, lgn, {'EKF2'});
            plot(S_SMOOTH(1, :), 'linewidth', 3); lgn = cat(2, lgn, {'EKS'});
            %     plot(S_SMOOTH2(1, :), 'linewidth', 3); lgn = cat(2, lgn, {'EKS2'});
            legend(lgn);%, 'interpreter', 'none');
            title(CountryAndRegionList(k), 'interpreter', 'none');
            xlabel('Days since 100th case');
            ylabel('Number of new cases');
            set(gca, 'fontsize', 18);
            set(gca, 'box', 'on');
            axis tight
            grid
            
            lgn = {};
            ksigma = 1.0;
            figure
            hold on
            %     errorbar(S_MINUS2(2, :), ksigma*sqrt(squeeze(P_MINUS2(2, 2, :)))); lgn = cat(2, lgn, {'S_MINUS2'});
            errorbar(S_PLUS(2, :), ksigma*sqrt(squeeze(P_PLUS(2, 2, :)))); lgn = cat(2, lgn, {'\pm 3\sigma EKF Envelopes'});
            %     errorbar(S_PLUS2(2, :), ksigma*sqrt(squeeze(P_PLUS2(2, 2, :)))); lgn = cat(2, lgn, {'\pm 3\sigma EKF2 Envelopes'});
            errorbar(S_SMOOTH(2, :), ksigma*sqrt(squeeze(P_SMOOTH(2, 2, :)))); lgn = cat(2, lgn, {'\pm 3\sigma EKS Envelopes'});
            %     errorbar(S_SMOOTH2(2, :), ksigma*sqrt(squeeze(P_SMOOTH2(2, 2, :)))); lgn = cat(2, lgn, {'\pm 3\sigma EKS2 Envelopes'});
            plot(S_PLUS(2, :), 'linewidth', 3); lgn = cat(2, lgn, {'EKF'});
            %     plot(S_PLUS2(2, :), 'linewidth', 3); lgn = cat(2, lgn, {'EKF2'});
            plot(S_SMOOTH(2, :), 'linewidth', 3); lgn = cat(2, lgn, {'EKS'});
            %     plot(S_SMOOTH2(2, :), 'linewidth', 3); lgn = cat(2, lgn, {'EKS2'});
            %     plot(Lambda_GeoGenRatiosSmoothed, 'linewidth', 3); lgn = cat(2, lgn, {'Lambda_GeoGenRatiosSmoothed'});
            legend(lgn);%, 'interpreter', 'none');
            title(CountryAndRegionList(k), 'interpreter', 'none');
            xlabel('Days since 100th case');
            ylabel('\lambda_k');
            set(gca, 'fontsize', 18);
            set(gca, 'box', 'on');
            axis tight
            grid
        end
    end
end

license('inuse') % list the used licenses (useful for building standalone packages)
