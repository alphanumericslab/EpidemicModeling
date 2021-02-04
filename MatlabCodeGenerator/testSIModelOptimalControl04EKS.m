
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

data_file = './../../covid-policy-tracker/data/OxCGRT_latest.csv'; % The data-file cloned from: https://github.com/OxCGRT/covid-policy-tracker/tree/master/data
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

for k = 1:10%3 : NumGeoLocations%index_ : index_%240 %122 : 125%225 %1 : NumGeoLocations
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
end

license('inuse') % list the used licenses (useful for building standalone packages)
