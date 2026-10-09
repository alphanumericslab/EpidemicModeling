function tests = test_growth_estimators()
    % TEST_GROWTH_ESTIMATORS Recover known growth and inspect forecast conventions.
    setup_once([]);
    tests = functiontests(localfunctions);
end

function setup_once(test_case)
    % SETUP_ONCE Add function library.
    addpath(fullfile(fileparts(fileparts(mfilename('fullpath'))), 'functions'));
end

function test_log_exponential(test_case)
    % TEST_LOG_EXPONENTIAL Recover exact amplitude and per-day growth.
    x = 25 * exp(.035 * (0:39));
    [~, a, g, fit] = rt_exp_fit_log_lin_reg(x, 7, 1);
    verifyEqual(test_case, g(7:end), .035 * ones(1, 34), 'AbsTol', 1e-12);
    verifyEqual(test_case, a(7:end), x(7:end), 'RelTol', 1e-12);
    verifyEqual(test_case, fit(7:end), x(7:end) * exp(.035), 'RelTol', 1e-12);
end

function test_generation_ratio(test_case)
    % TEST_GENERATION_RATIO Lagged exponential ratios recover constant growth.
    x = 25 * exp(.035 * (0:39));
    [~, g, ~, gs] = rt_exp_fit_gen_ratios(x, 7, 5, 1);
    verifyEqual(test_case, g(1:5), zeros(1, 5));
    verifyEqual(test_case, g(6:end), .035 * ones(1, 35), 'AbsTol', 1e-12);
    verifyEqual(test_case, gs(12:end), .035 * ones(1, 29), 'AbsTol', 1e-12);
end

function test_nonlinear_fit(test_case)
    % TEST_NONLINEAR_FIT Optional nlinfit recovers the known exponential.
    assumeTrue(test_case, exist('nlinfit', 'file') == 2, ...
        'Statistics and Machine Learning Toolbox required.');
    x = 25 * exp(.035 * (0:39));
    [~, a, g] = rt_exp_fit_nonlin_ls(x, 7, 1);
    verifyEqual(test_case, g(7:end), .035 * ones(1, 34), 'AbsTol', 1e-5);
    verifyEqual(test_case, a(7:end), x(7:end), 'RelTol', 1e-5);

    x(11) = 0;
    [~, a, g] = rt_exp_fit_nonlin_ls(x, 8, 1, false);
    verifyEqual(test_case, g(10), 0);
    verifyEqual(test_case, a(10), x(10));
end
