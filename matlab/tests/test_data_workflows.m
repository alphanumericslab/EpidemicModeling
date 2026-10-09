function tests = test_data_workflows()
    % TEST_DATA_WORKFLOWS Native CSV/JSON integration tests for the upgraded pipeline.
    % Author: Reza Sameni | Emory University
    addpath(fullfile(fileparts(fileparts(mfilename('fullpath'))), 'functions'));
    tests = functiontests(localfunctions);
end

function test_csv_training_prescription_and_scoring(test_case)
    % TEST_CSV_TRAINING_PRESCRIPTION_AND_SCORING Train without future cases and validate
    % daily output.
    folder = tempname;
    mkdir(folder);
    cleanup = onCleanup(@() rmdir(folder, 's'));
    day = 0:49;
    u = [2 * (day >= 25); 2 * (day >= 35)];

    truth.population = 1e6;
    truth.params = default_si_params([3; 3]);
    truth.params.a = [.06; .04];
    truth.params.b = .02;
    truth.state = [.9999; .0001; .3];

    truth.covariance = diag([1e-8 1e-8 .01]);
    [cases, ~] = forecast_npi(truth, u);
    data = table(repmat("Example", 50, 1), repmat("", 50, 1), string(datetime(2020, 1, 1) + ...
        days((0:49)'), 'yyyyMMdd'), cumsum(cases)', u(1, :)', u(2, :)', 'VariableNames', ...
        {'CountryName', 'RegionName', 'Date', 'ConfirmedCases', 'NPI A', 'NPI B'});
    geo = table("Example", "", 'VariableNames', {'CountryName', 'RegionName'});
    population = [geo table(1e6, 'VariableNames', {'Population2020'})];

    costs = [geo table(1, 1.5, 'VariableNames', {'NPI A', 'NPI B'})];
    data_file = fullfile(folder, 'data.csv');
    geo_file = fullfile(folder, 'geo.csv');
    pop_file = fullfile(folder, 'pop.csv');
    cost_file = fullfile(folder, 'cost.csv');

    model_file = fullfile(folder, 'model.json');
    writetable(data, data_file);
    writetable(geo, geo_file);
    writetable(population, pop_file);
    writetable(costs, cost_file);

    bundle = train_npi_prescriptor('2020-01-01', '2020-02-09', data_file, geo_file, ...
        pop_file, {'NPI A', 'NPI B'}, [3; 3], model_file);
    verifyEqual(test_case, numel(bundle.models), 1);
    verifyEqual(test_case, bundle.models{1}.training_days, 40);
    output = prescribe_npi('2020-02-10', '2020-02-14', geo_file, cost_file, fullfile(folder, ...
        'output.csv'), model_file);
    verifyEqual(test_case, height(output), 5);

    verifyEqual(test_case, output.Date(1), "2020-02-10");
    scores = forecast_quality_assessment([1; 1.5], .3, '2020-01-01', '2020-02-09', ...
        '2020-01-01', '2020-02-19', 10, data_file, geo_file, pop_file, {'NPI A', 'NPI B'}, ...
        [0; 0], [3; 3], model_file);
    verifyEqual(test_case, scores.days, 10);
    verifyTrue(test_case, isfinite(scores.rmse));
    results = train_predict_prescribe_npi([1; 1.5], .3, '2020-01-01', '2020-02-09', ...
        '2020-01-01', '2020-02-19', data_file, geo_file, pop_file, {'NPI A', 'NPI B'}, [0; 0], ...
        [3; 3], model_file);

    verifyEqual(test_case, numel(results), 1);
    verifySize(test_case, results{1}.fixed_cases, [1 10]);
end

function test_wide_country_aggregation(test_case)
    % TEST_WIDE_COUNTRY_AGGREGATION JHU metadata, absent countries and one-based thresholds.
    folder = tempname;
    mkdir(folder);
    cleanup = onCleanup(@() rmdir(folder, 's'));
    files = cell(1, 3);
    arrays = {[0 1 3; 0 2 4], [0 0 1; 0 0 0], [0 0 0; 0 1 1]};

    for k = 1:3
        values = arrays{k};
        table_data = table(["A"; "B"], ["Example"; "Example"], [0; 0], [0; 0], values(:, ...
            1), values(:, 2), values(:, 3), 'VariableNames', {'Province/State', ...
            'Country/Region', 'Lat', 'Long', '1/1/20', '1/2/20', '1/3/20'});
        files{k} = fullfile(folder, sprintf('wide_%d.csv', k));
        writetable(table_data, files{k});
    end

    [total, infected, recovered, deceased, first, threshold, ...
        count] = read_covid19_data(files{:}, ["Example", "Absent"], 5);
    verifyEqual(test_case, total, [0 3 7; 0 0 0]);
    verifyEqual(test_case, first, [2 0]);
    verifyEqual(test_case, threshold, [3 0]);
    verifyEqual(test_case, count, 3);

    verifyEqual(test_case, infected, total - recovered - deceased);
end
