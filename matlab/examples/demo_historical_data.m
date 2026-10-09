function demo_historical_data()
    % DEMO_HISTORICAL_DATA Example 05: historical Oxford data and held-out forecast.
    % Uses only the bundled CSVs; no external repositories or downloads are needed.
    % Author: Reza Sameni | Emory University
    root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
    data = read_oxford_data(fullfile(root, 'data', 'OxCGRT_latest.csv'));
    series = data(data.CountryName == "France" & data.RegionName == "", :);
    [daily, smoothed] = prepare_cases(series.ConfirmedCases);
    figure('Color', 'w');

    plot(series.Date, [daily; smoothed]', 'LineWidth', 1.5);
    grid on;
    xlabel('Date');
    ylabel('Cases/day');
    title('France: historical reports');

    legend('Daily reports', 'Causal 7-day mean');
    columns = {'C1_School closing', 'C2_Workplace closing', 'C3_Cancel public events', ...
        'C4_Restrictions on gatherings', 'C5_Close public transport', ...
        'C6_Stay at home requirements', 'C7_Restrictions on internal movement', ...
        'C8_International travel controls', 'H1_Public information campaigns', ...
        'H2_Testing policy', 'H3_Contact tracing', 'H6_Facial Coverings'};
    maximum = [3; 3; 2; 4; 2; 3; 2; 4; 2; 3; 2; 4];
    sample = series(end - 113:end, :);
    train = sample(1:end - 14, :);

    test = sample(end - 13:end, :);
    u = fillmissing(sample{:, columns}, 'previous');
    u(isnan(u)) = 0;
    populations = read_geo_table(fullfile(root, 'data', 'populations.csv'));
    selected = populations.CountryName == "France" & populations.RegionName == "";

    model = fit_npi_model(train.ConfirmedCases, u(1:end - 14, :)', ...
        populations.Population2020(selected), maximum);
    [predicted, ~] = forecast_npi(model, u(end - 13:end, :)');
    actual = max(0, diff([train.ConfirmedCases(end); test.ConfirmedCases]));
    error = predicted(:) - actual;
    fprintf('Held-out MAE: %.2f; RMSE: %.2f\n', mean(abs(error)), sqrt(mean(error.^2)));

    figure('Color', 'w');
    plot(test.Date, [actual predicted(:)], 'LineWidth', 2);
    grid on;
    xlabel('Date');
    ylabel('Cases/day');

    title('France: held-out forecast');
    legend('Held-out reports', 'Known-policy forecast');
end
