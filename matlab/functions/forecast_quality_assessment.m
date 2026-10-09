function output = forecast_quality_assessment(npi_weights, human_npi_cost_factor, ...
    start_train_date_str, end_train_date_str, start_regression_date_str, ...
    end_predict_prescribe_date_str, max_look_ahead_days, data_file, geo_file, populations_file, ...
    included_ip, npi_mins, npi_maxes, trained_model_params_file)
    % FORECAST_QUALITY_ASSESSMENT Score held-out daily forecasts without future-case leakage.
    % Arguments follow train_predict_prescribe_npi plus max_look_ahead_days.
    % Known held-out NPI schedules are scenario inputs; held-out cases are used
    % only for scoring. Returns country_name, region_name, days, mae, rmse, bias.
    % npi_weights, human_npi_cost_factor and npi_mins retain the historical call
    % signature; they do not affect evaluation under known intervention schedules.
    % Author: Reza Sameni | Emory University
    bundle = train_npi_prescriptor(start_train_date_str, end_train_date_str, data_file, ...
        geo_file, populations_file, included_ip, npi_maxes, trained_model_params_file, ...
        start_regression_date_str);
    data = read_oxford_data(data_file);
    output = table();
    included_ip = cellstr(included_ip);
    endpoint = min(datetime(end_predict_prescribe_date_str), datetime(end_train_date_str) + ...
        days(max_look_ahead_days));

    for k = 1:numel(bundle.models)
        model = bundle.models{k};
        geo = data.CountryName == string(model.country_name) & data.RegionName == ...
            string(model.region_name);
        segment = data(geo & data.Date > datetime(end_train_date_str) & data.Date <= endpoint, :);

        if isempty(segment)
            continue
        end

        assert(segment.Date(1) == datetime(end_train_date_str) + days(1) && ...
            all(diff(segment.Date) == days(1)), 'Held-out dates must be consecutive.');
        last = data.ConfirmedCases(geo & data.Date == datetime(end_train_date_str));
        actual = max(0, diff([last(end); segment.ConfirmedCases]));
        u = fillmissing(segment{:, included_ip}, 'previous');
        u(isnan(u)) = 0;

        [predicted, ~] = forecast_npi(model, u');
        error = predicted(:) - actual;
        row = table(string(model.country_name), string(model.region_name), numel(error), ...
            mean(abs(error)), sqrt(mean(error.^2)), mean(error), 'VariableNames', ...
            {'country_name', 'region_name', 'days', 'mae', 'rmse', 'bias'});
        output = [output; row]; %#ok<AGROW>
    end
end
