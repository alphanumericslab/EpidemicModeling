function output = prescribe_npi(start_date_str, end_date_str, ip_file, costs_file, output_file, ...
    model_file, epsilon)
    % PRESCRIBE_NPI Write XPRIZE-format daily NPI prescriptions from a JSON model.
    % output = prescribe_npi(start_date_str,end_date_str,ip_file,costs_file,...
    %     output_file,model_file,epsilon)
    % ip_file selects CountryName/RegionName pairs. costs_file holds per-NPI costs
    % in the model's original column names. model_file is explicit (required);
    % epsilon is the human/NPI mix (default .5). Returns the written table.
    % Prescription must start the day after each model's training endpoint.
    % Python uses keyword-only model_file and epsilon. No hardcoded paths are used.
    % Author: Reza Sameni | Emory University

    if nargin < 7
        epsilon = .5;
    end

    bundle = jsondecode(fileread(model_file));
    costs = read_geo_table(costs_file);
    plans = read_geo_table(ip_file);
    dates = (datetime(start_date_str):days(1):datetime(end_date_str))';

    assert(~isempty(dates), 'Invalid date interval.');
    output = table();
    models = bundle.models;

    if ~iscell(models)
        models = num2cell(models);
    end

    columns = cellstr(bundle.npi_columns);

    for k = 1:numel(models)
        model = models{k};
        country = string(model.country_name);
        region = string(model.region_name);

        if ~any(plans.CountryName == country & plans.RegionName == region)
            continue
        end

        assert(datetime(model.last_date) + days(1) == dates(1), ...
            'Prescription must follow training immediately.');
        selected = costs(costs.CountryName == country & costs.RegionName == region, :);

        assert(height(selected) == 1, 'Missing or duplicate intervention costs.');
        weights = selected{1, columns}';
        [controls, ~, ~] = optimal_npi(model, numel(dates), weights, epsilon);
        block = table(repmat(country, numel(dates), 1), repmat(region, numel(dates), 1), ...
            string(dates, 'yyyy-MM-dd'), 'VariableNames', {'CountryName', 'RegionName', 'Date'});
        block = [block array2table(controls', 'VariableNames', columns)];

        output = [output; block]; %#ok<AGROW>
    end

    assert(height(output) > 0, 'No models matched the intervention plan.');
    writetable(output, output_file);
end
