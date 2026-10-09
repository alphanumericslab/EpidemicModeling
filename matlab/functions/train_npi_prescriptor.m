function bundle = train_npi_prescriptor(start_date_str, end_date_str, data_file, geo_file, populations_file, included_ip, npi_maxes, trained_model_params_file, regression_start_date)
% TRAIN_NPI_PRESCRIPTOR Train selected Oxford geographies and write a JSON bundle.
% bundle = train_npi_prescriptor(start_date_str,end_date_str,data_file,geo_file,...
%     populations_file,included_ip,npi_maxes,trained_model_params_file)
% Dates are inclusive ISO strings. data_file is an Oxford CSV; geo_file lists
% CountryName/RegionName; populations_file requires Population2020. included_ip
% lists preserved NPI column names; npi_maxes supplies bounds. The optional
% regression_start_date restricts NNLS fitting; filtering uses all training days.
% Output JSON schema is shared with Python and contains models and npi_columns.
% Geographies with <14 observations are skipped; daily gaps cause an error.
% Author: Reza Sameni | Emory University
% Reference: Sameni (2022), doi:10.1109/JSTSP.2021.3129118.
if nargin < 9, regression_start_date = ''; end
data = read_oxford_data(data_file); geos = read_geo_table(geo_file); populations = read_geo_table(populations_file);
models = {}; included_ip = cellstr(included_ip);
for k = 1:height(geos)
    country = geos.CountryName(k); region = geos.RegionName(k);
    selected = data.CountryName == country & data.RegionName == region & data.Date >= datetime(start_date_str) & data.Date <= datetime(end_date_str);
    segment = data(selected,:); if height(segment) < 14, continue; end
    assert(all(diff(segment.Date) == days(1)), 'Training dates must be consecutive.');
    pop = populations(populations.CountryName == country & populations.RegionName == region,:);
    assert(height(pop) == 1, 'Need exactly one population per geography.');
    u = fillmissing(segment{:,included_ip},'previous'); u(isnan(u)) = 0; u = u';
    regression_start = 1;
    if ~isempty(regression_start_date)
        regression_start = find(segment.Date >= datetime(regression_start_date),1);
        assert(~isempty(regression_start), 'Regression start is after training data.');
    end
    model = fit_npi_model(segment.ConfirmedCases,u,pop.Population2020,npi_maxes,regression_start);
    model.country_name = char(country); model.region_name = char(region);
    model.last_date = char(string(segment.Date(end),'yyyy-MM-dd')); model.last_input = u(:,end);
    models{end+1} = model; %#ok<AGROW>
end
assert(~isempty(models), 'No selected geography has at least 14 training days.');
bundle.schema_version = 1; bundle.npi_columns = included_ip(:)'; bundle.models = models;
fid = fopen(trained_model_params_file,'w'); assert(fid >= 0, 'Cannot open model output.');
cleanup = onCleanup(@() fclose(fid)); fprintf(fid,'%s\n',jsonencode(bundle,'PrettyPrint',true));
end
