function [total, infected, recovered, deceased, first_index, threshold_index, num_days] = read_covid19_data(confirmed_datafile, death_datafile, recovered_datafile, region_list, min_cases)
% READ_COVID19_DATA Aggregate aligned Johns Hopkins wide-format country time series.
% CSVs have four metadata columns followed by aligned date columns. region_list
% contains country substrings; min_cases is the cumulative-count threshold.
% Outputs are regions-by-days. Index vectors are one-based, with 0 for a
% country that never reaches the first-case or threshold condition.
% Author: Reza Sameni | Emory University
paths = {confirmed_datafile,death_datafile,recovered_datafile}; aggregates = cell(1,3);
regions = string(region_list); reference_dates = {};
for j = 1:3
    table_data = readtable(paths{j},'VariableNamingRule','preserve','TextType','string');
    dates = table_data.Properties.VariableNames(5:end);
    if j == 1, reference_dates = dates; else, assert(isequal(dates,reference_dates), 'Date columns must align.'); end
    country = string(table_data{:,2}); country(ismissing(country)) = "";
    values = zeros(numel(regions),numel(dates));
    for k = 1:numel(regions), values(k,:) = sum(table_data{contains(country,regions(k)),5:end},1,'omitnan'); end
    aggregates{j} = values;
end
total = aggregates{1}; deceased = aggregates{2}; recovered = aggregates{3}; infected = total-deceased-recovered;
num_days = size(total,2); first_index = zeros(1,numel(regions)); threshold_index = first_index;
for k = 1:numel(regions)
    idx = find(total(k,:)>0,1); if ~isempty(idx), first_index(k) = idx; end
    idx = find(total(k,:)>=min_cases,1); if ~isempty(idx), threshold_index(k) = idx; end
end
end
