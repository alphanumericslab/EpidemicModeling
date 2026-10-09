function data = read_oxford_data(path)
    % READ_OXFORD_DATA Read, normalize and sort an offline Oxford/XPRIZE CSV.
    % data = read_oxford_data(path)
    % Required columns: CountryName, Date, ConfirmedCases. RegionName blanks
    % become empty strings; Date accepts YYYYMMDD or ISO dates. Duplicate daily
    % geography rows are rejected. Output is a table with preserved NPI names.
    % Author: Reza Sameni | Emory University
    opts = detectImportOptions(path, 'VariableNamingRule', 'preserve');
    text_columns = intersect({'CountryName', 'RegionName'}, opts.VariableNames);
    opts = setvartype(opts, text_columns, 'string');
    data = readtable(path, opts);

    assert(all(ismember({'CountryName', 'Date', 'ConfirmedCases'}, ...
        data.Properties.VariableNames)), 'Missing Oxford columns.');

    if ~ismember('RegionName', data.Properties.VariableNames)
        data.RegionName = repmat("", height(data), 1);
    end

    data.RegionName(ismissing(data.RegionName)) = "";

    if isnumeric(data.Date)
        data.Date = datetime(string(data.Date), 'InputFormat', 'yyyyMMdd');
    else
        date_strings = string(data.Date);

        if all(strlength(date_strings) == 8)
            format = 'yyyyMMdd';
        else
            format = 'yyyy-MM-dd';
        end

        data.Date = datetime(date_strings, 'InputFormat', format);
    end

    keys = strcat(data.CountryName, "|", data.RegionName, "|", string(data.Date, 'yyyy-MM-dd'));

    assert(numel(unique(keys)) == height(data), 'Duplicate country/region/date rows.');
    data = sortrows(data, {'CountryName', 'RegionName', 'Date'});
end
