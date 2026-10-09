function data = read_geo_table(path)
% READ_GEO_TABLE Read a geographic CSV with preserved names and normalized blanks.
% data = read_geo_table(path). CountryName and RegionName are string columns;
% missing regions become empty strings. Other columns retain inferred types.
% Author: Reza Sameni | Emory University
opts = detectImportOptions(path,'VariableNamingRule','preserve');
opts = setvartype(opts,intersect({'CountryName','RegionName'},opts.VariableNames),'string');
data = readtable(path,opts);
if ~ismember('RegionName',data.Properties.VariableNames), data.RegionName = repmat("",height(data),1); end
data.RegionName(ismissing(data.RegionName)) = "";
end
