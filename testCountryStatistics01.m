clear;
close all;
clc;


AllCases = importdata('Z:\csse_covid_19_data\csse_covid_19_time_series\time_series_covid19_confirmed_global.csv');
AllDeaths = importdata('Z:\csse_covid_19_data\csse_covid_19_time_series\time_series_covid19_deaths_global.csv');
% CountryList = {'Canada', 'US', 'China', 'Iran', 'France', 'Germany', 'Italy', 'Spain', 'Netherlands', 'Belgium', 'United Kingdom'};
CountryList = {'Korea', 'US', 'Spain', 'Italy', 'France', 'Germany', 'United Kingdom', 'China', 'Iran', 'Turkey', 'Belgium', 'Netherlands', 'Switzerland'};%,  'Canada', 'Brazil', 'Portugal'};
style = {'-', '.', 's', '|', '.*'};
NumDays = size(AllCases.data, 2) - 2;
min_cases = 100; % min number of cases
period = 55; % days

ratio = zeros(length(CountryList), NumDays);
total = zeros(length(CountryList), NumDays);
min_cases_date = zeros(1, length(CountryList));
for k = 1 : length(CountryList)
    CountryRows = find(contains(AllCases.textdata(:, 2) , CountryList(k)));
    cases = sum(AllCases.data(CountryRows - 1, 3:end), 1);
    deaths = sum(AllDeaths.data(CountryRows - 1, 3:end), 1);
    
    %     if(k == 7)%isequal(CountryList(k) , 'Iran'))
    %         deaths = deaths + round([zeros(1, 61), 0.05*cumsum(ones(1, 19))].*deaths);
    %     end
    
    %     if(k == 6)%isequal(CountryList(k) , 'China'))
    %         deaths = deaths + round([zeros(1, 1), 0.05*cumsum(ones(1, 79))].*deaths);
    %     end
    
    %     if(k == 9)%isequal(CountryList(k) , 'Turkey'))
    %         deaths = deaths + round([zeros(1, 1), 0.05*cumsum(ones(1, 79))].*deaths);
    %     end
    
    ratio(k, :) = 100*deaths ./cases;
    total(k, :) = cases;
    min_cases_date(k) = find(cases() >= min_cases, 1);
end
% ratio(ratio > 15.0) = nan;

figure
for k = 1 : length(CountryList)
    indexes = min_cases_date(k) : min(min_cases_date(k) + period, size(total, 2));
    semilogy(0 : length(indexes) - 1, total(k, indexes)', 'linewidth', 2);
    text(length(indexes)-1, total(k, indexes(end)), CountryList{k});
    if k == 1
        hold on
    end
end
grid
% legend(CountryList);
xlabel(['days since ' num2str(min_cases) ' case reports']);
ylabel('number of cases');
set(gca, 'fontsize', 14);
set(gca, 'box', 'on');

figure
for k = 1 : length(CountryList)
    indexes = min_cases_date(k) : min(min_cases_date(k) + period, size(total, 2));
    plot(1 : length(indexes) - 1, diff(log(total(k, indexes))), 'linewidth', 2);
    text(length(indexes)-1, total(k, indexes(end)), CountryList{k});
    if k == 1
        hold on
    end
end
grid
% legend(CountryList);
xlabel(['days since ' num2str(min_cases) ' case reports']);
ylabel('number of cases');
set(gca, 'fontsize', 14);
set(gca, 'box', 'on');

figure
hold on
for k = 1 : length(CountryList)
    indexes = min_cases_date(k) : min(min_cases_date(k) + period, size(total, 2));
    plot(0 : length(indexes) - 1, ratio(k, indexes)', 'linewidth', 2);
    text(length(indexes)-1, ratio(k, indexes(end)), CountryList{k});
    plot(0 : length(indexes) - 1, ratio(k, indexes)', '.');
end
grid
% legend(CountryList);
% xlabel(['days since ', AllCases.textdata(1, 5)]);
xlabel(['days since ' num2str(min_cases) ' case reports']);
ylabel('death to total cases ratio (%)');
set(gca, 'fontsize', 14);
set(gca, 'box', 'on');


% % % figure
% % % hold on
% % % % for k = 1 : length(CountryList)
% % % semilogy(total(:, 1:end)');
% % % % end
% % % grid
% % % legend(CountryList);
% % % xlabel(['days since ', AllCases.textdata(1, 5)]);
% % % ylabel('death to total cases %');
% % % set(gca, 'fontsize', 14);
% % % set(gca, 'box', 'on');