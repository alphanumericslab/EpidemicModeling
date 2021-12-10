clear;
close all;
clc;


AllCases = importdata('./../../COVID-19/csse_covid_19_data/csse_covid_19_time_series/time_series_covid19_confirmed_global.csv');
AllDeaths = importdata('./../../COVID-19/csse_covid_19_data/csse_covid_19_time_series/time_series_covid19_deaths_global.csv');

CountryList = {'US'}%, 'Brazil', 'Russia', 'Spain', 'United Kingdom', 'Italy', 'France', 'Germany', 'Turkey', 'India', 'Iran', 'Peru', 'Canada', 'China'};
ttl = 'Global COVID-19 mortality rate per official case reports, by 26 May 2020';

% CountryList = {'Russia', 'Spain', 'United Kingdom', 'Italy', 'France', 'Germany', 'Belgium', 'Netherlands', 'Belarus', 'Sweden', 'Portugal', 'Switzerland'};
% ttl = 'European COVID-19 mortality rate per official case reports, by 26 May 2020';

style = {'-', '.', 's', '|', '.*'};
NumDays = size(AllCases.data, 2) - 2;
min_cases = 100; % min number of cases
period = 300; % days

% load totaltests tests

ratio = zeros(length(CountryList), NumDays);
total = zeros(length(CountryList), NumDays);
death = zeros(length(CountryList), NumDays);
min_cases_date = zeros(1, length(CountryList));
for k = 1 : length(CountryList)
    CountryRows = find(contains(AllCases.textdata(:, 2) , CountryList(k)));
    cases = sum(AllCases.data(CountryRows - 1, 3:end), 1);
    deaths = sum(AllDeaths.data(CountryRows - 1, 3:end), 1);
    
    %     if(k == 9)%isequal(CountryList(k) , 'Iran'))
    %         deaths = deaths + round([zeros(1, 61), 0.05*cumsum(ones(1, 19))].*deaths);
    %     end
    
    %     if(k == 6)%isequal(CountryList(k) , 'China'))
    %         deaths = deaths + round([zeros(1, 1), 0.05*cumsum(ones(1, 79))].*deaths);
    %     end
    
    %     if(k == 9)%isequal(CountryList(k) , 'Turkey'))
    %         deaths = deaths + round([zeros(1, 1), 0.05*cumsum(ones(1, 79))].*deaths);
    %     end
    
    ratio(k, :) = 100*deaths ./cases;% .* tests(k);
    total(k, :) = cases;
    death(k, :) = deaths;
    min_cases_date(k) = find(cases() >= min_cases, 1);
end
% ratio(ratio > 15.0) = nan;

for k = 1 : length(CountryList)
    if k == 1
        figure
        days = 0 : 60;
        semilogy(days, min_cases*2.^(days/(1.0)), 'linewidth', 2, 'color', 0.7*ones(1, 3));
        hold on
        semilogy(days, min_cases*2.^(days/(2.0)), 'linewidth', 2, 'color', 0.7*ones(1, 3));
        semilogy(days, min_cases*2.^(days/(3.0)), 'linewidth', 2, 'color', 0.7*ones(1, 3));
        semilogy(days, min_cases*2.^(days/(4.0)), 'linewidth', 2, 'color', 0.7*ones(1, 3));
        semilogy(days, min_cases*2.^(days/(5.0)), 'linewidth', 2, 'color', 0.7*ones(1, 3));
        semilogy(days, min_cases*2.^(days/(6.0)), 'linewidth', 2, 'color', 0.7*ones(1, 3));
    end
    indexes = min_cases_date(k) : min(min_cases_date(k) + period, size(total, 2));
    semilogy(0 : length(indexes) - 1, total(k, indexes)', 'linewidth', 2);
    text(length(indexes)-1, total(k, indexes(end)), CountryList{k});
end
a = axis;
a(3) = min_cases;
a(4) = 1e6;
axis(a);


grid
% legend(CountryList);
xlabel(['days since ' num2str(min_cases) ' case reports']);
ylabel('number of cases');
set(gca, 'fontsize', 14);
set(gca, 'box', 'on');

figure
for k = 1 : length(CountryList)
    indexes = min_cases_date(k) : min(min_cases_date(k) + period, size(total, 2));
    infected = diff(total(k, indexes));
    [mx, peak] = max(infected);
    plot(1 : length(indexes) - 1, infected, 'linewidth', 2);
    %     text(length(indexes)-1, infected(end), CountryList{k});
    text(peak, mx, CountryList{k});
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
title(ttl);

% for k = 1 : length(CountryList)
%     indexes = min_cases_date(k) : min(min_cases_date(k) + period, size(total, 2));
%     CountryRows = find(contains(AllCases.textdata(:, 2) , CountryList(k)));
%     cases = sum(AllCases.data(CountryRows - 1, 3:end), 1);
%     deaths = sum(AllDeaths.data(CountryRows - 1, 3:end), 1);
%     daily_cases = diff(cases);
%     daily_deaths = diff(deaths);
%     
%     lag = finddelay(daily_cases(indexes(1:end-1)), daily_deaths(indexes(1:end-1)));
%     
%     figure
%     hold on
%     plot(daily_cases, 'linewidth', 2);
%     plot(daily_deaths, 'linewidth', 2);
% %     plot(deaths, 'linewidth', 2);
%     grid
%     title([CountryList(k) 'delay = ' num2str(lag)]);
%     set(gca, 'fontsize', 14);
%     set(gca, 'box', 'on');
%     legend('I(t)', 'P(t)');
%     xlabel('day');
% end


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