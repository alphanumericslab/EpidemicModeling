clear;
close all;
clc;

AllCasesFname = './../COVID-19/csse_covid_19_data/csse_covid_19_time_series/time_series_covid19_confirmed_global.csv';
AllDeathsFname = './../COVID-19/csse_covid_19_data/csse_covid_19_time_series/time_series_covid19_deaths_global.csv';
AllRecoveredFname = './../COVID-19/csse_covid_19_data/csse_covid_19_time_series/time_series_covid19_recovered_global.csv';

RegionList = {'US', 'Brazil', 'Russia', 'Spain', 'United Kingdom', 'Italy', 'France', 'Germany', 'Turkey', 'India', 'Iran', 'Peru', 'Canada', 'China'};
% RegionList = {'Russia', 'Spain', 'United Kingdom', 'Italy', 'France', 'Germany', 'Belgium', 'Netherlands', 'Belarus', 'Sweden', 'Portugal', 'Switzerland'};
% RegionList = {'Spain', 'United Kingdom', 'Italy', 'France', 'Germany'};
% RegionList = {'France', 'Iran'};

style = {'-', '.', 's', '|', '.*'};
min_cases = 100; % min number of cases
period = 300; % days
avgdays = 3;
R0period = 7;

% load COBID-19 data (https://github.com/CSSEGISandData/COVID-19.git)
[TotalCases, Infected, Recovered, Deceased, FirstCaseDateIndex, MinCaseDateIndex, NumDays] = ReadCOVID19Data(AllCasesFname, AllDeathsFname, AllRecoveredFname, RegionList, min_cases);

wlen = 7;
NewCases = diff(TotalCases, 1, 2);
NewCasesSmoothed = filter(ones(1, wlen), wlen, NewCases')';
NewCasesSmoothed(NewCasesSmoothed < 0) = 0;

InfectedSmoothed = filter(ones(1, wlen), wlen, Infected')';
InfectedSmoothed(InfectedSmoothed < 0) = 0;

lag = 0;
% Rate = NewCasesSmoothed(:, 1 + lag : end)./InfectedSmoothed(:, 1 : end-1-lag);
Rate = NewCases(:, 1 + lag : end)./Infected(:, 1 : end-1-lag);
Rate(Rate < 0) = 0;

figure;
% plot(NewCases');
plot(NewCasesSmoothed');
grid
title('NewCasesSmoothed');

figure;
% plot(NewCases');
plot(Infected');
grid
title('Infected');

country_code = 4;
lgn = {};
figure
hold on
plot(TotalCases(country_code, :)); lgn = cat(2, lgn, 'Total');
plot(Infected(country_code, :)); lgn = cat(2, lgn, 'Active Cases');
stem(NewCases(country_code, :)); lgn = cat(2, lgn, 'New cases');
plot(diff(Recovered(country_code, :))); lgn = cat(2, lgn, 'Recovered');
plot(diff(Deceased(country_code, :))); lgn = cat(2, lgn, 'Deceased');
grid
legend(lgn);
title(RegionList(country_code));

figure
plot(log(Rate(country_code, :)'));
grid
title(['log(Rate) for ' RegionList(country_code)]);

% % % % MortalityRate = 100*Deceased./TotalCases;% .* tests(k);
% % % % 
% % % % % Smooth the data using Tikhonov regularization
% % % % % DiffOrderOrFilterCoefs = [1 -2 1];
% % % % DiffOrderOrFilterCoefs = 2;
% % % % lambda = 3.0;
% % % % InfectedSmoothed = TikhonovRegularization(Infected, DiffOrderOrFilterCoefs, lambda);
% % % % % InfectedSmoothed = LPFilter(Infected, 0.1);
% % % % 
% % % % ratio = zeros(length(RegionList), NumDays);
% % % % R0 = zeros(length(RegionList), NumDays);
% % % % for k = 1 : length(RegionList)
% % % %     % ratio(k, FirstCaseDateIndex(k)+1:end) = InfectedSmoothed(k, FirstCaseDateIndex(k)+1:end)./InfectedSmoothed(k, FirstCaseDateIndex(k):end-1);
% % % %     % ratio(k, FirstCaseDateIndex(k)+1:end) = Infected(k, FirstCaseDateIndex(k)+1:end)./InfectedSmoothed(k, FirstCaseDateIndex(k):end-1);
% % % %     ratio(k, FirstCaseDateIndex(k) + R0period : end) = Infected(k, FirstCaseDateIndex(k) + R0period : end)./Infected(k, FirstCaseDateIndex(k) : end - R0period);
% % % %     
% % % %     %     R0(k, FirstCaseDateIndex(k)+1:end) = log(ratio(k, FirstCaseDateIndex(k)+1:end));
% % % %     R0(k, FirstCaseDateIndex(k)+1:end) = exp(filter(ones(1, avgdays), avgdays, log(ratio(k, FirstCaseDateIndex(k)+1:end))));
% % % % end
% % % % 
% % % % % Plot infected and smoothed infeced results
% % % % for k = 1 : length(RegionList)
% % % %     figure
% % % %     hold on
% % % %     plot(Infected(k, :), 'linewidth', 2);
% % % %     plot(InfectedSmoothed(k, :), 'linewidth', 2);
% % % %     legend('Infected', 'Filtered');
% % % %     title(['Active (infected) cases of ' RegionList{k}]);
% % % %     xlabel('days');
% % % %     ylabel('Number');
% % % %     grid;
% % % % end
% % % % 
% % % % 
% % % % days = 0 : 60;
% % % % figure
% % % % hold on
% % % % for k = 1 : length(RegionList)
% % % %     indexes = MinCaseDateIndex(k) : min(MinCaseDateIndex(k) + period, size(R0, 2));
% % % %     plot(0 : length(indexes) - 1, Infected(k, indexes), 'linewidth', 2);
% % % %     text(length(indexes)-1, Infected(k, indexes(end)), RegionList{k});
% % % % end
% % % % grid
% % % % title(['Number of active cases (infected) since the ' num2str(min_cases) ' case report']);
% % % % xlabel('days');
% % % % ylabel('number of cases');
% % % % set(gca, 'fontsize', 12);
% % % % 
% % % % figure
% % % % hold on
% % % % for k = 1 : length(RegionList)
% % % %     indexes = MinCaseDateIndex(k) : min(MinCaseDateIndex(k) + period, size(R0, 2));
% % % %     plot(0 : length(indexes) - 1, R0(k, indexes), 'linewidth', 2);
% % % %     if(R0(k, indexes(end)) > 1)
% % % %         clr = zeros(1, 3);
% % % %     else
% % % %         clr = 0.5 + zeros(1, 3);
% % % %     end
% % % %     text(length(indexes)-1, R0(k, indexes(end)), RegionList{k}, 'color', clr);
% % % % end
% % % % maxdays = 300;
% % % % plot((0:maxdays-1), ones(1, maxdays), 'k--', 'linewidth', 2);
% % % % grid
% % % % title(['Reproduction number for ' num2str(R0period) ' day generations geometrically averaged over ' num2str(avgdays) ' successive days']);
% % % % xlabel(['days since the ' num2str(min_cases) 'th case report']);
% % % % ylabel('R_0');
% % % % set(gca, 'fontsize', 12);
% % % % 
% % % % figure
% % % % hold on
% % % % for k = 1 : length(RegionList)
% % % %     indexes = MinCaseDateIndex(k) : min(MinCaseDateIndex(k) + period, size(MortalityRate, 2));
% % % %     plot(0 : length(indexes) - 1, MortalityRate(k, indexes)', 'linewidth', 2);
% % % %     text(length(indexes)-1, MortalityRate(k, indexes(end)), RegionList{k});
% % % %     plot(0 : length(indexes) - 1, MortalityRate(k, indexes)', '.');
% % % % end
% % % % grid
% % % % xlabel(['days since ' num2str(min_cases) 'th case reports']);
% % % % ylabel('death to total cases ratio (%)');
% % % % set(gca, 'fontsize', 12);
% % % % set(gca, 'box', 'on');
% % % % title('Mortality rate per official case reports');
% % % % 
