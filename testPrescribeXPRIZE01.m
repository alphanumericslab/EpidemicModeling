% a script for testing the prescription algorithm developed for the XPRIZE
% Pandemic Response Challenge

clear;
close all;
clc

START_DATE = "2020-05-01"; % start time
END_DATE = "2021-02-02"; % end time
LATEST_DATA_FILE = './../covid-policy-tracker/data/OxCGRT_latest.csv'; % The historic data file cloned from: https://github.com/OxCGRT/covid-policy-tracker/tree/master/data
GEO_FILE = "xprize-sample-data/countries_regions.csv"; % countries and regions to include
POPULATION_FILE = "xprize-sample-data/populations.csv"; % country and regional populations

TRAINED_MODEL_PARAMS_FILE = "xprize-sample-data/prescription_trained_params.csv"; % file to log the trained model parameters
INCLUDED_IP = {'C1_School closing',... % see full descriptions at: https://github.com/OxCGRT/covid-policy-tracker/blob/master/documentation/codebook.md%'C1_Flag',...
    'C2_Workplace closing',... %'C2_Flag',...
    'C3_Cancel public events',... %'C3_Flag',...
    'C4_Restrictions on gatherings',... %'C4_Flag',...
    'C5_Close public transport',... %'C5_Flag',...
    'C6_Stay at home requirements',... %'C6_Flag',...
    'C7_Restrictions on internal movement',... 'C7_Flag',...
    'C8_International travel controls',... %'E1_Income support', 'E1_Flag', 'E2_Debt/contract relief',... % 'E3_Fiscal measures', 'E4_International support',...
    'H1_Public information campaigns',... %'H1_Flag',...
    'H2_Testing policy',...
    'H3_Contact tracing',... % 'H4_Emergency investment in healthcare', 'H5_Investment in vaccines',...
    'H6_Facial Coverings',... %'H6_Flag',... % {'H7_Vaccination policy', 'H7_Flag',...%'Holidays'...
    };
IP_MAXES = [3, 3, 2, 4, 2, 3, 2, 4, 2, 3, 2, 4]'; % Maximum values for the 12 phase NPI values defined in the oxford data set




% % % CHECK the infections vs the JH uni dataset







% Train the model parameters
TrainNPIPrescriptor(START_DATE, END_DATE, LATEST_DATA_FILE, GEO_FILE, POPULATION_FILE, INCLUDED_IP, IP_MAXES, TRAINED_MODEL_PARAMS_FILE);


IP_FILE = "prescriptions/robojudge_test_scenario.csv";
costs_file = "xprize-sample-data/uniform_random_costs.csv";
output_file = "xprize-sample-data/prescriptor2020-08-01_2020-08-04.csv";

ip_file_historic = "xprize-sample-data/2020-09-30_historical_ip.csv";


PrescribeNPI(START_DATE, END_DATE, ip_file, costs_file, output_file);

license('inuse') % list the used licenses (useful for building standalone packages)