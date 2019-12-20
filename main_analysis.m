%% Define global variabiles:
global N_BINS EARLIEST_VALID_PT LATEST_VALID_PT SUCCESS_TH_ANGLE RED_COLOR...
    GREEN_COLOR BLUE_COLOR HOME_DIR MIN_N_PT_FOR_MEASURE

HOME_DIR = pwd;
N_BINS = 6;
SUCCESS_TH_ANGLE = 30;
%% 08/25/2019 - correct for system delay:
% EARLIEST_VALID_PT = -.1;
% LATEST_VALID_PT = .7;
EARLIEST_VALID_PT = -.16;
LATEST_VALID_PT = .64;
MIN_N_PT_FOR_MEASURE = 4;

RED_COLOR = [172, 59, 59]/255;
GREEN_COLOR = [85, 170, 85]/255; 
BLUE_COLOR = [86 85 149]/255;

%% Analysis of Experiment 1:
analyze_retention_group_v1;
min_pt_E1 = min_pt';

%% Analysis of Experiment 2:
analyze_statLearn_group_v3
min_pt_E2 = min_pt_all;

%% Meta analysis:
meta_analysis_v1;

%% Run stats in R:
disp('Run R script');
pause;

%% continue with main analysis:
post_main;



