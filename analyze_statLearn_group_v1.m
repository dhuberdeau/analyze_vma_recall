function [data_grp, err_grp] = analyze_statLearn_group_v1(subject_list, N_TARGETS, varargin)

%% define subject list:
% subject_list = {...
% 'Data_SP026_03222018_E2.mat',...
% 'Data_S048_04052018_E2.mat',...
% ...% 'Data_S049_04052018_E2.mat'... %poor TR compliance
% 'Data_S050_04052018_E2.mat',...
% 'Data_S051_04052018_E2.mat',...
% 'Data_S052_04052018_E2.mat',...
% 'Data_S053_04052018_E2.mat',...
% ...% 'Data_S054_04052018_E2.mat',... %poor TR compliance
% ...% 'Data_S055_04052018_E2.mat',... %poor TR compliance
% 'Data_S056_04062018_E2.mat',...
% 'Data_S057_04062018_E2.mat',...
% 'Data_S058_04062018_E2.mat',...
% 'Data_S059_04062018_E2.mat',...
% 'Data_S060_04062018_E2.mat',...
% 'Data_S061_04062018_E2.mat',...
% 'Data_S062_04062018_E2.mat',...
% 'Data_S063_04062018_E2.mat',...    
% 'Data_S064_04062018_E2.mat',...
% 'Data_S065_04062018_E2.mat',...
% 'Data_S066_04092018_E2.mat',...
% 'Data_S069_04092018_E2_fixed.mat',... 
% ...% 'Data_S070_04092018_E2.mat',... %poor TR compliance
% 'Data_S079_04112018_E2.mat',...
%     };

if nargin < 3
    plot_bool = 1;
else
    plot_bool = varargin{1};
end

load('individual_memory_test_results.mat');

switch N_TARGETS
    case 3
        N_SYMBOLS = N_TARGETS*2;
    case 4
        N_SYMBOLS = N_TARGETS*3;
    case 6
        N_SYMBOLS = N_TARGETS*2;
    otherwise
        error('Invalid number of targets specified.');
end

succ_block_1 = nan(4, length(subject_list));
succ_block_2 = nan(4, length(subject_list));
succ_block_1_supTH = nan(4, length(subject_list));
succ_block_2_supTH = nan(4, length(subject_list));

succ_block_0_res = nan(8, length(subject_list));
succ_block_1_res = nan(8, length(subject_list));
succ_block_2_res = nan(8, length(subject_list));
succ_block_0_supTH_res = nan(8, length(subject_list));
succ_block_1_supTH_res = nan(8, length(subject_list));
succ_block_2_supTH_res = nan(8, length(subject_list));

succ_block_1_targ = nan(N_TARGETS, 4, length(subject_list)); % separate per target
succ_block_2_targ = nan(N_TARGETS, 4, length(subject_list)); % separate per target
succ_block_1_supTH_targ = nan(N_TARGETS, 4, length(subject_list)); % separate per target
succ_block_2_supTH_targ = nan(N_TARGETS, 4, length(subject_list)); % separate per target

succ_block_1_symb = nan(N_SYMBOLS, 2, length(subject_list)); % separate per symbol
succ_block_2_symb = nan(N_SYMBOLS, 2, length(subject_list)); % separate per symbol
succ_block_1_supTH_symb = nan(N_SYMBOLS, 2, length(subject_list)); % separate per symbol
succ_block_2_supTH_symb = nan(N_SYMBOLS, 2, length(subject_list)); % separate per symbol

succ_block_1_symb_marginal = nan(N_SYMBOLS, length(subject_list));
succ_block_2_symb_marginal = nan(N_SYMBOLS, length(subject_list));
succ_block_1_supTH_symb_marginal = nan(N_SYMBOLS, length(subject_list));
succ_block_2_supTH_symb_marginal = nan(N_SYMBOLS, length(subject_list));

min_pt_all = nan(1, length(subject_list));
viewtime_all = nan(400, length(subject_list));
de_all = nan(400, length(subject_list));
type_all = nan(400, length(subject_list));
% conf_mat = nan(N_TARGETS,N_TARGETS,length(subject_list));
move_class_all = nan(400, length(subject_list));
targ_all = nan(400, length(subject_list));

h1 = figure;
individ_analysis_errors = {};
for i_sub = 1:length(subject_list)
    try
    load(subject_list{i_sub});
    
    [data_indiv, err_out] = analyze_statLearn_individual_v1(Data, N_TARGETS, 0);
    
    targ_temp = Data.Target;
    targ_temp(targ_temp == 4) = 1;
    targ_temp(targ_temp == 5) = 2;
    targ_temp(targ_temp == 6) = 3;
%     cm_2 = confusionmat(targ_temp(Data.Type == 2), data_indiv.Mov_class(Data.Type == 2));
%     conf_mat(:, :, i_sub) = cm_2;
    
    p_bin0 = data_indiv.p_bins{1};
    p_bin1 = data_indiv.p_bins{2};
    p_bin2 = data_indiv.p_bins{3};
    pt_z_bin0 = data_indiv.z_bins{1};
    pt_z_bin2 = data_indiv.z_bins{2};
    
    pr_corr_block_0 = data_indiv.prCorr_lowRes{1};
    pr_corr_block_1 = data_indiv.prCorr_lowRes{2};
    pr_corr_block_2 = data_indiv.prCorr_lowRes{3};
    pr_corr_supTH_block_0 = data_indiv.prCorr_supTH_lowRes{1};
    pr_corr_supTH_block_1 = data_indiv.prCorr_supTH_lowRes{2};
    pr_corr_supTH_block_2 = data_indiv.prCorr_supTH_lowRes{3};
    
    pr_corr_block_0_res = data_indiv.prCorr_hiRes{1};
    pr_corr_block_1_res = data_indiv.prCorr_hiRes{2};
    pr_corr_block_2_res = data_indiv.prCorr_hiRes{3};
    pr_corr_supTH_block_0_res = data_indiv.prCorr_supTH_hiRes{1};
    pr_corr_supTH_block_1_res = data_indiv.prCorr_supTH_hiRes{2};
    pr_corr_supTH_block_2_res = data_indiv.prCorr_supTH_hiRes{3};
    
    pr_corr_block_1_target = data_indiv.prCorr_lowRes_target{1};
    pr_corr_block_2_target = data_indiv.prCorr_lowRes_target{2};
    pr_corr_supTH_block_1_target = data_indiv.prCorr_supTH_lowRes_target{1};
    pr_corr_supTH_block_2_target = data_indiv.prCorr_supTH_lowRes_target{2};

    pr_corr_block_1_symbol = data_indiv.prCorr_lowRes_symbol{1};
    pr_corr_block_2_symbol = data_indiv.prCorr_lowRes_symbol{2};
    pr_corr_supTH_block_1_symbol = data_indiv.prCorr_supTH_lowRes_symbol{1};
    pr_corr_supTH_block_2_symbol = data_indiv.prCorr_supTH_lowRes_symbol{2};
    
    pr_corr_block_1_symbol_marginal = data_indiv.prCorr_symbol_marginal{1};
    pr_corr_block_2_symbol_marginal = data_indiv.prCorr_symbol_marginal{2};
    pr_corr_block_1_supTH_symbol_marginal = data_indiv.prCorr_supTH_symbol_marginal{1};
    pr_corr_block_2_supTH_symbol_marginal = data_indiv.prCorr_supTH_symbol_marginal{2};
    
    min_pt_all(i_sub) = data_indiv.min_pt;

    type0 = data_indiv.types{1};
    type1 = data_indiv.types{2};
    type2 = data_indiv.types{3};
    type3 = data_indiv.types{4};
    type4 = data_indiv.types{5};

    Dir_e = data_indiv.Dir_e;
    Dir_a = data_indiv.Dir_a;

    Kin_x = data_indiv.kinematics{1};
    Kin_y = data_indiv.kinematics{2};

    viewtime_all(1:length(data_indiv.ViewTime), i_sub) = data_indiv.ViewTime;
    de_all(1:length(Dir_e), i_sub) = Dir_e;
    type_all(1:length(Data.Type), i_sub) = Data.Type;
    targ_all(1:length(Data.Target), i_sub) = Data.Target;
    move_class_all(1:length(Data.Target), i_sub) = data_indiv.Mov_class;
    
%     close all;
    succ_block_1(:, i_sub) = pr_corr_block_1';
    succ_block_2(:, i_sub) = pr_corr_block_2';
    
    succ_block_1_supTH(:, i_sub) = pr_corr_supTH_block_1;
    succ_block_2_supTH(:, i_sub) = pr_corr_supTH_block_2;
    
    succ_block_0_res(:, i_sub) = pr_corr_block_0_res';
    succ_block_1_res(:, i_sub) = pr_corr_block_1_res';
    succ_block_2_res(:, i_sub) = pr_corr_block_2_res';
    
    succ_block_0_supTH_res(:, i_sub) = pr_corr_supTH_block_0_res;
    succ_block_1_supTH_res(:, i_sub) = pr_corr_supTH_block_1_res;
    succ_block_2_supTH_res(:, i_sub) = pr_corr_supTH_block_2_res;
    
    succ_block_1_targ(:, :, i_sub) = pr_corr_block_1_target;
    succ_block_2_targ(:, :, i_sub) = pr_corr_block_2_target;
    succ_block_1_supTH_targ(:, :, i_sub) = pr_corr_supTH_block_1_target;
    succ_block_2_supTH_targ(:, :, i_sub) = pr_corr_supTH_block_2_target;
    
    succ_block_1_symb(:, :, i_sub) = pr_corr_block_1_symbol;
    succ_block_2_symb(:, :, i_sub) = pr_corr_block_2_symbol;
    succ_block_1_supTH_symb(:, :, i_sub) = pr_corr_supTH_block_1_symbol;
    succ_block_2_supTH_symb(:, :, i_sub) = pr_corr_supTH_block_2_symbol;
    
    succ_block_1_symb_marginal(:, i_sub) = pr_corr_block_1_symbol_marginal;
    succ_block_2_symb_marginal(:, i_sub) = pr_corr_block_2_symbol_marginal;
    succ_block_1_supTH_symb_marginal(:, i_sub) = pr_corr_block_1_supTH_symbol_marginal;
    succ_block_2_supTH_symb_marginal(:, i_sub) = pr_corr_block_2_supTH_symbol_marginal;
    
    open_winds = findobj('type', 'figure');
    for i_wind = 1:(length(open_winds)-1)
        close(i_wind);
    end
    figure(h1);
    subplot(4,5,i_sub); hold on
    plot(Kin_x(:, type0), Kin_y(:, type0), 'r.-');
    plot(Kin_x(:, type1), Kin_y(:, type1), 'g.-');
    plot(Kin_x(:, type2), Kin_y(:, type2), 'b.-');
    
    
    catch indiv_err
        individ_analysis_errors{length(individ_analysis_errors) + 1} = ...
            indiv_err;
       warning('Subject analysis failed') 
    end

end

%%
EARLIEST_VALID_PT = -.2;
LATEST_VALID_PT = .85;
if plot_bool
    figure;
    for i_sub = 1:length(subject_list)
        subplot(5,5,i_sub); hold on;
        plot(viewtime_all(type_all(:, i_sub) == 0, i_sub), de_all(type_all(:, i_sub) == 0, i_sub), '.', 'MarkerSize', 16, 'Color', [172, 59, 59]/255)
        plot(viewtime_all(type_all(:, i_sub) == 1, i_sub), de_all(type_all(:, i_sub) == 1, i_sub), '.', 'MarkerSize', 16, 'Color', [85, 170, 85]/255)
        plot(viewtime_all(type_all(:, i_sub) == 2, i_sub), de_all(type_all(:, i_sub) == 2, i_sub), '.', 'MarkerSize', 16, 'Color', [86/255 85/255 149/255])
        plot([min_pt_all(i_sub), min_pt_all(i_sub)], [-200 200], 'k-')
        axis([EARLIEST_VALID_PT, LATEST_VALID_PT, -200 200])
    end
end

%%
valid_subs = min_pt_all > .15;
% valid_subs([4, 6, 7, 21, 22]) = 1; % include subs who need modified mPT
% valid_subs(10) = 0; % exclude subs who didn't do task

succ_block_1_valid = succ_block_1(:,valid_subs);
succ_block_2_valid = succ_block_2(:,valid_subs);
if plot_bool
    figure; hold on;
    errorbar(1:4, nanmean(succ_block_1_valid'), sqrt(nanvar(succ_block_1_valid')./length(subject_list)), 'go', 'MarkerSize', 12, 'LineWidth', 2);
    errorbar(1:4, nanmean(succ_block_2_valid'), sqrt(nanvar(succ_block_2_valid')./length(subject_list)), 'bo', 'MarkerSize', 12, 'LineWidth', 2);
    axis([1 4 .4 1])
end
valid_subs = min_pt_all > .15;
succ_supTH_block_1_valid = succ_block_1_supTH(:,valid_subs);
succ_supTH_block_2_valid = succ_block_2_supTH(:,valid_subs);
if plot_bool
    figure; hold on;
    errorbar(1:4, nanmean(succ_supTH_block_1_valid'), sqrt(nanvar(succ_supTH_block_1_valid')./length(subject_list)), 'go', 'MarkerSize', 12, 'LineWidth', 2);
    errorbar(1:4, nanmean(succ_supTH_block_2_valid'), sqrt(nanvar(succ_supTH_block_2_valid')./length(subject_list)), 'bo', 'MarkerSize', 12, 'LineWidth', 2);
    axis([1 4 .4 1])

    [a,b,c] = anova1(succ_block_2_valid');
end
[x,y,z,q] = ttest(succ_block_2_valid(1,:) - succ_block_2_valid(end,:));

%% same but with greater resolution
valid_subs = min_pt_all > .15;
% valid_subs([4, 6, 7, 21, 22]) = 1; % include subs who need modified mPT
% valid_subs(10) = 0; % exclude subs who didn't do task
succ_block_0_valid = succ_block_0_res(:,valid_subs);
succ_block_1_valid = succ_block_1_res(:,valid_subs);
succ_block_2_valid = succ_block_2_res(:,valid_subs);
if plot_bool
    figure; hold on;
    errorbar(.5:.5:4, nanmean(succ_block_0_valid'), sqrt(nanvar(succ_block_0_valid')./length(subject_list)), 'r-o', 'MarkerSize', 12, 'LineWidth', 2);
    errorbar(.5:.5:4, nanmean(succ_block_1_valid'), sqrt(nanvar(succ_block_1_valid')./length(subject_list)), 'g-o', 'MarkerSize', 12, 'LineWidth', 2);
    errorbar(.5:.5:4, nanmean(succ_block_2_valid'), sqrt(nanvar(succ_block_2_valid')./length(subject_list)), 'b-o', 'MarkerSize', 12, 'LineWidth', 2);
    axis([0 4.5 .4 1]);
end

valid_subs = min_pt_all > .15;
succ_supTH_block_0_valid = succ_block_0_supTH_res(:,valid_subs);
succ_supTH_block_1_valid = succ_block_1_supTH_res(:,valid_subs);
succ_supTH_block_2_valid = succ_block_2_supTH_res(:,valid_subs);
if plot_bool
    figure; hold on;
    errorbar(.5:.5:4, nanmean(succ_supTH_block_0_valid'), sqrt(nanvar(succ_supTH_block_0_valid')./length(subject_list)), 'r-o', 'MarkerSize', 12, 'LineWidth', 2);
    errorbar(.5:.5:4, nanmean(succ_supTH_block_1_valid'), sqrt(nanvar(succ_supTH_block_1_valid')./length(subject_list)), 'g-o', 'MarkerSize', 12, 'LineWidth', 2);
    errorbar(.5:.5:4, nanmean(succ_supTH_block_2_valid'), sqrt(nanvar(succ_supTH_block_2_valid')./length(subject_list)), 'b-o', 'MarkerSize', 12, 'LineWidth', 2);
    axis([0 4.5 .4 1]);

    [a,b,c] = anova1(succ_block_2_valid');
end
[x,y,z,q] = ttest(succ_block_2_valid(1,:) - succ_block_2_valid(end,:));

%% same but with each target individually
valid_subs = min_pt_all > .15;
% valid_subs([4, 6, 7, 21, 22]) = 1; % include subs who need modified mPT
% valid_subs(10) = 0; % exclude subs who didn't do task
succ_block_1_valid = succ_block_1_targ(:,:,valid_subs);
succ_block_2_valid = succ_block_2_targ(:,:,valid_subs);

rate_proxy_block_1 = reshape(nanmean(succ_block_1_valid, 2), N_TARGETS, sum(valid_subs));
rate_proxy_block_2 = reshape(nanmean(succ_block_2_valid, 2), N_TARGETS, sum(valid_subs));
[~, sort_ind_1] = sort(rate_proxy_block_1, 1);
[~, sort_ind_2] = sort(rate_proxy_block_2, 1);

succ_block_1_ordered = nan(size(succ_block_1_valid));
succ_block_2_ordered = nan(size(succ_block_2_valid));
for i_sub = 1:size(succ_block_1_valid,3)
    succ_block_1_ordered(:, :, i_sub) = succ_block_1_valid(sort_ind_1(:, i_sub), :, i_sub);
    succ_block_2_ordered(:, :, i_sub) = succ_block_2_valid(sort_ind_2(:, i_sub), :, i_sub);
end

symbol_list_targ = {'-o', '-s', '-x', '-.', '--', '.-'}; 
if plot_bool
    figure; hold on;
    for i_targ = 1:N_TARGETS
        errorbar(1:4, nanmean(succ_block_1_valid(i_targ, :, :), 3),...
            sqrt(nanvar(succ_block_1_valid(i_targ, :, :),[],3)./length(subject_list)),...
            ['g', symbol_list_targ{i_targ}], 'MarkerSize', 12, 'LineWidth', 2);

        errorbar(1:4, nanmean(succ_block_2_valid(i_targ, :, :), 3),...
            sqrt(nanvar(succ_block_2_valid(i_targ, :, :),[],3)./length(subject_list)),...
            ['b', symbol_list_targ{i_targ}], 'MarkerSize', 12, 'LineWidth', 2);
    end
    axis([1 4 .4 1])
end

% now for super-threshold PT
valid_subs = min_pt_all > .15;
succ_supTH_block_1_valid = succ_block_1_supTH_targ(:,:,valid_subs);
succ_supTH_block_2_valid = succ_block_2_supTH_targ(:,:,valid_subs);

rate_proxy_block_1 = reshape(nanmean(succ_supTH_block_1_valid, 2), N_TARGETS, sum(valid_subs));
rate_proxy_block_2 = reshape(nanmean(succ_supTH_block_2_valid, 2), N_TARGETS, sum(valid_subs));
[~, sort_ind_1] = sort(rate_proxy_block_1, 1);
[~, sort_ind_2] = sort(rate_proxy_block_2, 1);

succ_supTH_block_1_ordered = nan(size(succ_supTH_block_1_valid));
succ_supTH_block_2_ordered = nan(size(succ_supTH_block_2_valid));
for i_sub = 1:size(succ_supTH_block_1_valid,3)
    succ_supTH_block_1_ordered(:, :, i_sub) = succ_supTH_block_1_valid(sort_ind_1(:, i_sub), :, i_sub);
    succ_supTH_block_2_ordered(:, :, i_sub) = succ_supTH_block_2_valid(sort_ind_2(:, i_sub), :, i_sub);
end

if plot_bool
    figure; hold on;
    for i_targ = 1:N_TARGETS
        errorbar(1:4, nanmean(succ_supTH_block_1_ordered(i_targ, :, :), 3),...
            sqrt(nanvar(succ_supTH_block_1_ordered(i_targ, :, :),[],3)./length(subject_list)),...
            ['g', symbol_list_targ{i_targ}], 'MarkerSize', 12, 'LineWidth', 2);

        errorbar(1:4, nanmean(succ_supTH_block_2_ordered(i_targ, :, :), 3),...
            sqrt(nanvar(succ_supTH_block_2_ordered(i_targ, :, :),[],3)./length(subject_list)),...
            ['b', symbol_list_targ{i_targ}], 'MarkerSize', 12, 'LineWidth', 2);
    end
    axis([1 4 .4 1])
end
% [a,b,c] = anova1(succ_block_2_valid');
% [x,y,z,q] = ttest(succ_block_2_valid(1,:) - succ_block_2_valid(end,:));

%% same but with each symbol individually
valid_subs = min_pt_all > .15;
% valid_subs([4, 6, 7, 21, 22]) = 1; % include subs who need modified mPT
% valid_subs(10) = 0; % exclude subs who didn't do task
succ_block_1_valid = succ_block_1_symb(:,:,valid_subs);
succ_block_2_valid = succ_block_2_symb(:,:,valid_subs);

rate_proxy_block_1 = reshape(nanmean(succ_block_1_valid, 2), N_SYMBOLS, sum(valid_subs));
rate_proxy_block_2 = reshape(nanmean(succ_block_2_valid, 2), N_SYMBOLS, sum(valid_subs));
[~, sort_ind_1] = sort(rate_proxy_block_1, 1);
[~, sort_ind_2] = sort(rate_proxy_block_2, 1);

succ_block_1_ordered = nan(size(succ_block_1_valid));
succ_block_2_ordered = nan(size(succ_block_2_valid));
for i_sub = 1:size(succ_block_1_valid,3)
    succ_block_1_ordered(:, :, i_sub) = succ_block_1_valid(sort_ind_1(:, i_sub), :, i_sub);
    succ_block_2_ordered(:, :, i_sub) = succ_block_2_valid(sort_ind_2(:, i_sub), :, i_sub);
end

symbol_list_targ = {'s', 'o', 'x', '*', '+', 'd', '.', '-', '^', '<', '>', 'p'};
if plot_bool
    figure; hold on;
    for i_targ = 1:N_SYMBOLS
        errorbar([1.5 3.5], nanmean(succ_block_1_valid(i_targ, :, :), 3),...
            sqrt(nanvar(succ_block_1_valid(i_targ, :, :),[],3)./length(subject_list)),...
            ['g', symbol_list_targ{i_targ}], 'MarkerSize', 12, 'LineWidth', 2);

        errorbar([1.5 3.5], nanmean(succ_block_2_valid(i_targ, :, :), 3),...
            sqrt(nanvar(succ_block_2_valid(i_targ, :, :),[],3)./length(subject_list)),...
            ['b', symbol_list_targ{i_targ}], 'MarkerSize', 12, 'LineWidth', 2);
    end
    axis([1 4 0 1])
end

% now for super-threshold PT
valid_subs = min_pt_all > .15;
succ_supTH_block_1_valid = succ_block_1_supTH_symb(:,:,valid_subs);
succ_supTH_block_2_valid = succ_block_2_supTH_symb(:,:,valid_subs);

rate_proxy_block_1 = reshape(nanmean(succ_supTH_block_1_valid, 2), N_SYMBOLS, sum(valid_subs));
rate_proxy_block_2 = reshape(nanmean(succ_supTH_block_2_valid, 2), N_SYMBOLS, sum(valid_subs));
[~, sort_ind_1] = sort(rate_proxy_block_1, 1);
[~, sort_ind_2] = sort(rate_proxy_block_2, 1);

succ_supTH_block_1_ordered = nan(size(succ_supTH_block_1_valid));
succ_supTH_block_2_ordered = nan(size(succ_supTH_block_2_valid));
for i_sub = 1:size(succ_supTH_block_1_valid,3)
    succ_supTH_block_1_ordered(:, :, i_sub) = succ_supTH_block_1_valid(sort_ind_1(:, i_sub), :, i_sub);
    succ_supTH_block_2_ordered(:, :, i_sub) = succ_supTH_block_2_valid(sort_ind_2(:, i_sub), :, i_sub);
end

if plot_bool
    figure; hold on;
    for i_targ = 1:N_SYMBOLS
        errorbar([1.5 3.5], nanmean(succ_supTH_block_1_ordered(i_targ, :, :), 3),...
            sqrt(nanvar(succ_supTH_block_1_ordered(i_targ, :, :),[],3)./length(subject_list)),...
            ['g', symbol_list_targ{i_targ}], 'MarkerSize', 12, 'LineWidth', 2);

        errorbar([1.5 3.5], nanmean(succ_supTH_block_2_ordered(i_targ, :, :), 3),...
            sqrt(nanvar(succ_supTH_block_2_ordered(i_targ, :, :),[],3)./length(subject_list)),...
            ['b', symbol_list_targ{i_targ}], 'MarkerSize', 12, 'LineWidth', 2);
    end
    axis([1 4 0 1])
end
% [a,b,c] = anova1(succ_block_2_valid');
% [x,y,z,q] = ttest(succ_block_2_valid(1,:) - succ_block_2_valid(end,:));

%% same but with each symbol individually & computed marginally
valid_subs = min_pt_all > .15;
% valid_subs([4, 6, 7, 21, 22]) = 1; % include subs who need modified mPT
% valid_subs(10) = 0; % exclude subs who didn't do task
succ_block_1_valid = succ_block_1_symb_marginal(:,valid_subs);
succ_block_2_valid = succ_block_2_symb_marginal(:,valid_subs);

rate_proxy_block_1 = reshape(succ_block_1_valid, N_SYMBOLS, sum(valid_subs));
rate_proxy_block_2 = reshape(succ_block_2_valid, N_SYMBOLS, sum(valid_subs));
[~, sort_ind_1] = sort(rate_proxy_block_1, 1);
[~, sort_ind_2] = sort(rate_proxy_block_2, 1);

succ_block_1_ordered = nan(size(succ_block_1_valid));
succ_block_2_ordered = nan(size(succ_block_2_valid));
for i_sub = 1:size(succ_block_1_valid,2)
    succ_block_1_ordered(:, i_sub) = succ_block_1_valid(sort_ind_1(:, i_sub), i_sub);
    succ_block_2_ordered(:, i_sub) = succ_block_2_valid(sort_ind_2(:, i_sub), i_sub);
end

symbol_list_targ = {'s', 'o', 'x', '*', '+', 'd', '.', '-', '^', '<', '>', 'p'};
if plot_bool
    figure; hold on;
    offset_ = linspace(-.2, .2, N_SYMBOLS);
    offset_set = randperm(length(offset_));
    for i_targ = 1:N_SYMBOLS
        errorbar(2.5 + offset_(offset_set(i_targ)), nanmean(succ_block_1_valid(i_targ, :), 2),...
            sqrt(nanvar(succ_block_1_valid(i_targ, :),[],2)./length(subject_list)),...
            ['g', symbol_list_targ{i_targ}], 'MarkerSize', 12, 'LineWidth', 2);

        errorbar(2.5 + offset_(offset_set(i_targ)), nanmean(succ_block_2_valid(i_targ, :), 2),...
            sqrt(nanvar(succ_block_2_valid(i_targ, :),[],2)./length(subject_list)),...
            ['b', symbol_list_targ{i_targ}], 'MarkerSize', 12, 'LineWidth', 2);
    end
    axis([1 4 0 1])
end

% now for super-threshold PT
valid_subs = min_pt_all > .15;
succ_supTH_block_1_valid = succ_block_1_supTH_symb_marginal(:,valid_subs);
succ_supTH_block_2_valid = succ_block_2_supTH_symb_marginal(:,valid_subs);

rate_proxy_block_1 = reshape(succ_supTH_block_1_valid, N_SYMBOLS, sum(valid_subs));
rate_proxy_block_2 = reshape(succ_supTH_block_2_valid, N_SYMBOLS, sum(valid_subs));
[~, sort_ind_1] = sort(rate_proxy_block_1, 1);
[~, sort_ind_2] = sort(rate_proxy_block_2, 1);

succ_supTH_block_1_ordered = nan(size(succ_supTH_block_1_valid));
succ_supTH_block_2_ordered = nan(size(succ_supTH_block_2_valid));
for i_sub = 1:size(succ_supTH_block_1_valid,2)
    succ_supTH_block_1_ordered(:, i_sub) = succ_supTH_block_1_valid(sort_ind_1(:, i_sub), i_sub);
    succ_supTH_block_2_ordered(:, i_sub) = succ_supTH_block_2_valid(sort_ind_2(:, i_sub), i_sub);
end

if plot_bool
    figure; hold on;
    for i_targ = 1:N_SYMBOLS
        errorbar(2.5 + offset_(offset_set(i_targ)), nanmean(succ_supTH_block_1_ordered(i_targ, :), 2),...
            sqrt(nanvar(succ_supTH_block_1_ordered(i_targ, :),[],2)./length(subject_list)),...
            ['g', symbol_list_targ{i_targ}], 'MarkerSize', 12, 'LineWidth', 2);

        errorbar(2.5 + offset_(offset_set(i_targ)), nanmean(succ_supTH_block_2_ordered(i_targ, :), 2),...
            sqrt(nanvar(succ_supTH_block_2_ordered(i_targ, :),[],2)./length(subject_list)),...
            ['b', symbol_list_targ{i_targ}], 'MarkerSize', 12, 'LineWidth', 2);
    end
    axis([.5 4.5 0 1])
end
% [a,b,c] = anova1(succ_block_2_valid');
% [x,y,z,q] = ttest(succ_block_2_valid(1,:) - succ_block_2_valid(end,:));
%%
if plot_bool
    figure;
    hist(min_pt_all(valid_subs))
end

data_grp.viewtime_all = viewtime_all;
data_grp.Target_all = targ_all;
data_grp.DE_all = de_all;
data_grp.Type_all = type_all;
data_grp.minPT_all = min_pt_all;
data_grp.Mov_class = move_class_all;
data_grp.succ_block_1 = succ_block_1;
data_grp.succ_block_2 = succ_block_2;
data_grp.succ_block_1_supTH = succ_block_1_supTH;
data_grp.succ_block_2_supTH = succ_block_2_supTH;
data_grp.succ_block_1_hiRes = succ_block_1_res;
data_grp.succ_block_2_hiRes = succ_block_2_res;
data_grp.succ_block_1_supTH_hiRes = succ_block_1_supTH_res;
data_grp.succ_block_2_supTH_hiRes = succ_block_2_supTH_res;
data_grp.succ_block_1_targ = succ_block_1_targ;
data_grp.succ_block_2_targ = succ_block_2_targ;
data_grp.succ_block_1_supTH_targ = succ_block_1_supTH_targ;
data_grp.succ_block_2_supTH_targ = succ_block_2_supTH_targ;
data_grp.succ_block_1_symb = succ_block_1_symb;
data_grp.succ_block_2_symb = succ_block_2_symb;
data_grp.succ_block_1_supTH_symb = succ_block_1_supTH_symb;
data_grp.succ_block_2_supTH_symb = succ_block_2_supTH_symb;
data_grp.succ_block_1_symb_marginal = succ_block_1_symb_marginal;
data_grp.succ_block_2_symb_marginal = succ_block_2_symb_marginal;
data_grp.succ_block_1_supTH_symb_marginal = succ_block_1_supTH_symb_marginal;
data_grp.succ_block_2_supTH_symb_marginal = succ_block_2_supTH_symb_marginal;
% data_grp.conf_mat = conf_mat;

err_grp = individ_analysis_errors;