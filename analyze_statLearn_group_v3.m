% Main script for analyzing statistical learning variant of VMA. Mostly a
% wrapper and book-keeping for group and individual analyses. Also allows
% across-group analyses.
%
% Version Note:
% v1 - first full/stable analysis code across groups. Does: Finds individual 
% minPT, individual direction error, succ rates in windows across trials
% for above- and below-minPT threshold demonstrating learning of
% associations, computes probability of recall of each symbol (given below
% minPT) as a function of number of appearances of that symbol,
% investigates potential differences in recallability/learnability of
% symbols and of targerts, investigates order of learning symbols to check
% if they are independent or target-specific (i.e. are the symbols for the
% same target learned at similar rates or different ones?), investigates 
% learning of symbol associations for those individuals with perfect
% post-memory scores to see if, despite perfect recall scores, there is
% still a performance deficit which might indicate an effect of motivation.
%
% v2 - same as above, plus added fourth group that did the
% 3-target/6-symbol version that included a memory test at the end of the
% first full block. 

% define globals:

global N_BINS EARLIEST_VALID_PT LATEST_VALID_PT SUCCESS_TH_ANGLE ...
    GREEN_COLOR BLUE_COLOR RED_COLOR HOME_DIR

if isempty(N_BINS)
    % if these have not been set yet (e.g. only running this analysis and
    % not main_analysis.m
    N_BINS = 6;
    SUCCESS_TH_ANGLE = 30;
    EARLIEST_VALID_PT = -.1;
    LATEST_VALID_PT = .7;

    RED_COLOR = [172, 59, 59]/255;
    GREEN_COLOR = [85, 170, 85]/255; 
    BLUE_COLOR = [86 85 149]/255;
end

LIGHT_RED = [255, 200, 200]/255;
RELATIVE_PT_MINMAX = 0.5;
TRIAL_QUARTERS = 60;

subject_list_3T = {...
'Data_SP026_03222018_E2.mat',...
'Data_S048_04052018_E2.mat',...
...% 'Data_S049_04052018_E2.mat'... %poor TR compliance
'Data_S050_04052018_E2.mat',...
'Data_S051_04052018_E2.mat',...
'Data_S052_04052018_E2.mat',...
'Data_S053_04052018_E2.mat',...
...% 'Data_S054_04052018_E2.mat',... %poor TR compliance
...% 'Data_S055_04052018_E2.mat',... %poor TR compliance
'Data_S056_04062018_E2.mat',...
'Data_S057_04062018_E2.mat',...
'Data_S058_04062018_E2.mat',...
'Data_S059_04062018_E2.mat',...
'Data_S060_04062018_E2.mat',...
'Data_S061_04062018_E2.mat',...
'Data_S062_04062018_E2.mat',...
'Data_S063_04062018_E2.mat',...
'Data_S064_04062018_E2.mat',...
'Data_S065_04062018_E2.mat',...
'Data_S066_04092018_E2.mat',...
'Data_S069_04092018_E2_fixed.mat',... 
...% 'Data_S070_04092018_E2.mat',... %poor TR compliance
'Data_S079_04112018_E2.mat',...
    };
memory_test_result_inds_3T = ...
    [1, 2, 4, 5, 6, 7, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 23];

subject_list_4T = {...
% 'Data_S089_04132018_E2.mat',... %poor TR compliance
'Data_S090_04132018_E2.mat',...
'Data_S091_04132018_E2.mat',...
'Data_S092_04132018_E2.mat',...
'Data_S093_04162018_E2.mat',...
'Data_S094_04162018_E2.mat',...
'Data_S095_04162018_E2.mat',...
'Data_S096_04162018_E2.mat',...
...% 'Data_S098_04162018_E2.mat',... %poor TR compliance
...% 'Data_S099_04162018_E2.mat',... %poor TR compliance
'Data_S100_04162018_E2.mat',...
'Data_S101_04172018_E2.mat',...
'Data_S102_04172018_E2.mat',...
'Data_S103_04172018_E2.mat',...
'Data_S104_04172018_E2.mat',...
'Data_S105_04182018_E2.mat',...
'Data_S106_04182018_E2.mat',...
'Data_S107_04182018_E2.mat',...
...% 'Data_S108_04182018_E2.mat',... %poor TR compliance
'Data_S109_04182018_E2.mat',...
    };
memory_test_result_inds_4T = ...
    [2, 3, 4, 5, 6, 7, 8, 11, 12, 13, 14, 15, 16, 17, 18, 20];

subject_list_6T = {...
'Data_S072_04102018_E2.mat',...
'Data_S074_04102018_E2.mat',...
'Data_S075_04102018_E2.mat',...
'Data_S076_04102018_E2.mat',...
'Data_S077_04112018_E2.mat',...
'Data_S078_04112018_E2.mat',...
'Data_S080_04112018_E2.mat',...
'Data_S082_04112018_E2.mat',...
'Data_S083_04122018_E2.mat',...
'Data_S084_04122018_E2.mat',...
'Data_S085_04122018_E2.mat',...
...% 'Data_S086_04122018_E2.mat',... %poor TR compliance
'Data_S087_04122018_E2.mat',...
'Data_S088_04122018_E2.mat',...
    };
memory_test_result_inds_6T = ...
    [1 2 3 4 5 6 7 8 9 10 11 13 14];

subject_list_3T_B = {...
%     'Data_S146_12042018_E2b.mat',... %poor TR compliance
    'Data_S147_12042018_E2b.mat',...    
    'Data_S148_12052018_E2b.mat',...
    'Data_S149_12052018_E2b.mat',... 
    'Data_S150_12052018_E2b.mat',... 
    'Data_S151_12062018_E2b.mat',...
    'Data_S152_12072018_E2b.mat',...
    's153.mat',... 
    's154.mat',... nan
    's255.mat',... 
    's156_fixed.mat',... %sub missed last 8 trials of intro block..filled in with Nan so size is compatible with analysis.
    's157.mat',...
    's158.mat',...
    's159.mat',...
    's160.mat',...
    's161.mat',...
    's162.mat',...
    };
memory_test_result_inds_3T_B = ...
    [2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17];

load('individual_memory_test_results.mat');

%% Compute behavioral metrics for all participants across experiments

group_subjects = {subject_list_3T, subject_list_4T, subject_list_6T,...
    subject_list_3T_B};
group_mem_test_result_inds = {memory_test_result_inds_3T, ...
    memory_test_result_inds_4T, memory_test_result_inds_6T, memory_test_result_inds_3T_B};
re_order_symbols = {[1 4 2 5 3 6], [1 5 9 2 6 10 3 7 11 4 8 12],...
    [1 7 2 8 3 9 4 10 5 11 6 12], [1 4 2 5 3 6]};
N_TARGS_GRP = [3 4 6 3];
N_SYMBS_GRP = [6 12 12 6];

min_pt_all = nan(max([length(subject_list_3T),...
    length(subject_list_4T),...
    length(subject_list_6T), length(subject_list_3T_B)]),...
    length(group_subjects));
mov_class_all = cell(1, length(group_subjects));
target_all = cell(1, length(group_subjects));
type_all = cell(1, length(group_subjects));
pt_all = cell(1, length(group_subjects));
dir_err_all = cell(1, length(group_subjects));
grp_succ_hi_1 = cell(1, length(group_subjects));
grp_succ_hi_2 = cell(1, length(group_subjects));
memory_test_results_grp = cell(1, length(group_subjects));
memory_test_results_grp4_early = nan(1, length(subject_list_3T_B));

group_analysis_list = [1, 2, 3, 4];
% group_analysis_list = 4;
for i_grp = group_analysis_list
    [data_grp, er_grp] = analyze_statLearn_group_v1(group_subjects{i_grp}, N_TARGS_GRP(i_grp), 0);
    
    if iscell(individual_memory_test_results(i_grp).score_numerator)
        % need this conditional test because group 4 (3-targ B) did two tests and so results are stored in cell. 
        memory_test_results_grp{i_grp} =...
            individual_memory_test_results(i_grp).score_numerator{2}(group_mem_test_result_inds{i_grp})./...
            individual_memory_test_results(i_grp).score_denominator{2}(group_mem_test_result_inds{i_grp});
        memory_test_results_grp4_early = ...
            individual_memory_test_results(i_grp).score_numerator{1}(group_mem_test_result_inds{i_grp})./...
            individual_memory_test_results(i_grp).score_denominator{1}(group_mem_test_result_inds{i_grp});
    else
        memory_test_results_grp{i_grp} =...
            individual_memory_test_results(i_grp).score_numerator(group_mem_test_result_inds{i_grp})./...
            individual_memory_test_results(i_grp).score_denominator(group_mem_test_result_inds{i_grp});
    end
    
    grp_succ_hi_1{i_grp} = data_grp.succ_block_1_hiRes;
    grp_succ_hi_2{i_grp} = data_grp.succ_block_2_hiRes;
    
    min_pt_all(1:length(data_grp.minPT_all), i_grp) = data_grp.minPT_all;
%     conf_mat_all{i_grp} = data_grp.conf_mat;
    mov_class_all{i_grp} = data_grp.Mov_class;
    target_all{i_grp} = data_grp.Target_all;
    type_all{i_grp} = data_grp.Type_all;
    pt_all{i_grp} = data_grp.viewtime_all;
    dir_err_all{i_grp} = data_grp.DE_all;
    
    disp(['Group ', num2str(i_grp), ' completed']);
end

%% plot statistical learning curves with memory test overlain for each group

f1 = figure;
axis_lims = [0 4.5 -0.1 1];
cond_colors = {'g', 'b'};
table_out = cell(size(group_analysis_list));
memTest_table = cell(length(group_analysis_list)+1, 1);
group_titles = {'Three Targets', 'Four Targets', 'Six Targets', 'Three Targets'};
for i_grp = group_analysis_list
    subplot(2,4,i_grp);
    hold on;
    
    errorbar(.5:.5:4, nanmean(grp_succ_hi_1{i_grp}, 2),...
        sqrt(nanvar(grp_succ_hi_1{i_grp}, [], 2)./size(grp_succ_hi_1{i_grp},2)),...
        '.-', 'Color', GREEN_COLOR, 'LineWidth', 2, 'MarkerSize', 18);
    errorbar(.5:.5:4, nanmean(grp_succ_hi_2{i_grp}, 2),...
        sqrt(nanvar(grp_succ_hi_2{i_grp}, [], 2)./size(grp_succ_hi_2{i_grp},2)),...
        '.-', 'Color', BLUE_COLOR, 'LineWidth', 2, 'MarkerSize', 18);
    
    errorbar(4.25, nanmean(memory_test_results_grp{i_grp}),...
        sqrt(nanvar(memory_test_results_grp{i_grp})./length(memory_test_results_grp{i_grp})),...
        'ks', 'LineWidth', 2);
    axis(axis_lims)
    title(group_titles{i_grp});
    xlabel('Block')
    ylabel('Probability correct');
    
    if i_grp == 4
        % show early memory test results for group 4 only (the only group
        % that did this test)
        errorbar(1.25, nanmean(memory_test_results_grp4_early),...
        sqrt(nanvar(memory_test_results_grp4_early)./length(memory_test_results_grp4_early)),...
        'ks', 'LineWidth', 2); 
    end
    
    table_temp = [...
    repmat(reshape(repmat((1:8)', 1, size(grp_succ_hi_1{i_grp},2)), numel(grp_succ_hi_1{i_grp}), 1), 2, 1),...    
    repmat(reshape(repmat((1:size(grp_succ_hi_1{i_grp}, 2)), size(grp_succ_hi_1{i_grp}, 1), 1), numel(grp_succ_hi_1{i_grp}), 1), 2, 1),...
    [ones(numel(grp_succ_hi_1{i_grp}), 1); 2*ones(numel(grp_succ_hi_2{i_grp}), 1)],...
    [reshape(grp_succ_hi_1{i_grp}, numel(grp_succ_hi_1{i_grp}), 1); reshape(grp_succ_hi_2{i_grp}, numel(grp_succ_hi_2{i_grp}), 1)],...
    ];
    
    table_out{i_grp} = table_temp;
    
    if i_grp < 4
        memTest_table{i_grp} = [...
            [zeros(length(memory_test_results_grp{i_grp}),1);...
            ones(length(memory_test_results_grp{i_grp}),1)],...
            [memory_test_results_grp{i_grp}';...
            grp_succ_hi_2{i_grp}(end, :)']...
            ];
    else
        memTest_table{i_grp} = [...
            [zeros(length(memory_test_results_grp4_early),1);...
            ones(length(memory_test_results_grp4_early),1)],...
            [memory_test_results_grp4_early';...
            grp_succ_hi_2{i_grp}(2, :)']...
            ];
        memTest_table{i_grp+1} = [...
            [zeros(length(memory_test_results_grp{i_grp}),1);...
            ones(length(memory_test_results_grp{i_grp}),1)],...
            [memory_test_results_grp{i_grp}';...
            grp_succ_hi_2{i_grp}(end, :)']...
            ];
    end
    
        
end

set(f1, 'Position', [0 500 700 500])
saveas(f1, 'statLearning_PC.pdf');

% write out variables for proper stats in R:
csvwrite('prob_succ_E2_3T', table_out{1});
csvwrite('prob_succ_E2_4T', table_out{2});
csvwrite('prob_succ_E2_6T', table_out{3});
csvwrite('prob_succ_E2_3Tb', table_out{4});
csvwrite('memTest_score_E2_3T', memTest_table{1});
csvwrite('memTest_score_E2_4T', memTest_table{2});
csvwrite('memTest_score_E2_6T', memTest_table{3});
csvwrite('memTest_score_E2_3Tb_1', memTest_table{4});
csvwrite('memTest_score_E2_3Tb_2', memTest_table{5});

% run statistics in R:
% !#bin/bash Rscript [HOME_DIR, filesep, cross_situational_learning_stats_v3.R]
% results are in 'lme_results.mat' and 'memory_test_results_*.mat'
%% compute probability of recall as function of appearances of symbol
symbol_repeat_max = [18, 8, 8, 18];

f2 = figure; 
subplot(1,5,1); hold on;
% errorbar((1:max_lpt) + 0.1 , nanmean(lpt_all_none_E1(1:max_lpt, :),2),...
%         sqrt(nanvar(lpt_all_none_E1(1:max_lpt, :), [], 2)./size(lpt_all_none_E1, 2)),...
%         'Color', LIGHT_RED, 'MarkerSize', 18, 'LineWidth', 2);
errorbar((1:max_lpt) - 0.1, nanmean(lpt_all_symbolic_E1(1:max_lpt, :),2),...
    sqrt(nanvar(lpt_all_symbolic_E1(1:max_lpt, :), [], 2)./size(lpt_all_symbolic_E1, 2)),...
    'Color', BLUE_COLOR, 'MarkerSize', 18, 'LineWidth', 2);
errorbar((1:max_lpt) + 0.1 , nanmean(lpt_all_direct_E1(1:max_lpt, :),2),...
        sqrt(nanvar(lpt_all_direct_E1(1:max_lpt, :), [], 2)./size(lpt_all_direct_E1, 2)),...
        'Color', GREEN_COLOR, 'MarkerSize', 18, 'LineWidth', 2);
axis([0 symbol_repeat_max(1) 0 1]);
xlabel('Trial occurance');
ylabel('Success probability');
for i_grp = group_analysis_list
        % for no pre-cued trials
    [rec_lpt_none, rec_hpt_none] = compute_recall_probability_appearance_order(...
        mov_class_all{i_grp}, target_all{i_grp}, type_all{i_grp},...
        pt_all{i_grp}, min_pt_all(:,i_grp), 0);
    
    lpt_per_targ_none = reshape(nanmean(rec_lpt_none, 3), size(rec_lpt_none, 1), size(rec_lpt_none, 2));
    lpt_all_none = reshape(nanmean(rec_lpt_none, 1), size(rec_lpt_none,2), size(rec_lpt_none, 3));
    
    subplot(1,5,i_grp+1); hold on;
%     errorbar(1:symbol_repeat_max(i_grp), nanmean(lpt_all_none(1:symbol_repeat_max(i_grp), :),2),...
%         sqrt(nanvar(lpt_all_none(1:symbol_repeat_max(i_grp), :), [], 2)./size(lpt_all_none, 2)), 'Color', LIGHT_RED, ...
%         'MarkerSize', 18, 'LineWidth', 2);
%     axis([0 symbol_repeat_max(1) 0 1]);
    
    
    % for direct cued trials
    [rec_lpt_direct, rec_hpt_direct] = compute_recall_probability_appearance_order(...
        mov_class_all{i_grp}, target_all{i_grp}, type_all{i_grp},...
        pt_all{i_grp}, min_pt_all(:,i_grp), 1);
    
    lpt_per_targ_direct = reshape(nanmean(rec_lpt_direct, 3), N_SYMBS_GRP(i_grp), size(rec_lpt_direct, 2));
    lpt_all_direct = reshape(nanmean(rec_lpt_direct, 1), size(rec_lpt_direct,2), size(rec_lpt_direct, 3));
    
    subplot(1,5,i_grp+1); hold on;
    errorbar(1:symbol_repeat_max(i_grp), nanmean(lpt_all_direct(1:symbol_repeat_max(i_grp), :),2),...
        sqrt(nanvar(lpt_all_direct(1:symbol_repeat_max(i_grp), :), [], 2)./size(lpt_all_direct, 2)), 'Color', GREEN_COLOR,...
        'MarkerSize', 18, 'LineWidth', 2);
    axis([0 symbol_repeat_max(1) 0 1]);
    

    
    % for symbolicly cued trials
    [rec_lpt, rec_hpt] = compute_recall_probability_appearance_order(...
        mov_class_all{i_grp}, target_all{i_grp}, type_all{i_grp},...
        pt_all{i_grp}, min_pt_all(:,i_grp));
    
    lpt_per_targ = reshape(nanmean(rec_lpt, 3), N_SYMBS_GRP(i_grp), size(rec_lpt, 2));
    lpt_all = reshape(nanmean(rec_lpt, 1), size(rec_lpt,2), size(rec_lpt, 3));
    
    subplot(1,5,i_grp+1); hold on;
    errorbar(1:symbol_repeat_max(i_grp), nanmean(lpt_all(1:symbol_repeat_max(i_grp), :),2),...
        sqrt(nanvar(lpt_all(1:symbol_repeat_max(i_grp), :), [], 2)./size(lpt_all, 2)), 'Color', BLUE_COLOR, ...
        'MarkerSize', 18, 'LineWidth', 2);
    axis([0 symbol_repeat_max(1) 0 1]);
    
    xlabel('Trial occurance');
    ylabel('Success probability');
end
set(f2, 'Position', [1 197 900 300]);
set(f2, 'PaperOrientation', 'landscape')
saveas(f2, 'Symbole_learning_by_occurance.pdf');

%% compute pr(succ | pt) for first half and second half of blocks for each group:
% Doesn't go into the paper but is a nice confirmation of data.

inds_init = 1:(12+60);
inds_half_1 = inds_init(end) + (1:120);
inds_half_2 = inds_half_1(end) + (1:120);

figure;
sub_plot_inds = [1 2; 3 4; 5 6; 7 8];
hist_bins = linspace(EARLIEST_VALID_PT, LATEST_VALID_PT, N_BINS+1);
x_inds = hist_bins(1:(end-1)) + diff(hist_bins)/2;
for i_grp = group_analysis_list
    grp_pr_succ_0 = nan(N_BINS, size(pt_all{i_grp}, 2));
    grp_pr_succ_half1 = nan(N_BINS, size(pt_all{i_grp}, 2), 3);
    grp_pr_succ_half2 = nan(N_BINS, size(pt_all{i_grp}, 2), 3);
    for i_sub = 1:size(pt_all{i_grp},2)
        [grp_pr_succ_0(:, i_sub), ~, ~] = ...
         compute_succ_prob_across_pt(...
         pt_all{i_grp}(inds_init, i_sub),...
         dir_err_all{i_grp}(inds_init, i_sub),...
         type_all{i_grp}(inds_init, i_sub), N_BINS);
        
        [grp_pr_succ_half1(:, i_sub, 1),...
            grp_pr_succ_half1(:, i_sub, 2),... 
            grp_pr_succ_half1(:, i_sub, 3)] = ...
         compute_succ_prob_across_pt(...
            pt_all{i_grp}(inds_half_1, i_sub), ...
            dir_err_all{i_grp}(inds_half_1, i_sub), ...
            type_all{i_grp}(inds_half_1, i_sub), N_BINS);
        
        [grp_pr_succ_half2(:, i_sub, 1),...
            grp_pr_succ_half2(:, i_sub, 2),...
            grp_pr_succ_half2(:, i_sub, 3)] = ...
         compute_succ_prob_across_pt(...
            pt_all{i_grp}(inds_half_2, i_sub), ...
            dir_err_all{i_grp}(inds_half_2, i_sub), ...
            type_all{i_grp}(inds_half_2, i_sub), N_BINS);
    end
    subplot(4,2,sub_plot_inds(i_grp, 1)); hold on
    errorbar(x_inds, nanmean(grp_pr_succ_0,2),...
        sqrt(nanvar(grp_pr_succ_0,[],2)./sum(~isnan(grp_pr_succ_0),2)),...
        'rs');
    errorbar(x_inds, nanmean(grp_pr_succ_half1(:,:,2),2),...
        sqrt(nanvar(grp_pr_succ_half1(:,:,2),[],2)./sum(~isnan(grp_pr_succ_half1(:,:,2)),2)),...
        'gs');
    errorbar(x_inds, nanmean(grp_pr_succ_half1(:,:,3),2),...
        sqrt(nanvar(grp_pr_succ_half1(:,:,3),[],2)./sum(~isnan(grp_pr_succ_half1(:,:,3)),2)),...
        'bs');
    axis([min(hist_bins), max(hist_bins), 0, 1]);
    
    subplot(4,2,sub_plot_inds(i_grp, 2)); hold on
    errorbar(x_inds, nanmean(grp_pr_succ_0,2),...
        sqrt(nanvar(grp_pr_succ_0,[],2)./sum(~isnan(grp_pr_succ_0),2)),...
        'rs');
    errorbar(x_inds, nanmean(grp_pr_succ_half2(:,:,2),2),...
        sqrt(nanvar(grp_pr_succ_half2(:,:,2),[],2)./sum(~isnan(grp_pr_succ_half2(:,:,2)),2)),...
        'gs');
    errorbar(x_inds, nanmean(grp_pr_succ_half2(:,:,3),2),...
        sqrt(nanvar(grp_pr_succ_half2(:,:,3),[],2)./sum(~isnan(grp_pr_succ_half2(:,:,3)),2)),...
        'bs');
    axis([min(hist_bins), max(hist_bins), 0, 1]);
end

%% Align data to individual's minimum PT and combine across individuals. 
% Look at Pr(succ), response to catch trials, and variability; in quarters
% for each group

inds_init = 1:(12+60);
inds_q1 = inds_init(end) + (1:60);
inds_q2 = inds_q1(end) + (1:60);
inds_q3 = inds_q2(end) + (1:60);
inds_q4 = inds_q3(end) + (1:60);

% Plot directional error vs. PT, aligned to minPT
h1 = figure; h2 = figure;
rm_models = cell(1,4);
sub_plot_inds = [1 2 3 4 5; 6 7 8 9 10; 11 12 13 14 15; 16 17 18 19 20];
sub_plot_inds2 = [1 2 3 4; 5 6 7 8];
for i_grp = group_analysis_list
    n_subs_group = size(pt_all{i_grp}, 2);
    
    pt_q0 = pt_all{i_grp}(inds_init, :);
    pt_q1 = pt_all{i_grp}(inds_q1, :);
    pt_q2 = pt_all{i_grp}(inds_q2, :);
    pt_q3 = pt_all{i_grp}(inds_q3, :);
    pt_q4 = pt_all{i_grp}(inds_q4, :);
    type_q0 = type_all{i_grp}(inds_init, :);
    type_q1 = type_all{i_grp}(inds_q1, :);
    type_q2 = type_all{i_grp}(inds_q2, :);
    type_q3 = type_all{i_grp}(inds_q3, :);
    type_q4 = type_all{i_grp}(inds_q4, :);
    de_q0 = dir_err_all{i_grp}(inds_init, :);
    de_q1 = dir_err_all{i_grp}(inds_q1, :);
    de_q2 = dir_err_all{i_grp}(inds_q2, :);
    de_q3 = dir_err_all{i_grp}(inds_q3, :);
    de_q4 = dir_err_all{i_grp}(inds_q4, :);
    
    figure(h1)
    subplot(4,5,sub_plot_inds(i_grp, 1)); hold on;
    for i_sub = 1:size(pt_all{i_grp},2)
        plot(pt_all{i_grp}(inds_init, i_sub) - min_pt_all(i_sub, i_grp),...
            dir_err_all{i_grp}(inds_init, i_sub), 'r.');
    end
    plot([0 0], [-200 200], 'k-')
    axis([-.5 .5 -200 200])
    
    subplot(4,5,sub_plot_inds(i_grp, 2)); hold on;
    for i_sub = 1:size(pt_all{i_grp},2)
        plot(pt_q1(type_q1(:, i_sub) == 1, i_sub) - min_pt_all(i_sub, i_grp),...
            de_q1(type_q1(:, i_sub) == 1, i_sub),...
            'g.');
        plot(pt_q1(type_q1(:, i_sub) == 2, i_sub) - min_pt_all(i_sub, i_grp),...
            de_q1(type_q1(:, i_sub) == 2, i_sub),...
            'b.');
    end
    plot([0 0], [-200 200], 'k-')
    axis([-.5 .5 -200 200])
    
    subplot(4,5,sub_plot_inds(i_grp, 3)); hold on;
    for i_sub = 1:size(pt_all{i_grp},2)
        plot(pt_q2(type_q2(:, i_sub) == 1, i_sub) - min_pt_all(i_sub, i_grp),...
            de_q2(type_q2(:, i_sub) == 1, i_sub),...
            'g.');
        plot(pt_q2(type_q2(:, i_sub) == 2, i_sub) - min_pt_all(i_sub, i_grp),...
            de_q2(type_q2(:, i_sub) == 2, i_sub),...
            'b.');
    end
    plot([0 0], [-200 200], 'k-')
    axis([-.5 .5 -200 200])
    
    subplot(4,5,sub_plot_inds(i_grp, 4)); hold on;
    for i_sub = 1:size(pt_all{i_grp},2)
        plot(pt_q3(type_q3(:, i_sub) == 1, i_sub) - min_pt_all(i_sub, i_grp),...
            de_q3(type_q3(:, i_sub) == 1, i_sub),...
            'g.');
        plot(pt_q3(type_q3(:, i_sub) == 2, i_sub) - min_pt_all(i_sub, i_grp),...
            de_q3(type_q3(:, i_sub) == 2, i_sub),...
            'b.');
    end
    plot([0 0], [-200 200], 'k-')
    axis([-.5 .5 -200 200])
    
    subplot(4,5,sub_plot_inds(i_grp, 5)); hold on;
    for i_sub = 1:size(pt_all{i_grp},2)
        plot(pt_q4(type_q4(:, i_sub) == 1, i_sub) - min_pt_all(i_sub, i_grp),...
            de_q4(type_q4(:, i_sub) == 1, i_sub),...
            'g.');
        plot(pt_q4(type_q4(:, i_sub) == 2, i_sub) - min_pt_all(i_sub, i_grp),...
            de_q4(type_q4(:, i_sub) == 2, i_sub),...
            'b.');
    end
    plot([0 0], [-200 200], 'k-')
    axis([-.5 .5 -200 200])
    
    figure(h2);
    % compute variability by quartile for each trial type:
    pt_q = {pt_q1, pt_q2, pt_q3, pt_q4};
    de_q = {de_q1, de_q2, de_q3, de_q4};
    mpt = min_pt_all(:, i_grp);

    reach_var_1_set = nan(size(min_pt_all,1), 4);
    reach_var_2_set = nan(size(min_pt_all,1), 4);
    for i_q = 1:4
        this_pt = nan(30, size(min_pt_all,1), 3); %trials (more than necessary) x subjs x type
        this_de = nan(30, size(min_pt_all,1), 3);
        for i_type = 1:2
            for i_sub = 1:size(type_q1,2)
                pt_temp = pt_q{i_q}(type_q1(:, i_sub) == i_type, i_sub);
                de_temp = de_q{i_q}(type_q1(:, i_sub) == i_type, i_sub);

                this_pt(1:size(pt_temp,1), i_sub, i_type+1) = pt_temp;
                this_de(1:size(de_temp,1), i_sub, i_type+1) = de_temp;
            end
        end
        [reach_var_0, reach_var_1, reach_var_2] = ...
            compute_variability_by_pt(this_pt, this_de, mpt);
        
        reach_var_1_set(:, i_q) = reach_var_1(:,1);
        reach_var_2_set(:, i_q) = reach_var_2(:,1);

    end
    subplot(2,4,i_grp); hold on;
    errorbar(1:4, nanmean(reach_var_1_set, 1),...
        sqrt(nanvar(reach_var_1_set, [], 1)./sum(~isnan(reach_var_1_set),1)), 'g.-');
    errorbar(1:4, nanmean(reach_var_2_set, 1),...
        sqrt(nanvar(reach_var_2_set, [], 1)./sum(~isnan(reach_var_2_set),1)), 'b.-');
    axis([0.5 4.5 0 10])
    
    subplot(2,4,sub_plot_inds2(2,i_grp)); hold on;
    temp_var = reach_var_2_set - reach_var_1_set;
    errorbar(1:4, nanmean(temp_var, 1),...
        sqrt(nanvar(temp_var, [], 1)./sum(~isnan(temp_var),1)), 'k.-');
    axis([0.5 4.5 -5 5])
    
%     reach_var_2_valid = reach_var_2_set(sum(isnan(reach_var_2_set),2) < 1 & sum(isnan(reach_var_1_set),2) < 1, :);
%     reach_var_1_valid = reach_var_1_set(sum(isnan(reach_var_2_set),2) < 1 & sum(isnan(reach_var_1_set),2) < 1, :);
    reach_var_2_valid = reach_var_2_set(1:n_subs_group, :);
    reach_var_1_valid = reach_var_1_set(1:n_subs_group, :);
    quantile_valid = repmat(1:4, size(reach_var_2_valid,1), 1);
    subjects_valid = repmat((1:size(reach_var_1_valid,1))', 4, 1);
    var_table_2 = reshape(reach_var_2_valid, numel(reach_var_2_valid), 1);
    var_table_1 = reshape(reach_var_1_valid, numel(reach_var_1_valid), 1);
    quantile_table = reshape(quantile_valid, numel(quantile_valid), 1);
    t_table = table(quantile_table, var_table_1, var_table_2, ...
        'VariableNames', {'quantiles','var_1', 'var_2'});
    rm = fitrm(t_table, 'var_1-var_2~quantiles');
    rm_models{i_grp} = rm;
    
    % write out table to do analysis properly in R:
    csvwrite(['variability_', num2str(i_grp)], ...
        [[quantile_table; quantile_table], ...
        [subjects_valid; subjects_valid],...
        [ones(length(quantile_table), 1); 2*ones(length(quantile_table),1)],...
        [var_table_1; var_table_2]]);
end
set(h1, 'Position', [72 549 818 406]);
set(h1, 'PaperOrientation', 'landscape')
saveas(h1, 'Groups_by_quarter_PC_scatter.pdf');
%% 

% Plot success vs. PT and var. vs. PT, aligned to minPT
h1 = figure;
h2 = figure;
% sub_plot_inds = [1 2 3 4; 5 6 7 8; 9 10 11 12; 13 14 15 16];
sub_plot_inds = [1 2 3 4 5; 6 7 8 9 10; 11 12 13 14 15; 16 17 18 19 20];
% override globals:
n_bins = 12;
hist_bins = linspace(-RELATIVE_PT_MINMAX, RELATIVE_PT_MINMAX, n_bins+1);
x_inds = hist_bins(1:(end-1)) + diff(hist_bins)/2;
for i_grp = group_analysis_list
    grp_min_pt = min_pt_all(~isnan(min_pt_all(:, i_grp)), i_grp);
    pt_init = pt_all{i_grp}(inds_init, :) - repmat(grp_min_pt', length(inds_init), 1);
    pt_q1 = pt_all{i_grp}(inds_q1, :)  - repmat(grp_min_pt', length(inds_q1), 1);
    pt_q2 = pt_all{i_grp}(inds_q2, :) - repmat(grp_min_pt', length(inds_q2), 1);
    pt_q3 = pt_all{i_grp}(inds_q3, :)  - repmat(grp_min_pt', length(inds_q3), 1);
    pt_q4 = pt_all{i_grp}(inds_q4, :) - repmat(grp_min_pt', length(inds_q4), 1);
    type_init = type_all{i_grp}(inds_init, :);
    type_q1 = type_all{i_grp}(inds_q1, :);
    type_q2 = type_all{i_grp}(inds_q2, :);
    type_q3 = type_all{i_grp}(inds_q3, :);
    type_q4 = type_all{i_grp}(inds_q4, :);
    de_init = dir_err_all{i_grp}(inds_init, :);
    de_q1 = dir_err_all{i_grp}(inds_q1, :);
    de_q2 = dir_err_all{i_grp}(inds_q2, :);
    de_q3 = dir_err_all{i_grp}(inds_q3, :);
    de_q4 = dir_err_all{i_grp}(inds_q4, :);
    
    [ps0, ~, ~] = compute_succ_prob_across_pt(...
        pt_init(:), de_init(:), type_init(:), n_bins);
    
    var_inds_init = abs(de_init(:)) < SUCCESS_TH_ANGLE;
    [vr0, ~, ~] = compute_var_across_pt(...
        pt_init(var_inds_init),...
        de_init(var_inds_init),...
        type_init(var_inds_init),...
        n_bins);
    % plot initial block with only type 0 trials:
    figure(h1)
    subplot(4,5,sub_plot_inds(i_grp, 1)); hold on;
    plot(x_inds, ps0, 'r.-');
    plot([0 0], [0 1], 'k-')
    axis([-.5 .5 0 1])
    plot([0 0], [0 1], 'k-')
    
    figure(h2)
    subplot(4,5,sub_plot_inds(i_grp, 1)); hold on;
    plot(x_inds, vr0, 'r.-');
    axis([-.5 .5 0 20])
    plot([0 0], [0 20], 'k-')
    
    % compile data matrix to include each quarter:
    data_mat_pc_q0 = [abs(de_init(:)) < 30,... % directional error
    type_init(:),... % type
    pt_init(:),... % preparation time
    reshape(repmat((1:size(pt_init,2)), size(pt_init,1), 1), numel(pt_init), 1),... % subject number
    reshape(repmat(zeros(size(pt_init,1),1), 1, size(pt_init,2)), numel(pt_init), 1) ... % trial number
    ];
    
    data_mat_pc_q1 = [abs(de_q1(:)) < 30,... % directional error
    type_q1(:),... % type
    pt_q1(:),... % preparation time
    reshape(repmat((1:size(pt_q1,2)), size(pt_q1,1), 1), numel(pt_q1), 1),... % subject number
    reshape(repmat((1:size(pt_q1,1))', 1, size(pt_q1,2)), numel(pt_q1), 1) ... % trial number
    ];
    
    data_mat_pc_q2 = [abs(de_q2(:)) < 30,...
    type_q2(:),...
    pt_q2(:),...
    reshape(repmat((1:size(pt_q2,2)), size(pt_q2,1), 1), numel(pt_q2), 1),...
    reshape(repmat(TRIAL_QUARTERS + (1:size(pt_q2,1))', 1, size(pt_q2,2)), numel(pt_q2), 1)... % trial number
    ];

    data_mat_pc_q3 = [abs(de_q3(:)) < 30,...
    type_q3(:),...
    pt_q3(:),...
    reshape(repmat((1:size(pt_q3,2)), size(pt_q3,1), 1), numel(pt_q3), 1),...
    reshape(repmat(2*TRIAL_QUARTERS + (1:size(pt_q3,1))', 1, size(pt_q3,2)), numel(pt_q3), 1)... % trial number
    ];
    
    data_mat_pc_q4 = [abs(de_q4(:)) < 30,...
    type_q4(:),...
    pt_q4(:),...
    reshape(repmat((1:size(pt_q4,2)), size(pt_q4,1), 1), numel(pt_q4), 1),...
    reshape(repmat(3*TRIAL_QUARTERS + (1:size(pt_q4,1))', 1, size(pt_q4,2)), numel(pt_q4), 1)... % trial number
    ];

    % To use quarters of all trials: 
    data_mat_pc_quarters = [...
        data_mat_pc_q0(:,1:4), zeros(size(data_mat_pc_q0,1), 1);...
        data_mat_pc_q1(:,1:4), ones(size(data_mat_pc_q1,1), 1);...
        data_mat_pc_q2(:,1:4), 2*ones(size(data_mat_pc_q2,1), 1);...
        data_mat_pc_q3(:,1:4), 3*ones(size(data_mat_pc_q3,1), 1);...
        data_mat_pc_q4(:,1:4), 4*ones(size(data_mat_pc_q4,1), 1)];
    
    % To use numeric trial: 
    data_mat_pc_trials = [...
        data_mat_pc_q0(:,1:5);...
        data_mat_pc_q1(:,1:5);...
        data_mat_pc_q2(:,1:5);...
        data_mat_pc_q3(:,1:5);...
        data_mat_pc_q4(:,1:5)];

    csvwrite(['raw_data_mat_E2_pc_quarters', num2str(i_grp)], data_mat_pc_quarters);
    csvwrite(['raw_data_mat_E2_pc_trials', num2str(i_grp)], data_mat_pc_trials);
    
    [~, ps1, ps2] = compute_succ_prob_across_pt(...
        pt_q1(:), de_q1(:), type_q1(:), n_bins);
    
    var_inds_half1 = abs(de_q1(:)) < SUCCESS_TH_ANGLE;
    [~, vr1, vr2] = compute_var_across_pt(...
        pt_q1(var_inds_half1),...
        de_q1(var_inds_half1),...
        type_q1(var_inds_half1),...
        n_bins);    
    
    % plot block 1: type 1 & 2 trials
    figure(h1)
    subplot(4,5,sub_plot_inds(i_grp, 2)); hold on;
%     plot(x_inds, ps0, 'r.-');
    plot(x_inds, ps1, 'g.-');
    plot(x_inds, ps2, 'b.-');
    plot([0 0], [0 1], 'k-')
    axis([-.5 .5 0 1])
    plot([0 0], [0 1], 'k-')
    
    figure(h2)
    subplot(4,5,sub_plot_inds(i_grp, 2)); hold on;
%     plot(x_inds, vr0, 'r.-');
    plot(x_inds, vr1, 'g.-');
    plot(x_inds, vr2, 'b.-');
    axis([-.5 .5 0 20])
    plot([0 0], [0 20], 'k-')
    
    [~, ps1, ps2] = compute_succ_prob_across_pt(...
        pt_q2(:), de_q2(:), type_q2(:), n_bins);
    
    var_inds_half2 = abs(de_q2(:)) < SUCCESS_TH_ANGLE;
    [~, vr1, vr2] = compute_var_across_pt(...
        pt_q2(var_inds_half2),...
        de_q2(var_inds_half2),...
        type_q2(var_inds_half2),...
        n_bins);
    % plot block 2: type 1 & 2 trials
    figure(h1)
    subplot(4,5,sub_plot_inds(i_grp, 3)); hold on;
%     plot(x_inds, ps0, 'r.-');
    plot(x_inds, ps1, 'g.-');
    plot(x_inds, ps2, 'b.-');
    plot([0 0], [0 1], 'k-')
    axis([-.5 .5 0 1])
    plot([0 0], [0 1], 'k-')
    
    figure(h2)
    subplot(4,5,sub_plot_inds(i_grp, 3)); hold on;
%     plot(x_inds, vr0, 'r.-');
    plot(x_inds, vr1, 'g.-');
    plot(x_inds, vr2, 'b.-');
    axis([-.5 .5 0 20])
    plot([0 0], [0 20], 'k-')
    
    [~, ps1, ps2] = compute_succ_prob_across_pt(...
        pt_q3(:), de_q3(:), type_q3(:), n_bins);
    
    var_inds_half2 = abs(de_q3(:)) < SUCCESS_TH_ANGLE;
    [~, vr1, vr2] = compute_var_across_pt(...
        pt_q3(var_inds_half2),...
        de_q3(var_inds_half2),...
        type_q3(var_inds_half2),...
        n_bins);
    % plot block 3: type 1 & 2 trials
    figure(h1)
    subplot(4,5,sub_plot_inds(i_grp, 4)); hold on;
%     plot(x_inds, ps0, 'r.-');
    plot(x_inds, ps1, 'g.-');
    plot(x_inds, ps2, 'b.-');
    plot([0 0], [0 1], 'k-')
    axis([-.5 .5 0 1])
    plot([0 0], [0 1], 'k-')
    
    figure(h2)
    subplot(4,5,sub_plot_inds(i_grp, 4)); hold on;
%     plot(x_inds, vr0, 'r.-');
    plot(x_inds, vr1, 'g.-');
    plot(x_inds, vr2, 'b.-');
    axis([-.5 .5 0 20])
    plot([0 0], [0 20], 'k-')
    
    [~, ps1, ps2] = compute_succ_prob_across_pt(...
        pt_q4(:), de_q4(:), type_q4(:), n_bins);
    
    var_inds_half2 = abs(de_q4(:)) < SUCCESS_TH_ANGLE;
    [~, vr1, vr2] = compute_var_across_pt(...
        pt_q4(var_inds_half2),...
        de_q4(var_inds_half2),...
        type_q4(var_inds_half2),...
        n_bins);
    % plot block 4: type 1 & 2 trials
    figure(h1)
    subplot(4,5,sub_plot_inds(i_grp, 5)); hold on;
%     plot(x_inds, ps0, 'r.-');
    plot(x_inds, ps1, 'g.-');
    plot(x_inds, ps2, 'b.-');
    plot([0 0], [0 1], 'k-')
    axis([-.5 .5 0 1])
    plot([0 0], [0 1], 'k-')
    
    figure(h2)
    subplot(4,5,sub_plot_inds(i_grp, 5)); hold on;
%     plot(x_inds, vr0, 'r.-');
    plot(x_inds, vr1, 'g.-');
    plot(x_inds, vr2, 'b.-');
    axis([-.5 .5 0 20])
    plot([0 0], [0 20], 'k-')
end


%%
% Plot directional error vs. PT of catch trials , aligned to minPT
fig_catch = figure;
sub_plot_inds = [1 2 3 4 5; 6 7 8 9 10; 11 12 13 14 15; 16 17 18 19 20];
for i_grp = group_analysis_list
    pt_q1 = pt_all{i_grp}(inds_q1, :);
    pt_q2 = pt_all{i_grp}(inds_q2, :);
    pt_q3 = pt_all{i_grp}(inds_q3, :);
    pt_q4 = pt_all{i_grp}(inds_q4, :);
    type_q1 = type_all{i_grp}(inds_q1, :);
    type_q2 = type_all{i_grp}(inds_q2, :);
    type_q3 = type_all{i_grp}(inds_q3, :);
    type_q4 = type_all{i_grp}(inds_q4, :);
    de_q1 = dir_err_all{i_grp}(inds_q1, :);
    de_q2 = dir_err_all{i_grp}(inds_q2, :);
    de_q3 = dir_err_all{i_grp}(inds_q3, :);
    de_q4 = dir_err_all{i_grp}(inds_q4, :);   
    
    subplot(4,5,sub_plot_inds(i_grp, 2)); hold on;
    for i_sub = 1:size(pt_all{i_grp},2)
        plot(pt_q1(type_q1(:, i_sub) == 3, i_sub) - min_pt_all(i_sub, i_grp),...
            de_q1(type_q1(:, i_sub) == 3, i_sub),...
            'go');
        plot(pt_q1(type_q1(:, i_sub) == 4, i_sub) - min_pt_all(i_sub, i_grp),...
            de_q1(type_q1(:, i_sub) == 4, i_sub),...
            'bo');
    end
    plot([0 0], [-200 200], 'k-')
    axis([-.5 .5 -200 200])
    
    subplot(4,5,sub_plot_inds(i_grp, 3)); hold on;
    for i_sub = 1:size(pt_all{i_grp},2)
        plot(pt_q2(type_q2(:, i_sub) == 3, i_sub) - min_pt_all(i_sub, i_grp),...
            de_q2(type_q2(:, i_sub) == 3, i_sub),...
            'go');
        plot(pt_q2(type_q2(:, i_sub) == 4, i_sub) - min_pt_all(i_sub, i_grp),...
            de_q2(type_q2(:, i_sub) == 4, i_sub),...
            'bo');
    end
    plot([0 0], [-200 200], 'k-')
    axis([-.5 .5 -200 200])
    
    subplot(4,5,sub_plot_inds(i_grp, 4)); hold on;
    for i_sub = 1:size(pt_all{i_grp},2)
        plot(pt_q3(type_q3(:, i_sub) == 3, i_sub) - min_pt_all(i_sub, i_grp),...
            de_q3(type_q3(:, i_sub) == 3, i_sub),...
            'go');
        plot(pt_q3(type_q3(:, i_sub) == 4, i_sub) - min_pt_all(i_sub, i_grp),...
            de_q3(type_q3(:, i_sub) == 4, i_sub),...
            'bo');
    end
    plot([0 0], [-200 200], 'k-')
    axis([-.5 .5 -200 200])
    
    subplot(4,5,sub_plot_inds(i_grp, 5)); hold on;
    for i_sub = 1:size(pt_all{i_grp},2)
        plot(pt_q4(type_q4(:, i_sub) == 3, i_sub) - min_pt_all(i_sub, i_grp),...
            de_q4(type_q4(:, i_sub) == 3, i_sub),...
            'go');
        plot(pt_q4(type_q4(:, i_sub) == 4, i_sub) - min_pt_all(i_sub, i_grp),...
            de_q4(type_q4(:, i_sub) == 4, i_sub),...
            'bo');
    end
    plot([0 0], [-200 200], 'k-')
    axis([-.5 .5 -200 200])
end
set(fig_catch, 'Position', [72 549 818 406]);
set(fig_catch, 'PaperOrientation', 'landscape')
saveas(fig_catch, 'Groups_PC_catch_scatter.pdf');