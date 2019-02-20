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
global n_bins EARLIEST_VALID_PT LATEST_VALID_PT SUCCESS_TH_ANGLE
n_bins = 6;
SUCCESS_TH_ANGLE = 30;

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
% 'Data_S089_04132018_E2.mat',...
'Data_S090_04132018_E2.mat',...
'Data_S091_04132018_E2.mat',...
'Data_S092_04132018_E2.mat',...
'Data_S093_04162018_E2.mat',...
'Data_S094_04162018_E2.mat',...
'Data_S095_04162018_E2.mat',...
'Data_S096_04162018_E2.mat',...
...% 'Data_S098_04162018_E2.mat',...
...% 'Data_S099_04162018_E2.mat',...
'Data_S100_04162018_E2.mat',...
'Data_S101_04172018_E2.mat',...
'Data_S102_04172018_E2.mat',...
'Data_S103_04172018_E2.mat',...
'Data_S104_04172018_E2.mat',...
'Data_S105_04182018_E2.mat',...
'Data_S106_04182018_E2.mat',...
'Data_S107_04182018_E2.mat',...
...% 'Data_S108_04182018_E2.mat',...
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
...% 'Data_S086_04122018_E2.mat',...
'Data_S087_04122018_E2.mat',...
'Data_S088_04122018_E2.mat',...
    };
memory_test_result_inds_6T = ...
    [1 2 3 4 5 6 7 8 9 10 11 13 14];

subject_list_3T_B = {...
%     'Data_S146_12042018_E2b.mat',...    
    'Data_S147_12042018_E2b.mat',...    
    'Data_S148_12052018_E2b.mat',...
    'Data_S149_12052018_E2b.mat',... 
    'Data_S150_12052018_E2b.mat',... 
    'Data_S151_12062018_E2b.mat',...
    'Data_S152_12072018_E2b.mat',...
    's153.mat',...
    's255.mat',...
    's156.mat',...
    's157.mat',...
    's158.mat',...
    's159.mat',...
    's160.mat',...
    's161.mat',...
    };
memory_test_result_inds_3T_B = ...
    [2 3 4 5 6 7];

load('individual_memory_test_results.mat');

%%

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

group_analysis_list = [1, 2, 3, 4];
% group_analysis_list = 4;
for i_grp = group_analysis_list
    [data_grp, er_grp] = analyze_statLearn_group_v1(group_subjects{i_grp}, N_TARGS_GRP(i_grp), 0);
    
    if iscell(individual_memory_test_results(i_grp).score_numerator)
        % need this conditional test because group 4 (3-targ B) did two tests and so results are stored in cell. 
        memory_test_results_grp{i_grp} =...
            individual_memory_test_results(i_grp).score_numerator{2}(group_mem_test_result_inds{i_grp})./...
            individual_memory_test_results(i_grp).score_denominator{2}(group_mem_test_result_inds{i_grp});
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

    % compare symbol learn order
    strat = data_grp.succ_block_2_symb_marginal(re_order_symbols{i_grp}, :);
    [~, strat_sort] = sort(strat, 'descend');
    figure;
    k_order = 1;
    order_list = 1:2:N_SYMBS_GRP(i_grp);
    for i_order = order_list
        subplot(1,N_SYMBS_GRP(i_grp)/2,k_order)
        bar(1:N_SYMBS_GRP(i_grp), sum(strat_sort == i_order | ...
            strat_sort == (i_order + 1), 2));
        axis([.5 N_SYMBS_GRP(i_grp) + .5 0 10])
        k_order = k_order + 1;
    end
    
    disp(['Group ', num2str(i_grp), ' completed']);
end

%% plot statistical learning curves with memory test overlain for each group
figure;
axis_lims = [0 4.5 -0.1 1];
cond_colors = {'g', 'b'};
for i_grp = group_analysis_list
    subplot(2,4,i_grp);
    hold on;
    
    errorbar(.5:.5:4, nanmean(grp_succ_hi_1{i_grp}, 2),...
        sqrt(nanvar(grp_succ_hi_1{i_grp}, [], 2)./size(grp_succ_hi_1{i_grp},2)),...
        [cond_colors{1}, '.-']);
    errorbar(.5:.5:4, nanmean(grp_succ_hi_2{i_grp}, 2),...
        sqrt(nanvar(grp_succ_hi_2{i_grp}, [], 2)./size(grp_succ_hi_2{i_grp},2)),...
        [cond_colors{2}, '.-']);
    
    errorbar(4.5, nanmean(memory_test_results_grp{i_grp}),...
        sqrt(nanvar(memory_test_results_grp{i_grp})./length(memory_test_results_grp{i_grp})),...
        'ks');
    axis(axis_lims)
    title(['T=', num2str(N_TARGS_GRP(i_grp)), ', S=', num2str(N_SYMBS_GRP(i_grp))]);
end


% plot those from group 1 with perfect memory scores
if ismember(1, group_analysis_list)
    grp_succ_hi_2_perfect = grp_succ_hi_2{1}(:, memory_test_results_grp{1} == 1);
    grp_succ_hi_1_perfect = grp_succ_hi_1{1}(:, memory_test_results_grp{1} == 1);
    subplot(2,4,5); hold on;
    errorbar(.5:.5:4, nanmean(grp_succ_hi_1_perfect, 2), ...
        sqrt(nanvar(grp_succ_hi_1_perfect, [], 2)./size(grp_succ_hi_1_perfect,2)),...
        [cond_colors{1}, '.-']);
    errorbar(.5:.5:4, nanmean(grp_succ_hi_2_perfect, 2), ...
        sqrt(nanvar(grp_succ_hi_2_perfect, [], 2)./size(grp_succ_hi_2_perfect,2)),...
        [cond_colors{2}, '.-']);
    errorbar(.5:.5:4, nanmean(grp_succ_hi_1_perfect - grp_succ_hi_2_perfect, 2),...
        sqrt(nanvar(grp_succ_hi_1_perfect - grp_succ_hi_2_perfect, [], 2)./...
        size(grp_succ_hi_2_perfect,2)),...
        'k.-');
    plot(axis_lims(1:2), zeros(1,2), 'k-')
    axis(axis_lims)
    title(['perf. ', 'T=', num2str(N_TARGS_GRP(i_grp)), ', S=', num2str(N_SYMBS_GRP(i_grp))]);
end

% plot those from group 2 with perfect memory scores
if ismember(2, group_analysis_list)
    grp_succ_hi_2_perfect = grp_succ_hi_2{2}(:, memory_test_results_grp{2} == 1);
    grp_succ_hi_1_perfect = grp_succ_hi_1{2}(:, memory_test_results_grp{2} == 1);
    subplot(2,4,6); hold on;
    errorbar(.5:.5:4, nanmean(grp_succ_hi_1_perfect, 2), ...
        sqrt(nanvar(grp_succ_hi_1_perfect, [], 2)./size(grp_succ_hi_1_perfect,2)),...
        [cond_colors{1}, '.-']);
    errorbar(.5:.5:4, nanmean(grp_succ_hi_2_perfect, 2), ...
        sqrt(nanvar(grp_succ_hi_2_perfect, [], 2)./size(grp_succ_hi_2_perfect,2)),...
        [cond_colors{2}, '.-']);
    errorbar(.5:.5:4, nanmean(grp_succ_hi_1_perfect - grp_succ_hi_2_perfect, 2),...
        sqrt(nanvar(grp_succ_hi_1_perfect - grp_succ_hi_2_perfect, [], 2)./...
        size(grp_succ_hi_2_perfect,2)),...
        'k.-');
    plot(axis_lims(1:2), zeros(1,2), 'k-')
    axis(axis_lims)
    title(['perf. ', 'T=', num2str(N_TARGS_GRP(i_grp)), ', S=', num2str(N_SYMBS_GRP(i_grp))]);
end

% plot those from group 4 with perfect memory scores
if ismember(4, group_analysis_list)
    grp_succ_hi_2_perfect = grp_succ_hi_2{4}(:, memory_test_results_grp{4} == 1);
    grp_succ_hi_1_perfect = grp_succ_hi_1{4}(:, memory_test_results_grp{4} == 1);
    subplot(2,4,8); hold on;
    errorbar(.5:.5:4, nanmean(grp_succ_hi_1_perfect, 2), ...
        sqrt(nanvar(grp_succ_hi_1_perfect, [], 2)./size(grp_succ_hi_1_perfect,2)),...
        [cond_colors{1}, '.-']);
    errorbar(.5:.5:4, nanmean(grp_succ_hi_2_perfect, 2), ...
        sqrt(nanvar(grp_succ_hi_2_perfect, [], 2)./size(grp_succ_hi_2_perfect,2)),...
        [cond_colors{2}, '.-']);
    errorbar(.5:.5:4, nanmean(grp_succ_hi_1_perfect - grp_succ_hi_2_perfect, 2),...
        sqrt(nanvar(grp_succ_hi_1_perfect - grp_succ_hi_2_perfect, [], 2)./...
        size(grp_succ_hi_2_perfect,2)),...
        'k.-');
    plot(axis_lims(1:2), zeros(1,2), 'k-')
    axis(axis_lims)
    title(['perf. ', 'T=', num2str(N_TARGS_GRP(i_grp)), ', S=', num2str(N_SYMBS_GRP(i_grp))]);
end
%% compute probability of recall as function of appearances of symbol
symbol_repeat_max = [18, 8, 8, 18];
figure; 
for i_grp = group_analysis_list
    [rec_lpt, rec_hpt] = compute_recall_probability_appearance_order(...
        mov_class_all{i_grp}, target_all{i_grp}, type_all{i_grp},...
        pt_all{i_grp}, min_pt_all(:,i_grp));
    
    lpt_per_targ = reshape(nanmean(rec_lpt, 3), N_SYMBS_GRP(i_grp), size(rec_lpt, 2));
    lpt_all = reshape(nanmean(rec_lpt, 1), size(rec_lpt,2), size(rec_lpt, 3));
    
    subplot(1,4,i_grp); hold on;
    plot(1:size(rec_lpt,2), lpt_per_targ');
    axis([0 symbol_repeat_max(i_grp) 0 1]);
    
    [rec_lpt_direct, rec_hpt_direct] = compute_recall_probability_appearance_order(...
        mov_class_all{i_grp}, target_all{i_grp}, type_all{i_grp},...
        pt_all{i_grp}, min_pt_all(:,i_grp), 1);
    
    lpt_per_targ_direct = reshape(nanmean(rec_lpt_direct, 3), N_SYMBS_GRP(i_grp), size(rec_lpt_direct, 2));
    lpt_all_direct = reshape(nanmean(rec_lpt_direct, 1), size(rec_lpt_direct,2), size(rec_lpt_direct, 3));
    
    subplot(1,4,i_grp); hold on;
    plot(1:size(rec_lpt_direct,2), lpt_per_targ_direct');
    axis([0 symbol_repeat_max(i_grp) 0 1]);
end

figure; 
for i_grp = group_analysis_list
    % for direct cued trials
    [rec_lpt_direct, rec_hpt_direct] = compute_recall_probability_appearance_order(...
        mov_class_all{i_grp}, target_all{i_grp}, type_all{i_grp},...
        pt_all{i_grp}, min_pt_all(:,i_grp), 1);
    
    lpt_per_targ_direct = reshape(nanmean(rec_lpt_direct, 3), N_SYMBS_GRP(i_grp), size(rec_lpt_direct, 2));
    lpt_all_direct = reshape(nanmean(rec_lpt_direct, 1), size(rec_lpt_direct,2), size(rec_lpt_direct, 3));
    
    subplot(1,4,i_grp); hold on;
    errorbar(1:symbol_repeat_max(i_grp), nanmean(lpt_all_direct(1:symbol_repeat_max(i_grp), :),2),...
        sqrt(nanvar(lpt_all_direct(1:symbol_repeat_max(i_grp), :), [], 2)./size(lpt_all_direct, 2)), 'g.-',...
        'MarkerSize', 18, 'LineWidth', 3);
    axis([0 symbol_repeat_max(1) 0 1]);
    
    % for no pre-cued trials
    [rec_lpt_none, rec_hpt_none] = compute_recall_probability_appearance_order(...
        mov_class_all{i_grp}, target_all{i_grp}, type_all{i_grp},...
        pt_all{i_grp}, min_pt_all(:,i_grp), 0);
    
    lpt_per_targ_none = reshape(nanmean(rec_lpt_none, 3), size(rec_lpt_none, 1), size(rec_lpt_none, 2));
    lpt_all_none = reshape(nanmean(rec_lpt_none, 1), size(rec_lpt_none,2), size(rec_lpt_none, 3));
    
    subplot(1,4,i_grp); hold on;
    errorbar(1:symbol_repeat_max(i_grp), nanmean(lpt_all_none(1:symbol_repeat_max(i_grp), :),2),...
        sqrt(nanvar(lpt_all_none(1:symbol_repeat_max(i_grp), :), [], 2)./size(lpt_all_none, 2)), 'r.-',...
        'MarkerSize', 18, 'LineWidth', 3);
    axis([0 symbol_repeat_max(1) 0 1]);
    
    % for symbolicly cued trials
    [rec_lpt, rec_hpt] = compute_recall_probability_appearance_order(...
        mov_class_all{i_grp}, target_all{i_grp}, type_all{i_grp},...
        pt_all{i_grp}, min_pt_all(:,i_grp));
    
    lpt_per_targ = reshape(nanmean(rec_lpt, 3), N_SYMBS_GRP(i_grp), size(rec_lpt, 2));
    lpt_all = reshape(nanmean(rec_lpt, 1), size(rec_lpt,2), size(rec_lpt, 3));
    
    subplot(1,4,i_grp); hold on;
    errorbar(1:symbol_repeat_max(i_grp), nanmean(lpt_all(1:symbol_repeat_max(i_grp), :),2),...
        sqrt(nanvar(lpt_all(1:symbol_repeat_max(i_grp), :), [], 2)./size(lpt_all, 2)), 'b.-',...
        'MarkerSize', 18, 'LineWidth', 3);
    axis([0 symbol_repeat_max(1) 0 1]);
end
%% compute symbol confusion matrix: marginal (i.e. all)
conf_mat_all = {nan(6, 3, length(subject_list_3T)),...
    nan(12, 4, length(subject_list_4T)),...
    nan(12, 6, length(subject_list_6T)),...
    nan(6, 3, length(subject_list_3T_B))};
conf_mat_frac_all = {nan(6, 3, length(subject_list_3T)),...
    nan(12, 4, length(subject_list_4T)),...
    nan(12, 6, length(subject_list_6T)),...
    nan(6, 3, length(subject_list_3T_B))};
ideal_mat_all = {nan(6, 3, length(subject_list_3T)),...
    nan(12, 4, length(subject_list_4T)),...
    nan(12, 6, length(subject_list_6T)),...
    nan(6, 3, length(subject_list_3T_B))};

for i_grp = group_analysis_list
    targs_ = target_all{i_grp};
    mov_class_ = mov_class_all{i_grp};
    type_ = type_all{i_grp};
    pts_ = pt_all{i_grp};
    
    for i_sub = 1:size(targs_,2)
        targs = targs_(:, i_sub);
        type = type_(:, i_sub);
        mov_class = mov_class_(:, i_sub);
        pts = pts_(:, i_sub);
        
        targs = targs(type == 2);
        mov_class = mov_class(type == 2);
        pts = pts(type == 2);
        
        targs_lpt = targs(pts < min_pt_all(i_sub, i_grp));
        targs_hpt = targs(pts >= min_pt_all(i_sub, i_grp));
        move_class_lpt = mov_class(pts < min_pt_all(i_sub, i_grp));
        move_class_hpt = mov_class(pts >= min_pt_all(i_sub, i_grp));
%         hi_pt_cond = pts >= min_pt_all(i_sub, i_grp);
        
        % do for lpt now:
        targs = targs_lpt;
        mov_class = move_class_lpt;
        switch i_grp
            case 1
                targ_syms_1 = targs(targs <= (N_SYMBS_GRP(i_grp)/2));
                targ_syms_2 = targs(targs > (N_SYMBS_GRP(i_grp)/2));
                for i_targ = 1:N_TARGS_GRP(i_grp)
                    targ_syms_2(targ_syms_2 == (N_TARGS_GRP(i_grp) + i_targ)) = i_targ;
                end

                mov_class_1 = mov_class(targs <= (N_SYMBS_GRP(i_grp)/2));
                mov_class_2 = mov_class(targs > (N_SYMBS_GRP(i_grp)/2));

                conf_mat_1 = confusionmat(targ_syms_1, mov_class_1, 'order', 1:N_TARGS_GRP(i_grp));
                conf_mat_2 = confusionmat(targ_syms_2, mov_class_2, 'order', 1:N_TARGS_GRP(i_grp));

                ideal_mat_1 = confusionmat(targ_syms_1, targ_syms_1, 'order', 1:N_TARGS_GRP(i_grp));
                ideal_mat_2 = confusionmat(targ_syms_2, targ_syms_2, 'order', 1:N_TARGS_GRP(i_grp));

                combined_conf_mat = [conf_mat_1(1,:); conf_mat_2(1,:);...
                    conf_mat_1(2,:); conf_mat_2(2,:);...
                    conf_mat_1(3,:); conf_mat_2(3,:)];
                
                ideal_mat_indiv = [ideal_mat_1(1,:); ideal_mat_2(1,:);...
                    ideal_mat_1(2,:); ideal_mat_2(2,:);...
                    ideal_mat_1(3,:); ideal_mat_2(3,:)];
            case 2
                targ_syms_1 = targs(targs <= (N_SYMBS_GRP(i_grp)/3));
                targ_syms_2 = targs(targs > (N_SYMBS_GRP(i_grp)/3) & ...
                    targs <= (2*N_SYMBS_GRP(i_grp)/3));
                targ_syms_3 = targs(targs > (2*N_SYMBS_GRP(i_grp)/3));
                
                for i_targ = 1:N_TARGS_GRP(i_grp)
                    targ_syms_2(targ_syms_2 == (N_TARGS_GRP(i_grp) + i_targ)) = i_targ;
                end
                for i_targ = 1:N_TARGS_GRP(i_grp)
                    targ_syms_3(targ_syms_3 == (2*N_TARGS_GRP(i_grp) + i_targ)) = i_targ;
                end

                mov_class_1 = mov_class(targs <= (N_SYMBS_GRP(i_grp)/3));
                mov_class_2 = mov_class(targs > (N_SYMBS_GRP(i_grp)/3) & ...
                    targs <= (2*N_SYMBS_GRP(i_grp)/3));
                mov_class_3 = mov_class(targs > (2*N_SYMBS_GRP(i_grp)/3));

                conf_mat_1 = confusionmat(targ_syms_1, mov_class_1, 'order', 1:N_TARGS_GRP(i_grp));
                conf_mat_2 = confusionmat(targ_syms_2, mov_class_2, 'order', 1:N_TARGS_GRP(i_grp));
                conf_mat_3 = confusionmat(targ_syms_3, mov_class_3, 'order', 1:N_TARGS_GRP(i_grp));

                ideal_mat_1 = confusionmat(targ_syms_1, targ_syms_1, 'order', 1:N_TARGS_GRP(i_grp));
                ideal_mat_2 = confusionmat(targ_syms_2, targ_syms_2, 'order', 1:N_TARGS_GRP(i_grp));
                ideal_mat_3 = confusionmat(targ_syms_3, targ_syms_3, 'order', 1:N_TARGS_GRP(i_grp));
                
                combined_conf_mat = [conf_mat_1(1,:); conf_mat_2(1,:); conf_mat_3(1,:);...
                    conf_mat_1(2,:); conf_mat_2(2,:); conf_mat_3(2,:);...
                    conf_mat_1(3,:); conf_mat_2(3,:); conf_mat_3(3,:);...
                    conf_mat_1(4,:); conf_mat_2(4,:); conf_mat_3(4,:);...
                    ];
                
                ideal_mat_indiv = [ideal_mat_1(1,:); ideal_mat_2(1,:); ideal_mat_3(1,:);...
                    ideal_mat_1(2,:); ideal_mat_2(2,:); ideal_mat_3(2,:);...
                    ideal_mat_1(3,:); ideal_mat_2(3,:); ideal_mat_3(3,:);...
                    ideal_mat_1(4,:); ideal_mat_2(4,:); ideal_mat_3(4,:)];
            case 3
                targ_syms_1 = targs(targs <= (N_SYMBS_GRP(i_grp)/2));
                targ_syms_2 = targs(targs > (N_SYMBS_GRP(i_grp)/2));
                for i_targ = 1:N_TARGS_GRP(i_grp)
                    targ_syms_2(targ_syms_2 == (N_TARGS_GRP(i_grp) + i_targ)) = i_targ;
                end

                mov_class_1 = mov_class(targs <= (N_SYMBS_GRP(i_grp)/2));
                mov_class_2 = mov_class(targs > (N_SYMBS_GRP(i_grp)/2));

                conf_mat_1 = confusionmat(targ_syms_1, mov_class_1, 'order', 1:N_TARGS_GRP(i_grp));
                conf_mat_2 = confusionmat(targ_syms_2, mov_class_2, 'order', 1:N_TARGS_GRP(i_grp));

                ideal_mat_1 = confusionmat(targ_syms_1, targ_syms_1, 'order', 1:N_TARGS_GRP(i_grp));
                ideal_mat_2 = confusionmat(targ_syms_2, targ_syms_2, 'order', 1:N_TARGS_GRP(i_grp));

                combined_conf_mat = [conf_mat_1(1,:); conf_mat_2(1,:);...
                    conf_mat_1(2,:); conf_mat_2(2,:);...
                    conf_mat_1(3,:); conf_mat_2(3,:);...
                    conf_mat_1(4,:); conf_mat_2(4,:);...
                    conf_mat_1(5,:); conf_mat_2(5,:);...
                    conf_mat_1(6,:); conf_mat_2(6,:)];
                
                ideal_mat_indiv = [ideal_mat_1(1,:); ideal_mat_2(1,:);...
                    ideal_mat_1(2,:); ideal_mat_2(2,:);...
                    ideal_mat_1(3,:); ideal_mat_2(3,:);...
                    ideal_mat_1(4,:); ideal_mat_2(4,:);...
                    ideal_mat_1(5,:); ideal_mat_2(5,:);...
                    ideal_mat_1(6,:); ideal_mat_2(6,:)];
            case 4
                targ_syms_1 = targs(targs <= (N_SYMBS_GRP(i_grp)/2));
                targ_syms_2 = targs(targs > (N_SYMBS_GRP(i_grp)/2));
                for i_targ = 1:N_TARGS_GRP(i_grp)
                    targ_syms_2(targ_syms_2 == (N_TARGS_GRP(i_grp) + i_targ)) = i_targ;
                end

                mov_class_1 = mov_class(targs <= (N_SYMBS_GRP(i_grp)/2));
                mov_class_2 = mov_class(targs > (N_SYMBS_GRP(i_grp)/2));

                conf_mat_1 = confusionmat(targ_syms_1, mov_class_1, 'order', 1:N_TARGS_GRP(i_grp));
                conf_mat_2 = confusionmat(targ_syms_2, mov_class_2, 'order', 1:N_TARGS_GRP(i_grp));

                ideal_mat_1 = confusionmat(targ_syms_1, targ_syms_1, 'order', 1:N_TARGS_GRP(i_grp));
                ideal_mat_2 = confusionmat(targ_syms_2, targ_syms_2, 'order', 1:N_TARGS_GRP(i_grp));

                combined_conf_mat = [conf_mat_1(1,:); conf_mat_2(1,:);...
                    conf_mat_1(2,:); conf_mat_2(2,:);...
                    conf_mat_1(3,:); conf_mat_2(3,:)];
                
                ideal_mat_indiv = [ideal_mat_1(1,:); ideal_mat_2(1,:);...
                    ideal_mat_1(2,:); ideal_mat_2(2,:);...
                    ideal_mat_1(3,:); ideal_mat_2(3,:)];
            otherwise
                error('Illegal group number used.')
        end

        conf_mat_all{i_grp}(:,:,i_sub) = combined_conf_mat;
        
        ideal_mat_all{i_grp}(:,:,i_sub) = ideal_mat_indiv;
        
        for i_symb = 1:size(combined_conf_mat, 1)
            conf_mat_frac_all{i_grp}(i_symb,:,i_sub) = combined_conf_mat(i_symb,:)./sum(combined_conf_mat(i_symb,:),2);
        end
    end
end

%% compute symbol confusion matrix: split into halves:
conf_mat_all_blk = {nan(6, 3, length(subject_list_3T)), nan(6, 3, length(subject_list_3T));...
    nan(12, 4, length(subject_list_4T)), nan(12, 4, length(subject_list_4T));...
    nan(12, 6, length(subject_list_6T)), nan(12, 6, length(subject_list_6T));...
    nan(6, 3, length(subject_list_3T)), nan(6, 3, length(subject_list_3T))};

conf_mat_frac_all_blk = {nan(6, 3, length(subject_list_3T)), nan(6, 3, length(subject_list_3T));...
    nan(12, 4, length(subject_list_4T)), nan(12, 4, length(subject_list_4T));...
    nan(12, 6, length(subject_list_6T)), nan(12, 6, length(subject_list_6T));...
    nan(6, 3, length(subject_list_3T)), nan(6, 3, length(subject_list_3T))};

ideal_mat_all_blk = {nan(6, 3, length(subject_list_3T)), nan(6, 3, length(subject_list_3T));...
    nan(12, 4, length(subject_list_4T)), nan(12, 4, length(subject_list_4T));...
    nan(12, 6, length(subject_list_6T)), nan(12, 6, length(subject_list_6T));...
    nan(6, 3, length(subject_list_3T)), nan(6, 3, length(subject_list_3T))};

tr_blocks = [73, 193, 313]; %trial boundaries between first half and second

for i_grp = group_analysis_list
    for i_block = 1:(length(tr_blocks) - 1)
        targs_ = target_all{i_grp};
        mov_class_ = mov_class_all{i_grp};
        type_ = type_all{i_grp};
        pts_ = pt_all{i_grp};
        
        tr_inds = 1:length(targs_);
        inds = tr_inds(tr_inds > tr_blocks(i_block) & tr_inds < tr_blocks(i_block + 1));
        
        targs_ = targs_(inds, :);
        mov_class_ = mov_class_(inds, :);
        type_ = type_(inds, :);
        pts_ = pts_(inds, :);

        for i_sub = 1:size(targs_,2)
            try
                targs = targs_(:, i_sub);
                type = type_(:, i_sub);
                mov_class = mov_class_(:, i_sub);
                pts = pts_(:, i_sub);

                targs = targs(type == 2);
                mov_class = mov_class(type == 2);
                pts = pts(type == 2);

                targs_lpt = targs(pts < min_pt_all(i_sub, i_grp));
                targs_hpt = targs(pts >= min_pt_all(i_sub, i_grp));
                move_class_lpt = mov_class(pts < min_pt_all(i_sub, i_grp));
                move_class_hpt = mov_class(pts >= min_pt_all(i_sub, i_grp));
        %         hi_pt_cond = pts >= min_pt_all(i_sub, i_grp);

                % do for lpt now:
                targs = targs_lpt;
                mov_class = move_class_lpt;

                switch i_grp
                    case 1
                        targ_syms_1 = targs(targs <= (N_SYMBS_GRP(i_grp)/2));
                        targ_syms_2 = targs(targs > (N_SYMBS_GRP(i_grp)/2));
                        for i_targ = 1:N_TARGS_GRP(i_grp)
                            targ_syms_2(targ_syms_2 == (N_TARGS_GRP(i_grp) + i_targ)) = i_targ;
                        end

                        mov_class_1 = mov_class(targs <= (N_SYMBS_GRP(i_grp)/2));
                        mov_class_2 = mov_class(targs > (N_SYMBS_GRP(i_grp)/2));

                        conf_mat_1 = confusionmat(targ_syms_1, mov_class_1);
                        conf_mat_2 = confusionmat(targ_syms_2, mov_class_2);

                        ideal_mat_1 = confusionmat(targ_syms_1, targ_syms_1);
                        ideal_mat_2 = confusionmat(targ_syms_2, targ_syms_2);

                        combined_conf_mat = [conf_mat_1(1,:); conf_mat_2(1,:);...
                            conf_mat_1(2,:); conf_mat_2(2,:);...
                            conf_mat_1(3,:); conf_mat_2(3,:)];

                        ideal_mat_indiv = [ideal_mat_1(1,:); ideal_mat_2(1,:);...
                            ideal_mat_1(2,:); ideal_mat_2(2,:);...
                            ideal_mat_1(3,:); ideal_mat_2(3,:)];
                    case 2
                        targ_syms_1 = targs(targs <= (N_SYMBS_GRP(i_grp)/3));
                        targ_syms_2 = targs(targs > (N_SYMBS_GRP(i_grp)/3) & ...
                            targs <= (2*N_SYMBS_GRP(i_grp)/3));
                        targ_syms_3 = targs(targs > (2*N_SYMBS_GRP(i_grp)/3));

                        for i_targ = 1:N_TARGS_GRP(i_grp)
                            targ_syms_2(targ_syms_2 == (N_TARGS_GRP(i_grp) + i_targ)) = i_targ;
                        end
                        for i_targ = 1:N_TARGS_GRP(i_grp)
                            targ_syms_3(targ_syms_3 == (2*N_TARGS_GRP(i_grp) + i_targ)) = i_targ;
                        end

                        mov_class_1 = mov_class(targs <= (N_SYMBS_GRP(i_grp)/3));
                        mov_class_2 = mov_class(targs > (N_SYMBS_GRP(i_grp)/3) & ...
                            targs <= (2*N_SYMBS_GRP(i_grp)/3));
                        mov_class_3 = mov_class(targs > (2*N_SYMBS_GRP(i_grp)/3));

                        conf_mat_1 = confusionmat(targ_syms_1, mov_class_1, 'order', 1:N_TARGS_GRP(i_grp));
                        conf_mat_2 = confusionmat(targ_syms_2, mov_class_2, 'order', 1:N_TARGS_GRP(i_grp));
                        conf_mat_3 = confusionmat(targ_syms_3, mov_class_3, 'order', 1:N_TARGS_GRP(i_grp));

                        ideal_mat_1 = confusionmat(targ_syms_1, targ_syms_1, 'order', 1:N_TARGS_GRP(i_grp));
                        ideal_mat_2 = confusionmat(targ_syms_2, targ_syms_2, 'order', 1:N_TARGS_GRP(i_grp));
                        ideal_mat_3 = confusionmat(targ_syms_3, targ_syms_3, 'order', 1:N_TARGS_GRP(i_grp));

                        combined_conf_mat = [conf_mat_1(1,:); conf_mat_2(1,:); conf_mat_3(1,:);...
                            conf_mat_1(2,:); conf_mat_2(2,:); conf_mat_3(2,:);...
                            conf_mat_1(3,:); conf_mat_2(3,:); conf_mat_3(3,:);...
                            conf_mat_1(4,:); conf_mat_2(4,:); conf_mat_3(4,:);...
                            ];

                        ideal_mat_indiv = [ideal_mat_1(1,:); ideal_mat_2(1,:); ideal_mat_3(1,:);...
                            ideal_mat_1(2,:); ideal_mat_2(2,:); ideal_mat_3(2,:);...
                            ideal_mat_1(3,:); ideal_mat_2(3,:); ideal_mat_3(3,:);...
                            ideal_mat_1(4,:); ideal_mat_2(4,:); ideal_mat_3(4,:)];
                    case 3
                        targ_syms_1 = targs(targs <= (N_SYMBS_GRP(i_grp)/2));
                        targ_syms_2 = targs(targs > (N_SYMBS_GRP(i_grp)/2));
                        for i_targ = 1:N_TARGS_GRP(i_grp)
                            targ_syms_2(targ_syms_2 == (N_TARGS_GRP(i_grp) + i_targ)) = i_targ;
                        end

                        mov_class_1 = mov_class(targs <= (N_SYMBS_GRP(i_grp)/2));
                        mov_class_2 = mov_class(targs > (N_SYMBS_GRP(i_grp)/2));

                        conf_mat_1 = confusionmat(targ_syms_1, mov_class_1, 'order', 1:N_TARGS_GRP(i_grp));
                        conf_mat_2 = confusionmat(targ_syms_2, mov_class_2, 'order', 1:N_TARGS_GRP(i_grp));

                        ideal_mat_1 = confusionmat(targ_syms_1, targ_syms_1, 'order', 1:N_TARGS_GRP(i_grp));
                        ideal_mat_2 = confusionmat(targ_syms_2, targ_syms_2, 'order', 1:N_TARGS_GRP(i_grp));

                        combined_conf_mat = [conf_mat_1(1,:); conf_mat_2(1,:);...
                            conf_mat_1(2,:); conf_mat_2(2,:);...
                            conf_mat_1(3,:); conf_mat_2(3,:);...
                            conf_mat_1(4,:); conf_mat_2(4,:);...
                            conf_mat_1(5,:); conf_mat_2(5,:);...
                            conf_mat_1(6,:); conf_mat_2(6,:)];

                        ideal_mat_indiv = [ideal_mat_1(1,:); ideal_mat_2(1,:);...
                            ideal_mat_1(2,:); ideal_mat_2(2,:);...
                            ideal_mat_1(3,:); ideal_mat_2(3,:);...
                            ideal_mat_1(4,:); ideal_mat_2(4,:);...
                            ideal_mat_1(5,:); ideal_mat_2(5,:);...
                            ideal_mat_1(6,:); ideal_mat_2(6,:)];
                    case 4
                        targ_syms_1 = targs(targs <= (N_SYMBS_GRP(i_grp)/2));
                        targ_syms_2 = targs(targs > (N_SYMBS_GRP(i_grp)/2));
                        for i_targ = 1:N_TARGS_GRP(i_grp)
                            targ_syms_2(targ_syms_2 == (N_TARGS_GRP(i_grp) + i_targ)) = i_targ;
                        end

                        mov_class_1 = mov_class(targs <= (N_SYMBS_GRP(i_grp)/2));
                        mov_class_2 = mov_class(targs > (N_SYMBS_GRP(i_grp)/2));

                        conf_mat_1 = confusionmat(targ_syms_1, mov_class_1);
                        conf_mat_2 = confusionmat(targ_syms_2, mov_class_2);

                        ideal_mat_1 = confusionmat(targ_syms_1, targ_syms_1);
                        ideal_mat_2 = confusionmat(targ_syms_2, targ_syms_2);

                        combined_conf_mat = [conf_mat_1(1,:); conf_mat_2(1,:);...
                            conf_mat_1(2,:); conf_mat_2(2,:);...
                            conf_mat_1(3,:); conf_mat_2(3,:)];

                        ideal_mat_indiv = [ideal_mat_1(1,:); ideal_mat_2(1,:);...
                            ideal_mat_1(2,:); ideal_mat_2(2,:);...
                            ideal_mat_1(3,:); ideal_mat_2(3,:)];
                    otherwise
                        error('Illegal group number used.')
                end

                for i_symb = 1:size(combined_conf_mat, 1)
                    conf_mat_frac_all_blk{i_grp, i_block}(i_symb,:,i_sub) = combined_conf_mat(i_symb,:)./sum(combined_conf_mat(i_symb,:),2);
                end

                conf_mat_all_blk{i_grp, i_block}(:,:,i_sub) = combined_conf_mat;
    %             conf_mat_frac_all_blk{i_grp, i_block} = conf_mat_frac_all;
                ideal_mat_all_blk{i_grp, i_block}(:,:,i_sub) = ideal_mat_indiv;
            catch sub_err
                warning('Subject failed to compute confusion matrix');
            end
        end
    end
end


%% compare mPT across experiments (which had differing numbers of targets)
figure;
min_pt_corr = min_pt_all - .06;
errorbar([3 4 6 3], nanmean(min_pt_corr,1), sqrt(nanvar(min_pt_corr)./sum(~isnan(min_pt_corr))), 'ks');
axis([2.5 6.5 .2 .25])

%% compute pr(succ | pt) for first half and second half of blocks for each group:

inds_init = 1:(12+60);
inds_half_1 = inds_init(end) + (1:120);
inds_half_2 = inds_half_1(end) + (1:120);

figure;
sub_plot_inds = [1 2; 3 4; 5 6; 7 8];
hist_bins = linspace(EARLIEST_VALID_PT, LATEST_VALID_PT, n_bins+1);
x_inds = hist_bins(1:(end-1)) + diff(hist_bins)/2;
for i_grp = group_analysis_list
    grp_pr_succ_0 = nan(n_bins, size(pt_all{i_grp}, 2));
    grp_pr_succ_half1 = nan(n_bins, size(pt_all{i_grp}, 2), 3);
    grp_pr_succ_half2 = nan(n_bins, size(pt_all{i_grp}, 2), 3);
    for i_sub = 1:size(pt_all{i_grp},2)
        [grp_pr_succ_0(:, i_sub), ~, ~] = ...
         compute_succ_prob_across_pt(...
         pt_all{i_grp}(inds_init, i_sub),...
         dir_err_all{i_grp}(inds_init, i_sub),...
         type_all{i_grp}(inds_init, i_sub));
        
        [grp_pr_succ_half1(:, i_sub, 1),...
            grp_pr_succ_half1(:, i_sub, 2),... 
            grp_pr_succ_half1(:, i_sub, 3)] = ...
         compute_succ_prob_across_pt(...
            pt_all{i_grp}(inds_half_1, i_sub), ...
            dir_err_all{i_grp}(inds_half_1, i_sub), ...
            type_all{i_grp}(inds_half_1, i_sub));
        
        [grp_pr_succ_half2(:, i_sub, 1),...
            grp_pr_succ_half2(:, i_sub, 2),...
            grp_pr_succ_half2(:, i_sub, 3)] = ...
         compute_succ_prob_across_pt(...
            pt_all{i_grp}(inds_half_2, i_sub), ...
            dir_err_all{i_grp}(inds_half_2, i_sub), ...
            type_all{i_grp}(inds_half_2, i_sub));
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
% Look at Pr(succ), response to catch trials, and variability; for first
% and second half of blocks, for each group

inds_init = 1:(12+60);
inds_q1 = inds_init(end) + (1:60);
inds_q2 = inds_q1(end) + (1:60);
inds_q3 = inds_q2(end) + (1:60);
inds_q4 = inds_q3(end) + (1:60);

% Plot directional error vs. PT, aligned to minPT
figure;
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
end

%%
% Plot success vs. PT and var. vs. PT, aligned to minPT
h1 = figure;
h2 = figure;
% sub_plot_inds = [1 2 3 4; 5 6 7 8; 9 10 11 12; 13 14 15 16];
sub_plot_inds = [1 2 3 4 5; 6 7 8 9 10; 11 12 13 14 15; 16 17 18 19 20];
% override globals:
EARLIEST_VALID_PT = -0.5; LATEST_VALID_PT = 0.5; n_bins = 12;
hist_bins = linspace(EARLIEST_VALID_PT, LATEST_VALID_PT, n_bins+1);
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
        pt_init(:), de_init(:), type_init(:));
    
    var_inds_init = abs(de_init(:)) < SUCCESS_TH_ANGLE;
    [vr0, ~, ~] = compute_var_across_pt(...
        pt_init(var_inds_init),...
        de_init(var_inds_init), type_init(var_inds_init));
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
    
    
    [~, ps1, ps2] = compute_succ_prob_across_pt(...
        pt_q1(:), de_q1(:), type_q1(:));
    
    var_inds_half1 = abs(de_q1(:)) < SUCCESS_TH_ANGLE;
    [~, vr1, vr2] = compute_var_across_pt(...
        pt_q1(var_inds_half1),...
        de_q1(var_inds_half1), type_q1(var_inds_half1));    
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
        pt_q2(:), de_q2(:), type_q2(:));
    
    var_inds_half2 = abs(de_q2(:)) < SUCCESS_TH_ANGLE;
    [~, vr1, vr2] = compute_var_across_pt(...
        pt_q2(var_inds_half2),...
        de_q2(var_inds_half2), type_q2(var_inds_half2));
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
        pt_q3(:), de_q3(:), type_q3(:));
    
    var_inds_half2 = abs(de_q3(:)) < SUCCESS_TH_ANGLE;
    [~, vr1, vr2] = compute_var_across_pt(...
        pt_q3(var_inds_half2),...
        de_q3(var_inds_half2), type_q3(var_inds_half2));
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
        pt_q4(:), de_q4(:), type_q4(:));
    
    var_inds_half2 = abs(de_q4(:)) < SUCCESS_TH_ANGLE;
    [~, vr1, vr2] = compute_var_across_pt(...
        pt_q4(var_inds_half2),...
        de_q4(var_inds_half2), type_q4(var_inds_half2));
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
% % re-set globals:
% EARLIEST_VALID_PT = -0.2; LATEST_VALID_PT = 0.8; n_bins = 6;

%%
% Plot directional error vs. PT of catch trials , aligned to minPT
figure;
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