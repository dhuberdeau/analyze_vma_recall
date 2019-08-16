%% handle global variables
global EARLIEST_VALID_PT LATEST_VALID_PT SUCCESS_TH_ANGLE...
     RED_COLOR GREEN_COLOR BLUE_COLOR RELATIVE_PT_MINMAX

%% list subject data files:
group_subjects = {...
'Data_S027_03262018_E1.mat'...
'Data_S028_03262018_E1.mat'...
'Data_S029_03272018_E1.mat'...
'Data_S030_03272018_E1.mat'...
'Data_S031_03272018_E1.mat'...
'Data_S032_03272018_E1.mat'...
'Data_S033_03272018_E1.mat'...
'Data_S034_03272018_E1.mat'...
'Data_S035_03272018_E1.mat'...
...% 'Data_S036_03272018_E1.mat'...
'Data_S037_03272018_E1.mat'...
'Data_S038_03282018_E1.mat'...
'Data_S039_03282018_E1.mat'...
'Data_S040_03282018_E1.mat'...
'Data_S041_03282018_E1.mat'...
'Data_S042_03282018_E1.mat'...
'Data_S043_04042018_E1.mat'...
'Data_S044_04042018_E1.mat'...
'Data_S045_04042018_E1.mat'...
'Data_S071_04092018_E1.mat'...
'Data_S073_04102018_E1.mat'...
};

%% define some constants and parameters:
% EARLIEST_VALID_PT = -.1;
% LATEST_VALID_PT = .7;
% SUCCESS_TH_ANGLE = 30; % +/- 30-degrees around 0 directional error
% MIN_N_PT_FOR_MEASURE = 4;
% RED_COLOR = [172, 59, 59]/255;
% GREEN_COLOR = [85, 170, 85]/255; 
% BLUE_COLOR = [86 85 149]/255;

RELATIVE_PT_MINMAX = 0.5;
PT_BINS = 7;
%% Analyze each subject:

indiv_error_queue = {};

a = [];
all_timer = tic;
p_bins = nan(8, length(group_subjects), 3);
z_bins = nan(6, length(group_subjects), 2);
view_time = nan(70, length(group_subjects), 3);
min_pt = nan(1, length(group_subjects));
dir_error = nan(70, length(group_subjects), 3);
dir_absolute = nan(70, length(group_subjects), 3);
catch_tr_apt = nan(12, length(group_subjects), 2);
catch_tr_ppt = nan(12, length(group_subjects), 2);
catch_tr_de = nan(12, length(group_subjects), 2);
kin_x_catch3 = nan(25, 12, length(group_subjects));
kin_y_catch3 = nan(25, 12, length(group_subjects));
kin_x_catch4 = nan(25, 12, length(group_subjects));
kin_y_catch4 = nan(25, 12, length(group_subjects));
% kinematics_all = cell(2, length(group_subjects));

% setup matricies to collect all data:
N_trials = 192;
type_all = nan(N_trials, length(group_subjects));
pt_all = nan(N_trials, length(group_subjects));
move_all = nan(N_trials, length(group_subjects));
target_all = nan(N_trials, length(group_subjects));
direrror_all = nan(N_trials, length(group_subjects));
dir_rel_error_all = nan(N_trials, length(group_subjects)); % error relative to chosen target (even if wrong)
kin_PV = nan(N_trials, length(group_subjects));
kin_MT = nan(N_trials, length(group_subjects));
kin_VAR = nan(N_trials, length(group_subjects));


targ_angles = 0:90:300;
targ_angles(targ_angles > 180) = targ_angles(targ_angles > 180) - 360;
targ_angles_plus = [targ_angles, -180]; % 0 and -180 are same target

example_subject = 7;

trial_times = nan(length(group_subjects), 192);

h1 = figure; h2 = figure; h3 = figure; h4 = figure;
target_distances = nan(1, length(group_subjects));
for i_sub = 1:length(group_subjects)
    sub_timer = tic;
    try
        load(group_subjects{i_sub})
        data_indiv = analyze_retention_individual_v1(Data, 0);
        
        for i_tr = 1:size(trial_times,2)
            trial_times(i_sub, i_tr) = Data.Kinematics{i_tr}(end,1) - Data.Kinematics{i_tr}(1,1);
        end
        
        p_bin0 = data_indiv.p_bins{1};
        p_bin1 = data_indiv.p_bins{2};
        p_bin2 = data_indiv.p_bins{3};
        pt_z_bin0 = data_indiv.z_bins{1};
        pt_z_bin2 = data_indiv.z_bins{2};
        
        p_bins(:, i_sub, 1) = p_bin0;
        p_bins(:, i_sub, 2) = p_bin1;
        p_bins(:, i_sub, 3) = p_bin2;

        z_bins(:, i_sub, 1) = pt_z_bin0;
        z_bins(:, i_sub, 2) = pt_z_bin2;

        min_pt(i_sub) = data_indiv.min_pt;
        
        type0 = data_indiv.types{1};
        type1 = data_indiv.types{2};
        type2 = data_indiv.types{3};
        type3 = data_indiv.types{4};
        type4 = data_indiv.types{5};
        
        % assign data to appropriate indicies:
        type_all(:, i_sub) = Data.Type;
        pt_all(:, i_sub) = data_indiv.ViewTime;
        target_all(:, i_sub) = Data.Target;
        for i_tr = 1:length(data_indiv.Dir_a)
            [targ_off_mag, targ_class] = min(abs(data_indiv.Dir_a(i_tr) - targ_angles));
            if abs(targ_off_mag) <= SUCCESS_TH_ANGLE
                move_all(i_tr, i_sub) = targ_class;
            end
        end
        
        Dir_e = data_indiv.Dir_e;
        Dir_a = data_indiv.Dir_a;
        
        Kin_x = data_indiv.kinematics{1};
        Kin_y = data_indiv.kinematics{2};
        
        kin_PV(1:length(data_indiv.kin_summary{1}), i_sub) = data_indiv.kin_summary{1};
        kin_MT(1:length(data_indiv.kin_summary{2}), i_sub) = data_indiv.kin_summary{2};
        kin_VAR(1:length(data_indiv.kin_summary{3}), i_sub) = data_indiv.kin_summary{3};
        
        direrror_all(1:length(Dir_e), i_sub) = Dir_e;
        dir_rel_error_all(1:length(Dir_a), i_sub) = ...
            compute_relative_target_error(Dir_a, targ_angles_plus);
        
        for i_traj = 1:length(type3)
            kin_x_catch3(:, i_traj, i_sub) = Kin_x(:, type3(i_traj));
            kin_y_catch3(:, i_traj, i_sub) = Kin_y(:, type3(i_traj));
        end
        for i_traj = 1:length(type4)
            kin_x_catch4(:, i_traj, i_sub) = Kin_x(:, type4(i_traj));
            kin_y_catch4(:, i_traj, i_sub) = Kin_y(:, type4(i_traj));
        end
        
        vt_temp = data_indiv.ViewTime(Data.Type == 0);
        de_temp = Dir_e(Data.Type == 0);
        da_temp = Dir_a(Data.Type == 0);
        view_time(1:length(vt_temp), i_sub, 1) = vt_temp;
        dir_error(1:length(de_temp), i_sub, 1) = de_temp;
        dir_absolute(1:length(da_temp), i_sub, 1) = da_temp;

        vt_temp = data_indiv.ViewTime(Data.Type == 1);
        de_temp = Dir_e(Data.Type == 1);
        da_temp = Dir_a(Data.Type == 1);
        view_time(1:length(vt_temp), i_sub, 2) = vt_temp;
        dir_error(1:length(de_temp), i_sub, 2) = de_temp;
        dir_absolute(1:length(da_temp), i_sub, 2) = da_temp;

        vt_temp = data_indiv.ViewTime(Data.Type == 2);
        de_temp = Dir_e(Data.Type == 2);
        da_temp = Dir_a(Data.Type == 2);
        view_time(1:length(vt_temp), i_sub, 3) = vt_temp;
        dir_error(1:length(de_temp), i_sub, 3) = de_temp;
        dir_absolute(1:length(da_temp), i_sub, 3) = da_temp;
        
        catch_tr_apt(1:length(type3), i_sub, 1) = data_indiv.ViewTime(type3);
        catch_tr_apt(1:length(type4), i_sub, 2) = data_indiv.ViewTime(type4);
        
        catch_tr_ppt(1:length(type3), i_sub, 1) = Data.pPT(type3);
        catch_tr_ppt(1:length(type4), i_sub, 2) = Data.pPT(type4);
        
        catch_tr_de(1:length(type3), i_sub, 1) = Dir_e(type3);
        catch_tr_de(1:length(type4), i_sub, 2) = Dir_e(type4);

        open_winds = findobj('type', 'figure');
        for i_wind = 1:(length(open_winds)-4)
            close(i_wind + 2);
        end
        figure(h1); 
        subplot(ceil(sqrt(length(group_subjects))), ceil(sqrt(length(group_subjects))), i_sub); hold on;
        plot(data_indiv.ViewTime(type0), Dir_e(type0) , '.', 'MarkerSize', 16, 'Color', [172, 59, 59]/255)
        plot(data_indiv.ViewTime(type1), Dir_e(type1), '.', 'MarkerSize', 15, 'Color', [85, 170, 85]/255)
        plot(data_indiv.ViewTime(type2), Dir_e(type2), '.', 'MarkerSize', 14, 'Color', [86/255 85/255 149/255])
        plot(data_indiv.ViewTime(type3), Dir_e(type3), 'o', 'MarkerSize', 8, 'LineWidth', 2, 'Color', [225, 244, 162]/255)
        plot(data_indiv.ViewTime(type4), Dir_e(type4), 'o', 'MarkerSize', 8, 'LineWidth', 2, 'Color', [71, 113, 134]/255)
        plot([-.2 .8], [45 45], 'k')
        plot([-.2 .8], -[45 45], 'k')
        plot([min_pt(i_sub), min_pt(i_sub)], [-200 200], 'k-');
        axis([EARLIEST_VALID_PT, LATEST_VALID_PT, -200 200])
        
        figure(h2);
        subplot(ceil(sqrt(length(group_subjects))), ceil(sqrt(length(group_subjects))), i_sub); hold on;
        plot(Kin_x(:,type0), -Kin_y(:,type0), '-', 'LineWidth', 2, 'Color', [172, 59, 59]/255)
        plot(Kin_x(:,type1), -Kin_y(:,type1), '-', 'LineWidth', 1, 'Color', [85, 170, 85]/255)
        plot(Kin_x(:,type2), -Kin_y(:,type2), '-', 'LineWidth', .5, 'Color', [86/255 85/255 149/255])
        axis([-20 20 -20 20])
        
        if i_sub == example_subject
            figure(h3); hold on;
            plot(data_indiv.ViewTime(type0), Dir_e(type0) , '.', 'MarkerSize', 16, 'Color', [172, 59, 59]/255)
            plot(data_indiv.ViewTime(type1), Dir_e(type1), '.', 'MarkerSize', 15, 'Color', [85, 170, 85]/255)
            plot(data_indiv.ViewTime(type2), Dir_e(type2), '.', 'MarkerSize', 14, 'Color', [86/255 85/255 149/255])
            set(gca,'fontsize',20)
            axis([EARLIEST_VALID_PT, LATEST_VALID_PT, -200 200])
            
            figure(h4); hold on;
            plot(Kin_x(:,type0), -Kin_y(:,type0), '-', 'LineWidth', 2, 'Color', [172, 59, 59]/255)
            plot(Kin_x(:,type1), -Kin_y(:,type1), '-', 'LineWidth', 1, 'Color', [85, 170, 85]/255)
            plot(Kin_x(:,type2), -Kin_y(:,type2), '-', 'LineWidth', .5, 'Color', [86/255 85/255 149/255])
            axis([-20 20 -20 20])
        end
        
        target_distances(i_sub) = data_indiv.targ_len;
    catch individ_err
        warning(['Subject ', num2str(i_sub), ' failed.'])
        indiv_error_queue{length(indiv_error_queue) + 1} = individ_err;
    end
    a(i_sub) = toc(sub_timer);
end
toc(all_timer)

%% Test for differences among trial types on kinematic summary variables: 
% which include Movement time (MT), peak velocity (PV), and movement
% variability (VAR)
test_kin_summary_variables;

%% plot group results: pr correct.
table_pc = plot_behavior_over_pt(pt_all, abs(direrror_all) < SUCCESS_TH_ANGLE, type_all, PT_BINS);
f_ = gcf;
xlabel('Preparation Time (sec)')
ylabel('Probability correct');
set(f_, 'Position', [0 0 600 200])
saveas(f_, 'Probability_correct.pdf');

%% plot group results: absoluate direction error (de).
table_ade = plot_behavior_over_pt(pt_all, abs(direrror_all), type_all, PT_BINS);
f_ = gcf;
xlabel('Preparation Time (sec)')
ylabel('Directional error (degrees)');
set(f_, 'Position', [0 200 600 200])
saveas(f_, 'Direction_error.pdf');

%% plot group results: movement time (MT)
table_mt = plot_behavior_over_pt(pt_all, kin_MT, type_all, PT_BINS);
f_ = gcf;
xlabel('Preparation Time (sec)')
ylabel('Movement time (sec)');
set(f_, 'Position', [0 400 600 200])
saveas(f_, 'Movement_time.pdf');

%% plot group results: peak velocity (PV)
table_pv = plot_behavior_over_pt(pt_all, kin_PV, type_all, PT_BINS);
f_ = gcf;
xlabel('Preparation Time (sec)')
ylabel('Peak Velocity (cm/sec)');
set(f_, 'Position', [0 600 600 200])
saveas(f_, 'Peak_velocity.pdf');

%% plot group results: movement variability (computed from rde, relative direction error)
% table_mv = plot_behavior_over_pt(pt_all, kin_VAR, type_all, PT_BINS);
table_rde = plot_behavior_over_pt(pt_all, dir_rel_error_all, type_all, PT_BINS);

% plot mean variability (SD) of movement initiation directions binned by PT.
figure; hold on;
errorbar(nanmean(table_rde{1}) + 0.02, nanmean(table_rde{3}(:, :, 1)),...
sqrt(nanvar(table_rde{3}(:, :, 1))./sum(~isnan(table_rde{3}(:, :, 1)))), ...
'.-', 'LineWidth', 2, 'MarkerSize', 20, 'Color', RED_COLOR);

errorbar(nanmean(table_rde{1}) - 0.02, nanmean(table_rde{3}(:, :, 2)),...
sqrt(nanvar(table_rde{3}(:, :, 2))./sum(~isnan(table_rde{3}(:, :, 2)))), ...
'.-', 'LineWidth', 2, 'MarkerSize', 20, 'Color', GREEN_COLOR);

errorbar(nanmean(table_rde{1}), nanmean(table_rde{3}(:, :, 3)),...
sqrt(nanvar(table_rde{3}(:, :, 3))./sum(~isnan(table_rde{3}(:, :, 3)))), ...
'.-', 'LineWidth', 2, 'MarkerSize', 20, 'Color', BLUE_COLOR);

% errorbar(nanmean(table_rde{1}), nanmean(table_rde{3}(:, :, 4)),...
% sqrt(nanvar(table_rde{3}(:, :, 4))./sum(~isnan(table_rde{3}(:, :, 4)))), ...
% 'o-',  'LineWidth', 2, 'MarkerSize', 12, 'Color', GREEN_COLOR);
% 
% errorbar(nanmean(table_rde{1}), nanmean(table_rde{3}(:, :, 5)),...
% sqrt(nanvar(table_rde{3}(:, :, 5))./sum(~isnan(table_rde{3}(:, :, 5)))), ...
% 'o-',  'LineWidth', 2, 'MarkerSize', 12, 'Color', BLUE_COLOR);

axis([EARLIEST_VALID_PT, LATEST_VALID_PT, min(min(min(table_rde{3}))), max(max(max(table_rde{3})))])
set(gca,'fontsize',18)
legend('None', 'Direct', 'Symbolic', 'location', 'bestoutside')

f_ = gcf;
xlabel('Preparation Time (sec.)')
ylabel('Directional variability (degrees)');
set(f_, 'Position', [0 800 600 200])
saveas(f_, 'Direction_variability.pdf');

%% export relevant variables for analysis in R:

% export summary measures:
table_out_pc = make_data_table(table_pc);
table_out_var = make_data_table(table_rde, 3);
table_out_pv = make_data_table(table_pv);
table_out_err = make_data_table(table_ade);

csvwrite('table_E1_pc', table_out_pc);
csvwrite('table_E1_var', table_out_var);
csvwrite('table_E1_pv', table_out_pv);
csvwrite('table_E1_de', table_out_err);

% export raw measures:
temp_pc = abs(direrror_all(:)) < 30;
temp_type = type_all(:);
temp_pt_diff_ = pt_all - repmat(min_pt, size(pt_all,1), 1);
temp_pt_diff = temp_pt_diff_(:);
temp_sub = reshape(repmat((1:size(pt_all,2)), size(pt_all,1), 1), numel(pt_all), 1);
inds_valid_pt = temp_pt_diff > -RELATIVE_PT_MINMAX & temp_pt_diff < RELATIVE_PT_MINMAX;

data_mat_pc = [temp_pc(inds_valid_pt),...
    temp_type(inds_valid_pt),...
    temp_pt_diff(inds_valid_pt),...
    temp_sub(inds_valid_pt)];

csvwrite('raw_data_mat_E1_pc', data_mat_pc);

%% Align data to each participant's individual minPT and compute Pr(corr) for catch trials
view_time_all_0 = reshape(view_time(:, :, 1) - repmat(min_pt, size(view_time,1), 1),...
    size(view_time,1)*size(view_time,2), 1);
view_time_all_1 = reshape(view_time(:, :, 2) - repmat(min_pt, size(view_time,1), 1),...
    size(view_time,1)*size(view_time,2), 1);
view_time_all_2 = reshape(view_time(:, :, 3) - repmat(min_pt, size(view_time,1), 1),...
    size(view_time,1)*size(view_time,2), 1);

de_all_0 = reshape(dir_error(:, :, 1), size(dir_error,1)*size(dir_error,2), 1);
de_all_1 = reshape(dir_error(:, :, 2), size(dir_error,1)*size(dir_error,2), 1);
de_all_2 = reshape(dir_error(:, :, 3), size(dir_error,1)*size(dir_error,2), 1);

view_time_all_3 = reshape(catch_tr_apt(:, :, 1) - repmat(min_pt, size(catch_tr_apt,1), 1),...
    size(catch_tr_apt,1)*size(catch_tr_apt,2), 1);
view_time_all_4 = reshape(catch_tr_apt(:, :, 2) - repmat(min_pt, size(catch_tr_apt,1), 1),...
    size(catch_tr_apt,1)*size(catch_tr_apt,2), 1);

de_all_3 = reshape(catch_tr_de(:, :, 1), size(catch_tr_de,1)*size(catch_tr_de,2), 1);
de_all_4 = reshape(catch_tr_de(:, :, 2), size(catch_tr_de,1)*size(catch_tr_de,2), 1);

n_bins = 12;
hist_bins_aligned = linspace(-0.5, 0.5, n_bins + 1);
[n_pt_all, edge_pt_all] = histcounts(de_all_0, hist_bins_aligned);
x_ind = edge_pt_all(2:end) - diff(edge_pt_all)/2;

% setup probability of correct response per bin:
p_all_0 = nan(length(edge_pt_all) - 1, 1);
p_all_1 = nan(length(edge_pt_all) - 1, 1);
p_all_2 = nan(length(edge_pt_all) - 1, 1);
p_all_3 = nan(length(edge_pt_all) - 1, 1);
p_all_4 = nan(length(edge_pt_all) - 1, 1);

% setup variability of launch direction per bin:
var_all_0 = nan(length(edge_pt_all) - 1, 1);
var_all_1 = nan(length(edge_pt_all) - 1, 1);
var_all_2 = nan(length(edge_pt_all) - 1, 1);

for i_bin = 1:(length(edge_pt_all) - 1)
    inds_0 = view_time_all_0 >= edge_pt_all(i_bin) &...
        view_time_all_0 < edge_pt_all(i_bin + 1);
    inds_1 = view_time_all_1 >= edge_pt_all(i_bin) &...
        view_time_all_1 < edge_pt_all(i_bin + 1);
    inds_2 = view_time_all_2 >= edge_pt_all(i_bin) &...
        view_time_all_2 < edge_pt_all(i_bin + 1);
    inds_3 = view_time_all_3 >= edge_pt_all(i_bin) &...
        view_time_all_3 < edge_pt_all(i_bin + 1);
    inds_4 = view_time_all_4 >= edge_pt_all(i_bin) &...
        view_time_all_4 < edge_pt_all(i_bin + 1);
    
    % type 0 (no-precue)
    inds_0_valid = inds_0 & ~isnan(de_all_0);
    p_all_0(i_bin) = sum(abs(de_all_0(inds_0_valid)) < SUCCESS_TH_ANGLE)/sum(inds_0_valid);
    inds_0_var = inds_0_valid & abs(de_all_0) < SUCCESS_TH_ANGLE;
    var_all_0(i_bin) = sqrt(var(de_all_0(inds_0_var)));
    
    % type 1 (direct-precue)
    inds_1_valid = inds_1 & ~isnan(de_all_1);
    p_all_1(i_bin) = sum(abs(de_all_1(inds_1_valid)) < SUCCESS_TH_ANGLE)/sum(inds_1_valid);
    inds_1_var = inds_1_valid & abs(de_all_1) < SUCCESS_TH_ANGLE;
    var_all_1(i_bin) = sqrt(var(de_all_1(inds_1_var)));
    
    % type 2 (symbolic-precue)
    inds_2_valid = inds_2 & ~isnan(de_all_2);
    p_all_2(i_bin) = sum(abs(de_all_2(inds_2_valid)) < SUCCESS_TH_ANGLE)/sum(inds_2_valid);
    inds_2_var = inds_2_valid & abs(de_all_2) < SUCCESS_TH_ANGLE;
    var_all_2(i_bin) = sqrt(var(de_all_2(inds_2_var)));
    
    % type 3 (catch trial direct)
    inds_3_valid = inds_3 & ~isnan(de_all_3);
    p_all_3(i_bin) = sum(abs(de_all_3(inds_3_valid)) < SUCCESS_TH_ANGLE)/sum(inds_3_valid);
    
    % type 2 (catch trial symbolic)
    inds_4_valid = inds_4 & ~isnan(de_all_4);
    p_all_4(i_bin) = sum(abs(de_all_4(inds_4_valid)) < SUCCESS_TH_ANGLE)/sum(inds_4_valid);
    
    
end

%% plot probability of recall as function of appearance of symbol:
[rec_lpt, rec_hpt] = compute_recall_probability_appearance_order(...
    move_all, target_all, type_all,...
    pt_all, min_pt, 2);

max_lpt = 15;

lpt_per_targ = reshape(nanmean(rec_lpt, 3), 4, size(rec_lpt, 2));
lpt_all = reshape(nanmean(rec_lpt, 1), size(rec_lpt,2), size(rec_lpt, 3));

figure; hold on;
    errorbar((1:max_lpt) - 0.1, nanmean(lpt_all(1:max_lpt, :),2),...
        sqrt(nanvar(lpt_all(1:max_lpt, :), [], 2)./size(lpt_all, 2)), 'b.-',...
        'Color', [86 85 149]/255, 'MarkerSize', 18, 'LineWidth', 2);
axis([0 max_lpt 0 1]);

[rec_lpt_direct, rec_hpt_direct] = compute_recall_probability_appearance_order(...
    move_all, target_all, type_all,...
    pt_all, min_pt, 1);

lpt_per_targ_direct = reshape(nanmean(rec_lpt_direct, 3), 4, size(rec_lpt_direct, 2));
lpt_all_direct = reshape(nanmean(rec_lpt_direct, 1), size(rec_lpt_direct,2), size(rec_lpt_direct, 3));

errorbar((1:max_lpt) + 0.1 , nanmean(lpt_all_direct(1:max_lpt, :),2),...
        sqrt(nanvar(lpt_all_direct(1:max_lpt, :), [], 2)./size(lpt_all_direct, 2)), 'g.-',...
        'Color', [85, 170, 85]/255, 'MarkerSize', 18, 'LineWidth', 2);
axis([0, max_lpt + 0.1, 0, 1.1]);
set(gca,'fontsize',18)
xlabel('Symbol occurance');
ylabel('Probability correct');
f_ = gcf;
set(f_, 'Position', [0 0 300 420])
saveas(f_, 'Recall_probability.pdf');



