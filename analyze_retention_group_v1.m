%% handle global variables
global EARLIEST_VALID_PT LATEST_VALID_PT SUCCESS_TH_ANGLE

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
EARLIEST_VALID_PT = -.2;
LATEST_VALID_PT = .85;
SUCCESS_TH_ANGLE = 30; % +/- 30-degrees around 0 directional error
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

targ_angles = 0:90:300;
targ_angles(targ_angles > 180) = targ_angles(targ_angles > 180) - 360;

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

%% plot group aggregate results:
n_bins = 8;
vt = view_time(:);
hist_bins = linspace(EARLIEST_VALID_PT, LATEST_VALID_PT, n_bins+1);
[n_pt_all, edge_pt_all] = histcounts(vt, hist_bins);
x_ind = edge_pt_all(2:end) - diff(edge_pt_all)/2;

figure; hold on;
errorbar(x_ind + 0.02, nanmean(p_bins(:,:,1), 2), sqrt(nanvar(p_bins(:,:,1),0,2)./sum(~isnan(p_bins(:,:,1)),2)), '.-', 'MarkerSize', 20, 'Color', [172, 59, 59]/255);
errorbar(x_ind - 0.02, nanmean(p_bins(:,:,2), 2), sqrt(nanvar(p_bins(:,:,2),0,2)./sum(~isnan(p_bins(:,:,2)),2)), '.-', 'MarkerSize', 20, 'Color', [85, 170, 85]/255);
errorbar(x_ind, nanmean(p_bins(:,:,3), 2), sqrt(nanvar(p_bins(:,:,3),0,2)./sum(~isnan(p_bins(:,:,3)),2)), '.-', 'MarkerSize', 20, 'Color', [86/255 85/255 149/255]);
axis([EARLIEST_VALID_PT, LATEST_VALID_PT, 0 1])
set(gca,'fontsize',20)
legend('None', 'Direct', 'Symbolic')

n_bins = 6; hist_bins = linspace(-.005, .605, n_bins+1);
[pt_sort, i_sort] = sort(Data.pPT);
[n_ppt_all, edge_ppt_all] = histcounts(pt_sort, hist_bins);
x_ind = edge_ppt_all(2:end) - diff(edge_ppt_all)/2;
figure; hold on;
errorbar(x_ind + 0.01, nanmean(z_bins(:,:,1), 2), sqrt(nanvar(z_bins(:,:,1),0,2)./length(group_subjects)), '.-', 'MarkerSize', 20, 'Color', [172, 59, 59]/255);
errorbar(x_ind - 0.01, nanmean(z_bins(:,:,2), 2), sqrt(nanvar(z_bins(:,:,2),0,2)./length(group_subjects)), '.-', 'MarkerSize', 20, 'Color', [86/255 85/255 149/255]);
set(gca,'fontsize',20)
legend('None', 'Symbolic')

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

%% plot dir. err. aligned to minPT, and catch trial err., and pr(corr), and variability
figure;
subplot(4,1,1); hold on;
plot([0 0], [-200 200], 'k-')
plot(view_time_all_0, de_all_0, 'r.') %type 0
plot(view_time_all_1, de_all_1, 'g.') %type 1
plot(view_time_all_2, de_all_2, 'b.') %type 2
plot([-.5 .5], [30 30], '-', 'Color', [.5 .5 .5]);
plot([-.5 .5], -[30 30], '-', 'Color', [.5 .5 .5]);
axis([-0.5 0.5 -200 200]);

subplot(4,1,2); hold on;
plot([0 0], [-200 200], 'k-')
plot(view_time_all_3, de_all_3, 'go') %type 3 (catch trials, direct cue)
plot(view_time_all_4, de_all_4, 'bo') %type 4 (catch trials, symbol cue) 
plot([-.5 .5], [30 30], '-', 'Color', [.5 .5 .5]);
plot([-.5 .5], -[30 30], '-', 'Color', [.5 .5 .5]);
axis([-0.5 0.5 -200 200]);

subplot(4,1,3); hold on;
plot([0 0], [0 1], 'k-')
plot(x_ind, p_all_0, 'r.-');
plot(x_ind, p_all_1, 'g.-');
plot(x_ind, p_all_2, 'b.-');
plot(x_ind, p_all_3, 'go-');
plot(x_ind, p_all_4, 'bo-');
axis([-0.5 0.5 0 1]);
plot([-.5 .5], [.25 .25], '-', 'Color', [.5 .5 .5]);
plot([-.5 .5], [30 30], '--', 'Color', [.5 .5 .5]);
plot([-.5 .5], -[30 30], '--', 'Color', [.5 .5 .5]);

subplot(4,1,4); hold on;
plot([0 0], [0 12], 'k-')
plot(x_ind, var_all_0, 'r.-');
plot(x_ind, var_all_1, 'g.-');
plot(x_ind, var_all_2, 'b.-');

%% Plot variability below min pt and above min pt for each type, avg.ed across people..
[reach_var_persub_0, reach_var_persub_1, reach_var_persub_2] = ...
    compute_variability_by_pt(view_time, dir_error, min_pt);
% reach_var_persub_0 = nan(size(view_time,2), 2); %subject x low-PT or high-PT
% reach_var_persub_1 = nan(size(view_time,2), 2); %subject x low-PT or high-PT
% reach_var_persub_2 = nan(size(view_time,2), 2); %subject x low-PT or high-PT
% 
% pt_low_bound = -.4;
% pt_hgh_bound = .4;
% N_min_samples = 2;
% for i_sub = 1:size(view_time,2)
%     % select out this subject's view time and directional error for each type
%     this_view_time_0 = view_time(:, i_sub, 1) - min_pt(i_sub);
%     this_direrr_0 = dir_error(:, i_sub, 1);
%     
%     this_view_time_1 = view_time(:, i_sub, 2) - min_pt(i_sub);
%     this_direrr_1 = dir_error(:, i_sub, 2);
%     
%     this_view_time_2 = view_time(:, i_sub, 3) - min_pt(i_sub);
%     this_direrr_2 = dir_error(:, i_sub, 3);
%     
%  % compute the variability (st dev) for subset of PT's in prespecified range
%     this_inds_0 = this_view_time_0 > pt_low_bound & this_view_time_0 < 0 & abs(this_direrr_0) < 30;
%     if sum(this_inds_0) > N_min_samples
%         reach_var_persub_0(i_sub, 1) = sqrt(nanvar(this_direrr_0(this_inds_0)));
%     end
%     this_inds_0 = this_view_time_0 > 0 & this_view_time_0 < pt_hgh_bound;
%     if sum(this_inds_0) > N_min_samples
%         reach_var_persub_0(i_sub, 2) = sqrt(nanvar(this_direrr_0(this_inds_0)));
%     end
%     
%     this_inds_1 = this_view_time_1 > pt_low_bound & this_view_time_1 < 0 & abs(this_direrr_1) < 30;
%     if sum(this_inds_1) > N_min_samples
%         reach_var_persub_1(i_sub, 1) = sqrt(nanvar(this_direrr_1(this_inds_1)));
%     end
%     this_inds_1 = this_view_time_1 > 0 & this_view_time_1 < pt_hgh_bound;
%     if sum(this_inds_1) > N_min_samples
%         reach_var_persub_1(i_sub, 2) = sqrt(nanvar(this_direrr_1(this_inds_1)));
%     end
%     
%     this_inds_2 = this_view_time_2 > pt_low_bound & this_view_time_2 < 0 & abs(this_direrr_2) < 30;
%     if sum(this_inds_2) > N_min_samples
%         reach_var_persub_2(i_sub, 1) = sqrt(nanvar(this_direrr_2(this_inds_2)));
%     end
%     this_inds_2 = this_view_time_2 > 0 & this_view_time_2 < pt_hgh_bound;
%     if sum(this_inds_2) > N_min_samples
%         reach_var_persub_2(i_sub, 2) = sqrt(nanvar(this_direrr_2(this_inds_2)));
%     end
% end

reach_var_persub_0 = reach_var_persub_0(sum(isnan(reach_var_persub_0),2) < 1, :);
reach_var_persub_1 = reach_var_persub_1(sum(isnan(reach_var_persub_1),2) < 1, :);
reach_var_persub_2 = reach_var_persub_2(sum(isnan(reach_var_persub_2),2) < 1, :);

figure; hold on;
errorbar([1 2], nanmedian(reach_var_persub_0, 1),...
    sqrt(nanvar(reach_var_persub_0, [], 1)./sum(~isnan(reach_var_persub_0),1)), 'r.');
errorbar([1 2], nanmedian(reach_var_persub_1, 1),...
    sqrt(nanvar(reach_var_persub_1, [], 1)./sum(~isnan(reach_var_persub_1),1)), 'g.');
errorbar([1 2], nanmedian(reach_var_persub_2, 1),...
    sqrt(nanvar(reach_var_persub_2, [], 1)./sum(~isnan(reach_var_persub_2),1)), 'b.');
axis([0 3 0 30])
%% plot probability of recall as function of appearance of symbol:
[rec_lpt, rec_hpt] = compute_recall_probability_appearance_order(...
    move_all, target_all, type_all,...
    pt_all, min_pt, 2);

max_lpt = 15;

lpt_per_targ = reshape(nanmean(rec_lpt, 3), 4, size(rec_lpt, 2));
lpt_all = reshape(nanmean(rec_lpt, 1), size(rec_lpt,2), size(rec_lpt, 3));

figure; hold on;
    errorbar(1:max_lpt, nanmean(lpt_all(1:max_lpt, :),2),...
        sqrt(nanvar(lpt_all(1:max_lpt, :), [], 2)./size(lpt_all, 2)), 'b.-',...
        'MarkerSize', 18, 'LineWidth', 3);
axis([0 max_lpt 0 1]);

[rec_lpt_direct, rec_hpt_direct] = compute_recall_probability_appearance_order(...
    move_all, target_all, type_all,...
    pt_all, min_pt, 1);

lpt_per_targ_direct = reshape(nanmean(rec_lpt_direct, 3), 4, size(rec_lpt_direct, 2));
lpt_all_direct = reshape(nanmean(rec_lpt_direct, 1), size(rec_lpt_direct,2), size(rec_lpt_direct, 3));

errorbar(1:max_lpt, nanmean(lpt_all_direct(1:max_lpt, :),2),...
        sqrt(nanvar(lpt_all_direct(1:max_lpt, :), [], 2)./size(lpt_all_direct, 2)), 'g.-',...
        'MarkerSize', 18, 'LineWidth', 3);
axis([0 max_lpt 0 1]);

