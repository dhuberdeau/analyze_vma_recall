function [data_indiv, varargout] = analyze_statLearn_individual_v1(Data, N_TARGS, varargin)


% get global variables
global n_bins EARLIEST_VALID_PT LATEST_VALID_PT SUCCESS_TH_ANGLE
EARLIEST_VALID_PT = -.2;
LATEST_VALID_PT = .85;
SUCCESS_TH_ANGLE = 30;

%% supress unneccessary warnings:
warn_id = 'MATLAB:interp1:NaNstrip';
warning('off', warn_id);

%% check if Data output included the mmperpix field (some didn't by mistake)
if isfield(Data, 'mmperpix')
    mm_pix = Data.mmperpix;
%     screen_dims = [1920, 1080];
else
    load('mm_per_pix.mat')
%     screen_dims = [1680, 1050]; % old monitor
%     screen_dims = [1600 900]; % correct monitor (12/3/2018)
end

%% decide whether to output plots
if nargin > 2
    plot_out = varargin{1};
else
    plot_out = 1;
end

%% compute view time
Data.ViewTime = Data.pPT + Data.RT; % correct error in retention_TR_experiment_v3

%% define target locations and other parameters:
% screen_dims = [1680, 1050];
screen_dims = [1920, 1080];
home_position = screen_dims/2*mm_pix;

if N_TARGS == 3
    TARG_LEN = 350*mm_pix;
    targ_angles = (0:120:300)+15;
    targ_angles(targ_angles > 180) = targ_angles(targ_angles > 180) - 360;
% targ_coords_base = TARG_LEN*[cosd(targ_angles)', sind(targ_angles)'] + home_position;
elseif N_TARGS == 4
    TARG_LEN = 350*mm_pix;
    targ_angles = 0:90:300;
    targ_angles(targ_angles > 180) = targ_angles(targ_angles > 180) - 360;
elseif N_TARGS == 6
    TARG_LEN = 300*mm_pix;
    targ_angles = (0:60:300)+15;
    targ_angles(targ_angles > 180) = targ_angles(targ_angles > 180) - 360;
else
    error('Invalid number of targets specified.')
end

% SUCCESS_TH_SD = .95;

inds_temp = 1:length(Data.Type);
type0 = inds_temp(Data.Type == 0);
type1 = inds_temp(Data.Type == 1);
type2 = inds_temp(Data.Type == 2);
type3 = inds_temp(Data.Type == 3);
type4 = inds_temp(Data.Type == 4);

targ1 = inds_temp(Data.Target == 1);
targ2 = inds_temp(Data.Target == 2);
targ3 = inds_temp(Data.Target == 3);
targ4 = inds_temp(Data.Target == 4);

% targets double up in these studies. But which version matters:
temp_targ_ind = Data.Target;
if N_TARGS == 3
    temp_targ_ind(Data.Target == 4) = 1;
    temp_targ_ind(Data.Target == 5) = 2;
    temp_targ_ind(Data.Target == 6) = 3;
elseif N_TARGS == 4
    temp_targ_ind(Data.Target == 5) = 1;
    temp_targ_ind(Data.Target == 6) = 2;
    temp_targ_ind(Data.Target == 7) = 3;
    temp_targ_ind(Data.Target == 8) = 4;
    temp_targ_ind(Data.Target == 9) = 1;
    temp_targ_ind(Data.Target == 10) = 2;
    temp_targ_ind(Data.Target == 11) = 3;
    temp_targ_ind(Data.Target == 12) = 4;
elseif N_TARGS == 6
    temp_targ_ind(Data.Target == 7) = 1;
    temp_targ_ind(Data.Target == 8) = 2;
    temp_targ_ind(Data.Target == 9) = 3;
    temp_targ_ind(Data.Target == 10) = 4;
    temp_targ_ind(Data.Target == 11) = 5;
    temp_targ_ind(Data.Target == 12) = 6;
else
   error('Invalid number of targets specified.') 
end

%% define output data matricies:
Kin_x = nan(25, length(Data.Kinematics)); %structure to hold all kinematics
Kin_y = nan(25, length(Data.Kinematics)); %structure to hold all kinematics
Dir_e = nan(1, length(Data.Kinematics)); %directional error from target.
Dir_a = nan(1, length(Data.Kinematics));
Mov_class = nan(1, length(Data.Kinematics)); %target apparently chosen.

%% Trial analysis:
vel_prelim = cell(1, length(Data.Kinematics));

H = 60;
T = 1/H;
disc_time = 4*H;
veloc_TH = 10; %cm/sec
veloc_TH_1 = 20; TH_window_1 = round([-0.075, 0.15].*H);
dist_TH = TARG_LEN/2;

if plot_out
    figure; subplot(2,2,1); 
    hold on; subplot(2,2,2); 
    hold on; subplot(2,2,3); hold on;
end

error_queue = {};
for i_tr = 1:length(Data.Kinematics)
    try
        t1 = Data.Kinematics{i_tr}(:,1) - Data.Kinematics{i_tr}(1,1);
        x1 = Data.Kinematics{i_tr}(:,2)*mm_pix - home_position(1);
        y1 = Data.Kinematics{i_tr}(:,3)*mm_pix - home_position(2);

        dt = t1(end) - t1(1);
        t = linspace(0, dt, dt*H);

        x = interp1(t1, x1, t, 'pchip');
        y = interp1(t1, y1, t, 'pchip');

        dx1 = [0 diff(x)./diff(t)];
        dy1 = [0 diff(y)./diff(t)];

        dx = sgolayfilt(dx1, 3, 5);
        dy = sgolayfilt(dy1, 3, 5);

        v = sqrt(dx.^2 + dy.^2);
    
        if plot_out
            subplot(2,2,1);
            plot(t, v)
        end

        % isolate movement...
        % discard early time during retention period
        t1 = t(disc_time:end);
        v1 = v(disc_time:end);
        x1 = x(disc_time:end);
        y1 = y(disc_time:end);
        dx1 = dx(disc_time:end);
        dy1 = dy(disc_time:end);
        
        vel_prelim{i_tr} = [t1(:), v1(:)];

        % forward-search for mvmt start to limit search to init. movement
        k0 = 1;
        v0 = v1(k0);
        while v0 < veloc_TH_1 && (k0 + 1) < length(v1)
            k0 = k0 + 1;
            v0 = v1(k0);
        end

        t2 = t1(k0+(TH_window_1(1):TH_window_1(2)));
        v2 = v1(k0+(TH_window_1(1):TH_window_1(2)));
        x2 = x1(k0+(TH_window_1(1):TH_window_1(2)));
        y2 = y1(k0+(TH_window_1(1):TH_window_1(2)));
        dx2 = dx1(k0+(TH_window_1(1):TH_window_1(2)));
        dy2 = dy1(k0+(TH_window_1(1):TH_window_1(2)));

        [v0, k0] = max(v2);
        k_max = k0;
        v_max = v0;
        while v0 > veloc_TH && k0 >1
            k0 = k0 - 1;
            v0 = v2(k0);
        end
        vf = v_max; kf = k_max;
        while vf > veloc_TH && kf < length(v2)
            kf = kf + 1;
            if kf < length(v2)
                vf = v2(kf);
            else
                break
            end
        end

        t_sub = t2(k0:kf);
        v_sub = v2(k0:kf);
        x_sub = x2(k0:kf);
        y_sub = y2(k0:kf);
        dx_sub = dx2(k0:kf);
        dy_sub = dy2(k0:kf);

        if plot_out
            subplot(2,2,2);
            plot(t_sub - t_sub(1), v_sub);
            subplot(2,2,3); hold on;
            plot(x_sub, y_sub);
            plot(x_sub(1), y_sub(1), 'g.')
        end

        Kin_x(1:min([length(x_sub), size(Kin_x,1)]), i_tr) = x_sub(1:min([length(x_sub), size(Kin_x,1)]));
        Kin_y(1:min([length(y_sub), size(Kin_y,1)]), i_tr) = y_sub(1:min([length(y_sub), size(Kin_y,1)]));

        dist = sqrt(x_sub.^2 + y_sub.^2);
        [~, k_th] = min((dist - dist_TH).^2);

        dir_abs = rad2deg(atan2(dy_sub(k_th), dx_sub(k_th)));
        targ_dir = targ_angles(temp_targ_ind(i_tr));

        [targ_off_mag, targ_class] = min(abs(dir_abs - targ_angles));
        if abs(targ_off_mag) <= SUCCESS_TH_ANGLE
            Mov_class(i_tr) = targ_class;
        end
        
        Dir_a(i_tr) = dir_abs;

        dir_err = targ_dir - dir_abs;
        dir_err(dir_err < -180) = dir_err(dir_err < -180) + 360;
        dir_err(dir_err > 180) = dir_err(dir_err > 180) - 360;
        Dir_e(i_tr) = dir_err;

    catch err_
        warning(['Trial failed. Omitting trial ', num2str(i_tr)]);
        error_queue{length(error_queue) + 1} = err_;
    end
end

%%
% Plot all kinematics to check for basic acceptability of movements
if plot_out
    figure; hold on;
    plot(Kin_x(:, Data.Type == 0), -Kin_y(:, Data.Type == 0), '-', 'Color', [172, 59, 59]/255, 'LineWidth', 2)
    plot(Kin_x(:, Data.Type == 1), -Kin_y(:, Data.Type == 1), '-', 'Color', [85, 170, 85]/255, 'LineWidth', 1)
    plot(Kin_x(:, Data.Type == 2), -Kin_y(:, Data.Type == 2), '-', 'Color', [86/255 85/255 149/255], 'LineWidth', .5)
end
%%
% View directional error and compute minimum PT
if plot_out
    figure; hold on;
    plot(Data.ViewTime(type0), Dir_e(type0), '.', 'MarkerSize', 20, 'Color', [172, 59, 59]/255)
    plot(Data.ViewTime(type1), Dir_e(type1), '.', 'MarkerSize', 20, 'Color', [85, 170, 85]/255)
    plot(Data.ViewTime(type2), Dir_e(type2), '.', 'MarkerSize', 20, 'Color', [86/255 85/255 149/255])
    
    plot(Data.ViewTime(type3), Dir_e(type3), 'o', 'MarkerSize', 8, 'LineWidth', 2, 'Color', [225, 244, 162]/255)
    plot(Data.ViewTime(type4), Dir_e(type4), 'o', 'MarkerSize', 8, 'LineWidth', 2, 'Color', [71, 113, 134]/255)
    
    axis([EARLIEST_VALID_PT, LATEST_VALID_PT, -200 200])
    set(gca,'fontsize',20)
    legend('None', 'Direct', 'Symbolic')
end

[~, pt_s, ~, pt_m] = compute_min_pt(Data.ViewTime, Data.Type, Dir_e);
% 
valid_vt = Data.ViewTime > EARLIEST_VALID_PT & Data.ViewTime < LATEST_VALID_PT;
vt = Data.ViewTime(valid_vt);
pt0 = Data.ViewTime(valid_vt & Data.Type == 0);
pt1 = Data.ViewTime(valid_vt & Data.Type == 1);
pt2 = Data.ViewTime(valid_vt & Data.Type == 2);
de0 = Dir_e(valid_vt & Data.Type == 0);
de1 = Dir_e(valid_vt & Data.Type == 1);
de2 = Dir_e(valid_vt & Data.Type == 2);

%% Compute z-scores of directional error for each condition in PT-bins
m1 = nanmean(de1);
s1 = nanstd(de1);
z0 = (de0 - m1)./s1;
z2 = (de2 - m1)./s1;
z1 = (de1 - m1)./s1;

n_bins = 6;
hist_bins = linspace(EARLIEST_VALID_PT, LATEST_VALID_PT, n_bins+1);
[n_pt_all, edge_pt_all] = histcounts(vt, hist_bins);

z_bin0 = nan(length(n_pt_all), 1);
z_bin0_sd = z_bin0;
p_bin0 = z_bin0;
psd_bin0 = z_bin0;
for i_edge = 2:length(edge_pt_all)
    this_z = pt0 > edge_pt_all(i_edge - 1) & pt0 <= edge_pt_all(i_edge);
    z_bin0(i_edge-1) = nanmean(abs(z0(this_z)));
    z_bin0_sd(i_edge - 1) = nanstd(abs(z0(this_z)));
    p_bin0(i_edge - 1) = nanmean(abs(de0(this_z)) < SUCCESS_TH_ANGLE);
    psd_bin0(i_edge - 1) = nanstd(abs(de0(this_z)) < SUCCESS_TH_ANGLE);
end

p_bin1 = nan(length(n_pt_all), 1);
psd_bin1 = p_bin1;
for i_edge = 2:length(edge_pt_all)
    this_z = pt1 > edge_pt_all(i_edge - 1) & pt1 <= edge_pt_all(i_edge);
    p_bin1(i_edge - 1) = nanmean(abs(de1(this_z)) < SUCCESS_TH_ANGLE);
    psd_bin1(i_edge - 1) = nanstd(abs(de1(this_z)) < SUCCESS_TH_ANGLE);
end

z_bin2 = nan(length(n_pt_all), 1);
z_bin2_sd = z_bin2;
p_bin2 = z_bin2;
psd_bin2 = z_bin2;
for i_edge = 2:length(edge_pt_all)
    this_z = pt2 > edge_pt_all(i_edge - 1) & pt2 <= edge_pt_all(i_edge);
    z_bin2(i_edge-1) = nanmean(abs(z2(this_z)));
    z_bin2_sd(i_edge - 1) = nanstd(abs(z2(this_z)));
    p_bin2(i_edge - 1) = nanmean(abs(de2(this_z)) < SUCCESS_TH_ANGLE);
    psd_bin2(i_edge - 1) = nanstd(abs(de2(this_z)) < SUCCESS_TH_ANGLE);
end
%% Key Graph: Directional error vs. PT
x_ind = edge_pt_all(2:end) - diff(edge_pt_all)/2;

if plot_out
    figure; hold on;
    errorbar(x_ind + 0.02, z_bin0, z_bin0_sd, 'r.-', 'MarkerSize', 18);
    errorbar(x_ind - 0.02, z_bin2, z_bin2_sd, 'b.-', 'MarkerSize', 18);

    figure; hold on;
    plot(x_ind + .02, p_bin0, '.-', 'MarkerSize', 18, 'Color', [172, 59, 59]/255);
    plot(x_ind - .02, p_bin1, '.-', 'MarkerSize', 18, 'Color', [85, 170, 85]/255);
    plot(x_ind, p_bin2, '.-', 'MarkerSize', 18, 'Color', [86/255 85/255 149/255]);
    axis([EARLIEST_VALID_PT, LATEST_VALID_PT, 0 1])
    set(gca,'fontsize',20)
    legend('None', 'Direct', 'Symbolic')
end
%% plot pPT vs. ViewTime

[pt_sort, i_sort] = sort(Data.pPT);

if plot_out
    figure; 
    subplot(2,1,1);hold on
    plot(pt_sort, Data.ViewTime(i_sort), '.', 'MarkerSize', 10);
    plot(pt_sort, pt_sort);
    axis([-.1 .7 -1.5 1.5])
    set(gca,'fontsize',20)
    subplot(2,1,2); hold on
    pt_diff = pt_sort - Data.ViewTime(i_sort);
    plot(pt_sort, pt_diff, '.', 'MarkerSize', 20, 'Color', [.5 .5 .5])
    plot(pt_sort(Data.Type(i_sort) == 0), pt_diff(Data.Type(i_sort) == 0), '.', 'MarkerSize', 20, 'Color', [172, 59, 59]/255)
    plot(pt_sort(Data.Type(i_sort) == 1), pt_diff(Data.Type(i_sort) == 1), '.', 'MarkerSize', 20, 'Color', [85, 170, 85]/255)
    plot(pt_sort(Data.Type(i_sort) == 2), pt_diff(Data.Type(i_sort) == 2), '.', 'MarkerSize', 20, 'Color', [86/255 85/255 149/255])
    axis([-.1 .7 -1 1])
    set(gca,'fontsize',20)

    figure; 
    subplot(2,1,1);hold on
    plot(pt_sort, Data.ViewTime(i_sort), '.', 'MarkerSize', 10);
    plot(pt_sort, pt_sort);
    axis([-.1 .7 -1.5 1.5])
    set(gca,'fontsize',20)
    subplot(2,1,2); hold on
    pt_diff = pt_sort - Data.ViewTime(i_sort);
    plot(pt_sort, pt_diff, '.', 'MarkerSize', 20, 'Color', [.5 .5 .5])
    plot(pt_sort(Data.Type(i_sort) == 0), pt_diff(Data.Type(i_sort) == 0), '.', 'MarkerSize', 20, 'Color', [.5 .5 .5])
    plot(pt_sort(Data.Type(i_sort) == 1), pt_diff(Data.Type(i_sort) == 1), '.', 'MarkerSize', 20, 'Color', [.5 .5 .5])
    axis([-.1 .7 -1 1])
    set(gca,'fontsize',20)
end

hist_bins = linspace(-.005, .605, n_bins+1);
[n_ppt_all, edge_ppt_all] = histcounts(pt_sort, hist_bins);

ppt0 = Data.pPT(Data.Type == 0);
ppt1 = Data.pPT(Data.Type == 1);
ppt2 = Data.pPT(Data.Type == 2);

vt_resid = Data.pPT - Data.ViewTime;

pt_z0 = (vt_resid(Data.Type == 0) - nanmean(vt_resid(Data.Type == 1)))./nanstd(vt_resid(Data.Type == 1));
pt_z2 = (vt_resid(Data.Type == 2) - nanmean(vt_resid(Data.Type == 1)))./nanstd(vt_resid(Data.Type == 1));
pt_z_bin0 = nan(length(n_ppt_all), 1);
pt_z_bin0_sd = pt_z_bin0;
pt_z_bin2 = nan(length(n_ppt_all), 1);
pt_z_bin2_sd = pt_z_bin2;
for i_edge = 2:length(edge_ppt_all)
    this_z = ppt0 >= edge_ppt_all(i_edge - 1) & ppt0 < edge_ppt_all(i_edge);
    pt_z_bin0(i_edge-1) = nanmean(abs(pt_z0(this_z)));
    pt_z_bin0_sd(i_edge - 1) = nanstd(abs(pt_z0(this_z)));
    
    this_z = ppt2 >= edge_ppt_all(i_edge - 1) & ppt2 < edge_ppt_all(i_edge);
    pt_z_bin2(i_edge-1) = nanmean(abs(pt_z2(this_z)));
    pt_z_bin2_sd(i_edge - 1) = nanstd(abs(pt_z2(this_z)));
end

x_ind = edge_ppt_all(2:end) - diff(edge_ppt_all)/2;

if plot_out
    figure; hold on
    errorbar(x_ind+.01, pt_z_bin0, pt_z_bin0_sd, '.', 'MarkerSize', 20, 'Color', [172, 59, 59]/255);
    errorbar(x_ind-.01, pt_z_bin2, pt_z_bin2_sd, '.', 'MarkerSize', 20, 'Color', [86/255 85/255 149/255]);
    axis([-.1 .7 -.5 3])
    set(gca,'fontsize',20)
end
%% 
ind = 1:length(Dir_e);
if plot_out
    figure; hold on;
    plot(ind, Data.pPT, '-');
    plot(ind, Data.ViewTime, '.-', 'LineWidth', 2, 'MarkerSize', 12);
    plot(ind(Data.Type == 0), Data.ViewTime(Data.Type == 0), 'o', 'MarkerSize', 8);
    plot(ind(Data.Type == 0 & abs(Dir_e') < SUCCESS_TH_ANGLE), Data.ViewTime(Data.Type == 0 & abs(Dir_e') < SUCCESS_TH_ANGLE), '.', 'MarkerSize', 20);

    % plot([ind(1), ind(end)], [pt_s(pt_m), pt_s(pt_m)], 'k-', 'LineWidth', 2)
    % title('Trial Type 0: No Pre-Cue');
    % a = gca;
    % set(a, 'FontSize', 20);

    %
    % ind = 1:length(Dir_e);
    % figure; hold on;
    % plot(ind, Data.pPT, '-');
    % plot(ind, Data.ViewTime, '.-',  'LineWidth', 2, 'MarkerSize', 12);
    plot(ind(Data.Type == 2), Data.ViewTime(Data.Type == 2), 'o', 'MarkerSize', 8);
    plot(ind(Data.Type == 2 & abs(Dir_e') < SUCCESS_TH_ANGLE), Data.ViewTime(Data.Type == 2 & abs(Dir_e') < SUCCESS_TH_ANGLE), '.', 'MarkerSize', 20);
    plot([ind(1), ind(end)], [pt_s(pt_m), pt_s(pt_m)], 'k-', 'LineWidth', 2)
    % title('Trial Type 2: Symbolic Pre-Cue');
    a = gca;
    set(a, 'FontSize', 20);
end

%% find proportion successful symbolic and direct trials for each block above and below minPT:
bin_edges = [73, 133, 193, 253, length(Dir_e)+1];%trials separating blocks

pr_corr_block_1 = nan(1, length(bin_edges) - 1); %direct pre-cue
pr_corr_block_2 = nan(1, length(bin_edges) - 1); %symbolic pre-cue
pr_corr_block_0 = nan(1, length(bin_edges) - 1); %no pre-cue
pr_corr_supTH_block_1 = nan(1, length(bin_edges) - 1); %direct pre-cue (above minPT)
pr_corr_supTH_block_2 = nan(1, length(bin_edges) - 1); %symbolic pre-cue (above minPT)
pr_corr_supTH_block_0 = nan(1, length(bin_edges) - 1); %no pre-cue (above minPT)

pr_corr_block_1_target = nan(3, length(bin_edges) - 1); %direct/above minPT/all targets
pr_corr_block_2_target = nan(3, length(bin_edges) - 1); %direct/above minPT/all targets
pr_corr_supTH_block_1_target = nan(3, length(bin_edges) - 1); %direct pre-cue (above minPT)
pr_corr_supTH_block_2_target = nan(3, length(bin_edges) - 1); %symbolic pre-cue (above minPT)

if N_TARGS == 3
    target_list = [1 2 3; 4 5 6]; %symbols 1=4, 2=5, 3=6 for targets.
elseif N_TARGS == 4
    target_list = [1 2 3 4; 5 6 7 8; 9 10 11 12];
elseif N_TARGS == 6
    target_list = [1 2 3 4 5 6; 7 8 9 10 11 12];
else
    error('Invalid number of targets specified');
end
for i_bin = 1:(length(bin_edges) - 1)
    bin_inds = bin_edges(i_bin):(bin_edges(i_bin + 1) - 1);
    this_de = Dir_e(bin_inds);
    this_type = Data.Type(bin_inds);
    this_pt = Data.ViewTime(bin_inds);
    pr_corr_block_0(i_bin) = nanmean(abs(this_de(this_type == 0 & this_pt < pt_s(pt_m))) < SUCCESS_TH_ANGLE);
    pr_corr_block_1(i_bin) = nanmean(abs(this_de(this_type == 1 & this_pt < pt_s(pt_m))) < SUCCESS_TH_ANGLE);
    pr_corr_block_2(i_bin) = nanmean(abs(this_de(this_type == 2 & this_pt < pt_s(pt_m))) < SUCCESS_TH_ANGLE);
    
    pr_corr_supTH_block_0(i_bin) = nanmean(abs(this_de(this_type == 0 & this_pt > pt_s(pt_m))) < SUCCESS_TH_ANGLE);
    pr_corr_supTH_block_1(i_bin) = nanmean(abs(this_de(this_type == 1 & this_pt > pt_s(pt_m))) < SUCCESS_TH_ANGLE);
    pr_corr_supTH_block_2(i_bin) = nanmean(abs(this_de(this_type == 2 & this_pt > pt_s(pt_m))) < SUCCESS_TH_ANGLE);
    
    this_target = Data.Target(bin_inds);
    % do the same but for each target location separately.
    for i_targ = 1:size(target_list,2)
        cond1 = this_type == 1;
        cond2 = this_pt < pt_s(pt_m);
        if N_TARGS == 3
            cond3 = (this_target == target_list(1,i_targ) |...
            this_target == target_list(2,i_targ));
        elseif N_TARGS == 4
            cond3 = (this_target == target_list(1,i_targ) |...
            this_target == target_list(2,i_targ) |...
            this_target == target_list(3,i_targ));
        elseif N_TARGS == 6
            cond3 = (this_target == target_list(1,i_targ) |...
            this_target == target_list(2,i_targ));
        end
        pr_corr_block_1_target(i_targ, i_bin) = nanmean(abs(this_de(...
            cond1 & cond2 & cond3)) < SUCCESS_TH_ANGLE);
        
        cond1 = this_type == 2;
        pr_corr_block_2_target(i_targ, i_bin) = nanmean(abs(this_de(...
            cond1 & cond2 & cond3)) < SUCCESS_TH_ANGLE);
    end
    for i_targ = 1:size(target_list,2)
        cond1 = this_type == 1;
        cond2 = this_pt >= pt_s(pt_m);
        if N_TARGS == 3
            cond3 = (this_target == target_list(1,i_targ) |...
            this_target == target_list(2,i_targ));
        elseif N_TARGS == 4
            cond3 = (this_target == target_list(1,i_targ) |...
            this_target == target_list(2,i_targ) |...
            this_target == target_list(3,i_targ));
        elseif N_TARGS == 6
            cond3 = (this_target == target_list(1,i_targ) |...
            this_target == target_list(2,i_targ));
        end
        pr_corr_supTH_block_1_target(i_targ, i_bin) = nanmean(abs(this_de(...
            cond1 & cond2 & cond3)) < SUCCESS_TH_ANGLE);
        
        cond1 = this_type == 2;
        pr_corr_supTH_block_2_target(i_targ, i_bin) = nanmean(abs(this_de(...
            cond1 & cond2 & cond3)) < SUCCESS_TH_ANGLE);
    end
end

if plot_out
%     figure; hold on;
%     plot(pr_corr_supTH_block_1, 'gs', 'MarkerSize', 12, 'LineWidth', 3);
%     plot(pr_corr_supTH_block_2, 'bs', 'LineWidth', 3);
%     axis([.5 4.5 0 1]);
%     title('Above min-PT')
% 
%     figure; hold on;
%     plot(pr_corr_block_1, 'go', 'MarkerSize', 12, 'LineWidth', 3);
%     plot(pr_corr_block_2, 'bo', 'LineWidth', 3);
%     axis([.5 4.5 0 1]);
%     title('Below min-PT')
end

%% find proportion successful symbolic and direct trials for each block above and below minPT:
bin_edges = [73, 193, length(Dir_e)+1];%trials separating blocks

pr_corr_block_1_symbol = nan(6, length(bin_edges) - 1); %direct/above minPT/all targets
pr_corr_block_2_symbol = nan(6, length(bin_edges) - 1); %direct/above minPT/all targets
pr_corr_supTH_block_1_symbol = nan(6, length(bin_edges) - 1); %direct pre-cue (above minPT)
pr_corr_supTH_block_2_symbol = nan(6, length(bin_edges) - 1); %symbolic pre-cue (above minPT)

if N_TARGS == 3
    symbol_list = [1 2 3 4 5 6]; %symbols 1=4, 2=5, 3=6 for targets.
elseif N_TARGS == 4
    symbol_list = [1 2 3 4 5 6 7 8 9 10 11 12];
elseif N_TARGS == 6
    symbol_list = [1 2 3 4 5 6 7 8 9 10 11 12];
else
    error('Invalid number of targets specified');
end
for i_bin = 1:(length(bin_edges) - 1)
    bin_inds = bin_edges(i_bin):(bin_edges(i_bin + 1) - 1);
    this_de = Dir_e(bin_inds);
    this_type = Data.Type(bin_inds);
    this_pt = Data.ViewTime(bin_inds);
    
    this_target = Data.Target(bin_inds);
    
    % do the same but for each symbol separately.
    for i_sym = 1:length(symbol_list)
        pr_corr_block_1_symbol(i_sym, i_bin) = nanmean(abs(this_de(...
            this_type == 1 &...
            this_pt < pt_s(pt_m) &...
            this_target == symbol_list(i_sym))) < SUCCESS_TH_ANGLE);
        
        pr_corr_block_2_symbol(i_sym, i_bin) = nanmean(abs(this_de(...
            this_type == 2 &...
            this_pt < pt_s(pt_m) &...
            this_target == symbol_list(i_sym))) < SUCCESS_TH_ANGLE);
    end
    for i_sym = 1:length(symbol_list)
        pr_corr_supTH_block_1_symbol(i_sym, i_bin) = nanmean(abs(this_de(...
            this_type == 1 &...
            this_pt >= pt_s(pt_m) &...
            this_target == symbol_list(i_sym))) < SUCCESS_TH_ANGLE);
        
        pr_corr_supTH_block_2_symbol(i_sym, i_bin) = nanmean(abs(this_de(...
            this_type == 2 &...
            this_pt >= pt_s(pt_m) &...
            this_target == symbol_list(i_sym))) < SUCCESS_TH_ANGLE);
    end
end

%% find proportion successful symbolic and direct trials above and below minPT: 
% marginally across experiment
bin_edges = [73, length(Dir_e)+1];%trials separating blocks

pr_corr_block_1_symbol_marginal = nan(6, length(bin_edges) - 1); %direct/above minPT/all targets
pr_corr_block_2_symbol_marginal = nan(6, length(bin_edges) - 1); %direct/above minPT/all targets
pr_corr_supTH_block_1_symbol_marginal = nan(6, length(bin_edges) - 1); %direct pre-cue (above minPT)
pr_corr_supTH_block_2_symbol_marginal = nan(6, length(bin_edges) - 1); %symbolic pre-cue (above minPT)

if N_TARGS == 3
    symbol_list = [1 2 3 4 5 6]; %symbols 1=4, 2=5, 3=6 for targets.
elseif N_TARGS == 4
    symbol_list = [1 2 3 4 5 6 7 8 9 10 11 12];
elseif N_TARGS == 6
    symbol_list = [1 2 3 4 5 6 7 8 9 10 11 12];
else
    error('Invalid number of targets specified');
end
for i_bin = 1:(length(bin_edges) - 1)
    bin_inds = bin_edges(i_bin):(bin_edges(i_bin + 1) - 1);
    this_de = Dir_e(bin_inds);
    this_type = Data.Type(bin_inds);
    this_pt = Data.ViewTime(bin_inds);
    
    this_target = Data.Target(bin_inds);
    
    % do the same but for each symbol separately.
    for i_sym = 1:length(symbol_list)
        pr_corr_block_1_symbol_marginal(i_sym, i_bin) = nanmean(abs(this_de(...
            this_type == 1 &...
            this_pt < pt_s(pt_m) &...
            this_target == symbol_list(i_sym))) < SUCCESS_TH_ANGLE);
        
        pr_corr_block_2_symbol_marginal(i_sym, i_bin) = nanmean(abs(this_de(...
            this_type == 2 &...
            this_pt < pt_s(pt_m) &...
            this_target == symbol_list(i_sym))) < SUCCESS_TH_ANGLE);
    end
    for i_sym = 1:length(symbol_list)
        pr_corr_supTH_block_1_symbol_marginal(i_sym, i_bin) = nanmean(abs(this_de(...
            this_type == 1 &...
            this_pt >= pt_s(pt_m) &...
            this_target == symbol_list(i_sym))) < SUCCESS_TH_ANGLE);
        
        pr_corr_supTH_block_2_symbol_marginal(i_sym, i_bin) = nanmean(abs(this_de(...
            this_type == 2 &...
            this_pt >= pt_s(pt_m) &...
            this_target == symbol_list(i_sym))) < SUCCESS_TH_ANGLE);
    end
end

%% find distribution among symbols of the order of first appearance in exp.


%% plot subject's target confusion matrix:
conf_mat = confusionmat(temp_targ_ind', Mov_class);


%% plot target individual learning curves

target_symbols = {'s', 'o', 'x', '*', '+', 'd'};

if plot_out
    figure; hold on;
    for i_targ = 1:N_TARGS
        plot(1:4, pr_corr_block_1_target(i_targ,:), ...
            ['g', target_symbols{i_targ}, '-'],...
            'MarkerSize', 10,  'LineWidth', 3);
        
        plot(1:4, pr_corr_block_2_target(i_targ,:),...
            ['b', target_symbols{i_targ}, '-'],...
            'MarkerSize', 10,  'LineWidth', 3);
    end
    axis([0.5 4.5 0 1]);
    title('Below min-PT')

    figure; hold on;
    for i_targ = 1:N_TARGS
        plot(1:4, pr_corr_supTH_block_1_target(i_targ,:), ...
            ['g', target_symbols{i_targ}, '-'],...
            'MarkerSize', 10,  'LineWidth', 3);
        
        plot(1:4, pr_corr_supTH_block_2_target(i_targ,:),...
            ['b', target_symbols{i_targ}, '-'],...
            'MarkerSize', 10,  'LineWidth', 3);
    end
    axis([0.5 4.5 0 1]);
    title('Above min-PT')
end

%% plot symbol individual learning curves
symbol_symbols = {'s', 'o', 'x', '*', '+', 'd', '.', '-', '^', '<', '>', 'p'};
switch N_TARGS
    case 3
        n_symbs = 6;
    case 4
        n_symbs = 12;
    case 6
        n_symbs = 12;
    otherwise
        errro('Invalid number of targets specified');
end
            
if plot_out
    figure; hold on;
    for i_sym = 1:n_symbs
        plot([1.5, 3.5], pr_corr_block_1_symbol(i_sym,:),...
            ['g', symbol_symbols{i_sym}, '-'],...
            'MarkerSize', 10, 'LineWidth', 1.5);
    end
    
    for i_sym = 1:n_symbs
        plot([1.5, 3.5], pr_corr_block_2_symbol(i_sym,:),...
            ['b', symbol_symbols{i_sym}, '-'],...
            'MarkerSize', 10, 'LineWidth', 1.5);
    end
    axis([0.5 4.5 0 1]);
    title('Below min-PT')

    figure; hold on;
    for i_sym = 1:n_symbs
        plot([1.5, 3.5], pr_corr_supTH_block_1_symbol(i_sym,:),...
            ['g', symbol_symbols{i_sym}, '-'],...
            'MarkerSize', 10, 'LineWidth', 1.5);
    end
    
    for i_sym = 1:n_symbs
        plot([1.5, 3.5], pr_corr_supTH_block_2_symbol(i_sym,:),...
            ['b', symbol_symbols{i_sym}, '-'],...
            'MarkerSize', 10, 'LineWidth', 1.5);
    end

    axis([0.5 4.5 0 1]);
    title('Above min-PT')
end
%% find proportion successful symbolic and direct trials for each block above and below minPT,
% with higher resolution (half-blocks) than above:
bin_edges = 73:30:313;%trials
pr_corr_block_1_res = nan(1, length(bin_edges) - 1); %direct pre-cue
pr_corr_block_2_res = nan(1, length(bin_edges) - 1); %symbolic pre-cue
pr_corr_block_0_res = nan(1, length(bin_edges) - 1); %no pre-cue
pr_corr_supTH_block_1_res = nan(1, length(bin_edges) - 1); %direct pre-cue (above minPT)
pr_corr_supTH_block_2_res = nan(1, length(bin_edges) - 1); %symbolic pre-cue (above minPT)
pr_corr_supTH_block_0_res = nan(1, length(bin_edges) - 1); % no pre-cue
for i_bin = 1:(length(bin_edges) - 1)
    bin_inds = bin_edges(i_bin):(bin_edges(i_bin + 1) - 1);
    this_de = Dir_e(bin_inds);
    this_type = Data.Type(bin_inds);
    this_pt = Data.ViewTime(bin_inds);
    pr_corr_block_0_res(i_bin) = nanmean(abs(this_de(this_type == 0 & this_pt < pt_s(pt_m))) < SUCCESS_TH_ANGLE);
    pr_corr_block_1_res(i_bin) = nanmean(abs(this_de(this_type == 1 & this_pt < pt_s(pt_m))) < SUCCESS_TH_ANGLE);
    pr_corr_block_2_res(i_bin) = nanmean(abs(this_de(this_type == 2 & this_pt < pt_s(pt_m))) < SUCCESS_TH_ANGLE);
    
    pr_corr_supTH_block_0_res(i_bin) = nanmean(abs(this_de(this_type == 0 & this_pt > pt_s(pt_m))) < SUCCESS_TH_ANGLE);
    pr_corr_supTH_block_1_res(i_bin) = nanmean(abs(this_de(this_type == 1 & this_pt > pt_s(pt_m))) < SUCCESS_TH_ANGLE);
    pr_corr_supTH_block_2_res(i_bin) = nanmean(abs(this_de(this_type == 2 & this_pt > pt_s(pt_m))) < SUCCESS_TH_ANGLE);
end

if plot_out
    figure; hold on;
    plot(pr_corr_supTH_block_1_res, 'gs', 'MarkerSize', 12, 'LineWidth', 3);
    plot(pr_corr_supTH_block_2_res, 'bs', 'LineWidth', 3);
    plot(pr_corr_supTH_block_0_res, 'rs', 'LineWidth', 3);
    axis([.5 8.5 0 1]);
    title('Above min-PT')
    
    figure; hold on;
    plot(pr_corr_block_1_res, 'go', 'MarkerSize', 12, 'LineWidth', 3);
    plot(pr_corr_block_2_res, 'bo', 'LineWidth', 3);
    plot(pr_corr_block_0_res, 'ro', 'LineWidth', 3);
    axis([.5 8.5 0 1]);
    title('Below min-PT')
end

%% assign output:
data_indiv.p_bins = {p_bin0, p_bin1, p_bin2};
data_indiv.z_bins = {pt_z_bin0, pt_z_bin2};
data_indiv.z_bin_sd = {pt_z_bin0_sd, pt_z_bin2_sd};
data_indiv.min_pt = pt_s(pt_m);
data_indiv.Dir_e = Dir_e;
data_indiv.Dir_a = Dir_a;
data_indiv.Mov_class = Mov_class;
data_indiv.ViewTime = Data.ViewTime;
data_indiv.types = {type0, type1, type2, type3, type4};
data_indiv.kinematics = {Kin_x, Kin_y};
data_indiv.targ_len = TARG_LEN;
data_indiv.prCorr_lowRes = {pr_corr_block_0, pr_corr_block_1,...
    pr_corr_block_2};
data_indiv.prCorr_supTH_lowRes = {pr_corr_supTH_block_0, pr_corr_supTH_block_1,...
    pr_corr_supTH_block_2};
data_indiv.prCorr_hiRes = {pr_corr_block_0_res, ...
    pr_corr_block_1_res, pr_corr_block_2_res};
data_indiv.prCorr_supTH_hiRes = {pr_corr_supTH_block_0_res, ...
    pr_corr_supTH_block_1_res, pr_corr_supTH_block_2_res};
data_indiv.prCorr_lowRes_target = {pr_corr_block_1_target,...
    pr_corr_block_2_target};
data_indiv.prCorr_supTH_lowRes_target = {pr_corr_supTH_block_1_target,...
    pr_corr_supTH_block_2_target};
data_indiv.prCorr_lowRes_symbol = {pr_corr_block_1_symbol,...
    pr_corr_block_2_symbol};
data_indiv.prCorr_supTH_lowRes_symbol = {pr_corr_supTH_block_1_symbol,...
    pr_corr_supTH_block_2_symbol};
data_indiv.prCorr_symbol_marginal = {pr_corr_block_1_symbol_marginal,...
    pr_corr_block_2_symbol_marginal};
data_indiv.prCorr_supTH_symbol_marginal = {pr_corr_supTH_block_1_symbol_marginal,...
    pr_corr_supTH_block_2_symbol_marginal};
data_indiv.confusion_matrix = conf_mat;

varargout{1} = error_queue;