function [data_indiv, varargout] = analyze_retention_individual_v1(Data, varargin)
% Analysis of individual participant's data in recall variant of the
% visuomotor association (VMA) task with timed response.
%
% Input: Data - a data structure saved by the psychtoolbox experiment.
%
% Output: data_indiv - a structure of results from analysis of this
% individual for use in plotting aggregated results.
%
% David Huberdeau, last updated: 12/4/2018

% set global variables
global EARLIEST_VALID_PT LATEST_VALID_PT SUCCESS_TH_ANGLE
% EARLIEST_VALID_PT = -.2;
% LATEST_VALID_PT = .85;
% SUCCESS_TH_ANGLE = 30;

%% supress unneccessary warnings:
warn_id = 'MATLAB:interp1:NaNstrip';
warning('off', warn_id);

%% define camera calibration:
if isfield(Data, 'mmperpix')
    mm_pix = Data.mmperpix;
    screen_dims = [1920, 1080];
else
    load('mm_per_pix.mat')
%     screen_dims = [1680, 1050];
    screen_dims = [1600 900];
end
    
%% decide whether to output plots
if nargin > 1
    plot_out = varargin{1};
else
    plot_out = 1;
end
%% define some variables and constants:
Data.ViewTime = Data.pPT + Data.RT; % compute ViewTime from RT and prescribed PT

home_position = screen_dims/2*mm_pix;
TARG_LEN = 400*mm_pix;
targ_angles = 0:90:300;
targ_angles(targ_angles > 180) = targ_angles(targ_angles > 180) - 360;
targ_coords_base = TARG_LEN*[cosd(targ_angles)', sind(targ_angles)'] + home_position;

SUCCESS_TH_SD = .95;
% SUCCESS_TH_ANGLE = 45;

% EARLIEST_VALID_PT = -.2;
% LATEST_VALID_PT = .85;

MIN_N_PTS_FOR_MEASURE = 4;

TH_err_TH = 1;
START_DIST_TH = 1;


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

Kin_x = nan(25, length(Data.Kinematics)); %matrix to hold all kinematics
Kin_y = nan(25, length(Data.Kinematics)); %matrix to hold all kinematics
kin_PV = nan(1, length(Data.Kinematics)); %matrix to hold peak velocity
kin_MT = nan(1, length(Data.Kinematics)); %matrix to hold Movement time
kin_VAR = nan(1, length(Data.Kinematics)); %matrix to hold within-mvmt var.
Dir_e = nan(1, length(Data.Kinematics)); %directional error from target.
Dir_a = nan(1, length(Data.Kinematics)); %absolute direction of initial reach.

vel_prelim = cell(1, length(Data.Kinematics));

%% analyze each trial:
H = 60;
T = 1/H;
disc_time = 4*H;
veloc_TH = 10; %cm/sec
veloc_TH_1 = 20; 
TH_window_1 = round([-0.075, 0.20].*H); %limit search to [-.075, .2] sec window
% dist_TH = TARG_LEN/4;
dist_TH = 3;
if plot_out
    figure; subplot(2,2,1); hold on; 
    subplot(2,2,2); hold on; 
    subplot(2,2,3); hold on;
end
error_queue = {};
for i_tr = 1:length(Data.Kinematics)
    try
        t1 = Data.Kinematics{i_tr}(:,1) - Data.Kinematics{i_tr}(1,1);
%         x1 = (Data.Kinematics{i_tr}(:,2) - home_position(1))*mm_pix;
%         y1 = (Data.Kinematics{i_tr}(:,3) - home_position(2))*mm_pix;

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
        
        ax1 = [0 diff(dx)./diff(t)];
        ay1 = [0 diff(dy)./diff(t)];
        
        ax = sgolayfilt(ax1, 3, 5);
        ay = sgolayfilt(ay1, 3, 5);

        v = sqrt(dx.^2 + dy.^2);
        a = sqrt(ax.^2 + ay.^2); 

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
        v1 = v(disc_time:end);
        ax1 = ax(disc_time:end);
        ay1 = ay(disc_time:end);
        a1 = a(disc_time:end);
        
        vel_prelim{i_tr} = [t1(:), v1(:)];

        % forward-search for mvmt start to limit search to init. movement
        k0 = 1;
        v0 = v1(k0);
        while v0 < veloc_TH_1 && (k0 + 1) < length(v1)
            k0 = k0 + 1;
            v0 = v1(k0);
        end

        % take a relatively small window of time around the movement
        % initiation as determined by forwards search:
        t2 = t1(k0+(TH_window_1(1):TH_window_1(2)));
        x2 = x1(k0+(TH_window_1(1):TH_window_1(2)));
        y2 = y1(k0+(TH_window_1(1):TH_window_1(2)));
        v2 = v1(k0+(TH_window_1(1):TH_window_1(2)));
        a2 = a1(k0+(TH_window_1(1):TH_window_1(2)));
        dx2 = dx1(k0+(TH_window_1(1):TH_window_1(2)));
        dy2 = dy1(k0+(TH_window_1(1):TH_window_1(2)));
        ax2 = ax1(k0+(TH_window_1(1):TH_window_1(2)));
        ay2 = ay1(k0+(TH_window_1(1):TH_window_1(2)));

        % search backward in time from the peak velocity to find the time
        % at which the velocity first exceeded a low velocity threshold.
        [v0, k0] = max(v2);
        k_max = k0;
        v_max = v0;
        while v0 > veloc_TH && k0 >1
            k0 = k0 - 1;
            v0 = v2(k0);
        end
        
        % search forward in time from the peak velocity to find the time at
        % which the velocity last exceeded the low velocity threshold to
        % find movement end 
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
        x_sub = x2(k0:kf);
        y_sub = y2(k0:kf);
        v_sub = v2(k0:kf);
        a_sub = a2(k0:kf);
        dx_sub = dx2(k0:kf);
        dy_sub = dy2(k0:kf);
        ax_sub = ax2(k0:kf);
        ay_sub = ay2(k0:kf);

        dist = sqrt(x_sub.^2 + y_sub.^2);
        [TH_err, k_th] = min(abs(dist - dist_TH));
        
        if dist(1) < START_DIST_TH
            if plot_out
                subplot(2,2,2);
                plot(t_sub - t_sub(1), v_sub);

                subplot(2,2,3); hold on;
                plot(x_sub, y_sub);
                plot(x_sub(1), y_sub(1), 'g.')
            end
            
            % record movement kinematics and kinematic measures:
            Kin_x(1:min([length(x_sub), size(Kin_x,1)]), i_tr) = x_sub(1:min([length(x_sub), size(Kin_x,1)]));
            Kin_y(1:min([length(y_sub), size(Kin_y,1)]), i_tr) = y_sub(1:min([length(y_sub), size(Kin_y,1)]));
            
            kin_PV(i_tr) = v_max;
            kin_MT(i_tr) = t_sub(end) - t_sub(1);
            
            [f_sub, p_sub] = simple_psd(a_sub, H);
            f_targ = f_sub > 4 & f_sub < 8;
            kin_VAR(i_tr) = nanmean(p_sub(f_targ));
            
            if TH_err < TH_err_TH

                dir_abs = rad2deg(atan2(dy_sub(k_th), dx_sub(k_th)));
                targ_dir = targ_angles(Data.Target(i_tr));

                Dir_a(i_tr) = dir_abs;

                dir_err = targ_dir - dir_abs;
                dir_err(dir_err < -180) = dir_err(dir_err < -180) + 360;
                dir_err(dir_err > 180) = dir_err(dir_err > 180) - 360;
                Dir_e(i_tr) = dir_err;
            else
                Dir_e(i_tr) = nan;
                Dir_a(i_tr) = nan;
            end
            
        else
            error('Trajectory measured as starting outside of Threshold region. Trial Failure.');
        end

    catch err_
        warning(['Trial failed. Omitting trial ', num2str(i_tr)]);
        error_queue{length(error_queue) + 1} = err_;
    end
end

if plot_out
    %% plot trajectories to assess raw behavior.
    % Type 0 - No pre-cue
    % Type 1 - direct pre-cue
    % Type 2 - Symbolic pre-cue
    figure; hold on;
    plot(Kin_x(:, Data.Type == 0), -Kin_y(:, Data.Type == 0), '-', 'Color', [172, 59, 59]/255, 'LineWidth', 2)
    plot(Kin_x(:, Data.Type == 1), -Kin_y(:, Data.Type == 1), '-', 'Color', [85, 170, 85]/255, 'LineWidth', 2)
    plot(Kin_x(:, Data.Type == 2), -Kin_y(:, Data.Type == 2), '-', 'Color', [86/255 85/255 149/255], 'LineWidth', 2)

    %% Plot directional error
    figure; hold on;
    % plot(Data.ViewTime([type0, type1, type2]), Dir_e([type0, type1, type2]) + 10*rand(1, length([type0, type1, type2])), '.', 'MarkerSize', 18);
    plot(Data.ViewTime(type0), Dir_e(type0), '.', 'MarkerSize', 20, 'Color', [172, 59, 59]/255)
    plot(Data.ViewTime(type1), Dir_e(type1), '.', 'MarkerSize', 20, 'Color', [85, 170, 85]/255)
    plot(Data.ViewTime(type2), Dir_e(type2), '.', 'MarkerSize', 20, 'Color', [86/255 85/255 149/255])

    plot(Data.ViewTime(type3), Dir_e(type3), 'o', 'MarkerSize', 8, 'LineWidth', 2, 'Color', [225, 244, 162]/255)
    plot(Data.ViewTime(type4), Dir_e(type4), 'o', 'MarkerSize', 8, 'LineWidth', 2, 'Color', [71, 113, 134]/255)

    plot([-.2 .8], [45 45], 'k')
    plot([-.2 .8], -[45 45], 'k')
    axis([EARLIEST_VALID_PT, LATEST_VALID_PT, -200 200])
    set(gca,'fontsize',20)
    legend('None', 'Direct', 'Symbolic')
end

%% compute minimum PT for this individual:

[min_pt, pt_s, fun, pt_m] = compute_min_pt(Data.ViewTime, Data.Type, Dir_e);

% valid_vt = Data.ViewTime > EARLIEST_VALID_PT & Data.ViewTime < LATEST_VALID_PT;
% vt = Data.ViewTime(valid_vt);
% pt0 = Data.ViewTime(valid_vt & Data.Type == 0);
% pt1 = Data.ViewTime(valid_vt & Data.Type == 1);
% pt2 = Data.ViewTime(valid_vt & Data.Type == 2);
% de0 = Dir_e(valid_vt & Data.Type == 0);
% de1 = Dir_e(valid_vt & Data.Type == 1);
% de2 = Dir_e(valid_vt & Data.Type == 2);
% 
% % pt0 = Data.ViewTime(type0);
% % de0 = Dir_e(type0);
% pe0 = de0 <30 & de0 > -30;
% 
% pt_s = linspace(min(pt0), max(pt0), 30);
% fun = nan(size(pt_s));
% for i_s = 1:length(pt_s)
% %     fun(i_s) = sum(pe0(pt0 < .4 & pt0 > pt_s(i_s)))/sum(pt0 < .4 & pt0 > pt_s(i_s)) + sum(1 - pe0(pt0 > 0.1 & pt0 <= pt_s(i_s)))/sum(pt0 > 0.1 & pt0 <= pt_s(i_s));
%     fun(i_s) = sum(pe0(pt0 > pt_s(i_s))) + sum(1 - pe0(pt0 <= pt_s(i_s)));
% end
% [fun_m, pt_m] = max(fun);
% closeness_to_maximum = fun_m*.05; % 5% of maximum
% if sum(abs(fun - fun_m) < closeness_to_maximum) > 0 %there is a point nearly equal to the max
%     % take the later point's index as the min-PT
%     temp_inds = 1:length(fun);
%     equal_max_inds = temp_inds(abs(fun - fun_m) < closeness_to_maximum);
%     pt_m = max(equal_max_inds); %take the latest index among (near) equals
% end
% 
% % Comment below out: 12/16/2018 - make consistent with superior method used
% % in statistical learning paradigm.
% %
% % min_pt_est_set = nan(1,1000);
% % for i_set = 1:length(min_pt_est_set)
% %     
% %     %select subsample of pt0:
% %     rand_select = randperm(length(pt0));
% %     pt0_ = pt0(rand_select(1:(ceil(length(pt0)*.8))));
% %     pe0_ = pe0(rand_select(1:(ceil(length(pt0)*.8))));
% %     
% %     pt_s = linspace(min(pt0_), max(pt0_), 30);
% %     fun = nan(size(pt_s));
% %     for i_s = 1:length(pt_s)
% %     %     fun(i_s) = sum(pe0(pt0 < .4 & pt0 > pt_s(i_s)))/sum(pt0 < .4 & pt0 > pt_s(i_s)) + sum(1 - pe0(pt0 > 0.1 & pt0 <= pt_s(i_s)))/sum(pt0 > 0.1 & pt0 <= pt_s(i_s));
% %         fun(i_s) = sum(pe0_(pt0_ > pt_s(i_s))) + sum(1 - pe0_(pt0_ <= pt_s(i_s)));
% %     end
% % %     plot(pt_s, fun,'k')
% % 
% %     %%%%%% MAIN + MOD: find max but look for equal or near equal max's later
% %     [fun_m, pt_m] = max(fun);
% %     if sum(abs(fun - fun_m) < 1) > 0 %there is a point nearly equal to the max
% %         % take the later point's index as the min-PT
% %         temp_inds = 1:length(fun);
% %         equal_max_inds = temp_inds(abs(fun - fun_m) < 1);
% %         pt_m = max(equal_max_inds); %take the latest index among (near) equals
% %     end
% % %     plot([pt_s(pt_m), pt_s(pt_m)], [-200 200], 'k-', 'LineWidth', 2)
% %     min_pt_est_set(i_set) = pt_s(pt_m);
% % end
% % pt_s = mean(min_pt_est_set(min_pt_est_set > 0));
% % pt_m = 1;
% 
% if plot_out 
%     plot([pt_s(pt_m), pt_s(pt_m)], [-200 200], 'k-', 'LineWidth', 2)
% end

%% Plot probability of success at PT bins
% figure; hold on;

valid_vt = Data.ViewTime > EARLIEST_VALID_PT & Data.ViewTime < LATEST_VALID_PT;
vt = Data.ViewTime(valid_vt);
pt0 = Data.ViewTime(valid_vt & Data.Type == 0);
pt1 = Data.ViewTime(valid_vt & Data.Type == 1);
pt2 = Data.ViewTime(valid_vt & Data.Type == 2);
de0 = Dir_e(valid_vt & Data.Type == 0);
de1 = Dir_e(valid_vt & Data.Type == 1);
de2 = Dir_e(valid_vt & Data.Type == 2);

m1 = nanmean(de1);
s1 = nanstd(de1);
z0 = (de0 - m1)./s1;
z2 = (de2 - m1)./s1;
z1 = (de1 - m1)./s1;

n_bins = 8;
hist_bins = linspace(EARLIEST_VALID_PT, LATEST_VALID_PT, n_bins+1);
[n_pt_all, edge_pt_all] = histcounts(vt, hist_bins);

z_bin0 = nan(length(n_pt_all), 1);
z_bin0_sd = z_bin0;
p_bin0 = z_bin0;
psd_bin0 = z_bin0;
for i_edge = 2:length(edge_pt_all)
    this_z = pt0 > edge_pt_all(i_edge - 1) & pt0 <= edge_pt_all(i_edge);
    if sum(this_z) >= MIN_N_PTS_FOR_MEASURE
        z_bin0(i_edge-1) = nanmean(abs(z0(this_z)));
        z_bin0_sd(i_edge - 1) = nanstd(abs(z0(this_z)));
        p_bin0(i_edge - 1) = sum(abs(de0(this_z)) < SUCCESS_TH_ANGLE)/sum(~isnan(de0(this_z)));
        psd_bin0(i_edge - 1) = nanstd(abs(de0(this_z)) < SUCCESS_TH_ANGLE);
    end
end

p_bin1 = nan(length(n_pt_all), 1);
psd_bin1 = p_bin1;
for i_edge = 2:length(edge_pt_all)
    this_z = pt1 > edge_pt_all(i_edge - 1) & pt1 <= edge_pt_all(i_edge);
    if sum(this_z) >= MIN_N_PTS_FOR_MEASURE
        p_bin1(i_edge - 1) = sum(abs(de1(this_z)) < SUCCESS_TH_ANGLE)/sum(~isnan(de1(this_z)));
        psd_bin1(i_edge - 1) = nanstd(abs(de1(this_z)) < SUCCESS_TH_ANGLE);
    end
end

z_bin2 = nan(length(n_pt_all), 1);
z_bin2_sd = z_bin2;
p_bin2 = z_bin2;
psd_bin2 = z_bin2;
for i_edge = 2:length(edge_pt_all)
    this_z = pt2 > edge_pt_all(i_edge - 1) & pt2 <= edge_pt_all(i_edge);
    if sum(this_z) >= MIN_N_PTS_FOR_MEASURE
        z_bin2(i_edge-1) = nanmean(abs(z2(this_z)));
        z_bin2_sd(i_edge - 1) = nanstd(abs(z2(this_z)));
        p_bin2(i_edge - 1) = sum(abs(de2(this_z)) < SUCCESS_TH_ANGLE)/sum(~isnan(de2(this_z)));
        psd_bin2(i_edge - 1) = nanstd(abs(de2(this_z)) < SUCCESS_TH_ANGLE);
    end
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

n_bins = 6; hist_bins = linspace(-.005, .605, n_bins+1);
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

if plot_out
    %% 
    ind = 1:length(Dir_e);
    figure; hold on;
    plot(ind, Data.pPT, '-');
    plot(ind, Data.ViewTime, '.-');
    plot(ind(Data.Type == 0), Data.ViewTime(Data.Type == 0), 'o');
    plot(ind(Data.Type == 0 & abs(Dir_e') < 45), Data.ViewTime(Data.Type == 0 & abs(Dir_e') < 45), 'x');
    title('Trial Type 0: No Pre-Cue');

    %% 
    ind = 1:length(Dir_e);
    figure; hold on;
    plot(ind, Data.pPT, '-');
    plot(ind, Data.ViewTime, '.-');
    plot(ind(Data.Type == 2), Data.ViewTime(Data.Type == 2), 'o');
    plot(ind(Data.Type == 2 & abs(Dir_e') < 45), Data.ViewTime(Data.Type == 2 & abs(Dir_e') < 45), 'x');
    title('Trial Type 2: Symbolic Pre-Cue');
end

%% assign output:
data_indiv.p_bins = {p_bin0, p_bin1, p_bin2};
data_indiv.z_bins = {pt_z_bin0, pt_z_bin2};
data_indiv.z_bin_sd = {pt_z_bin0_sd, pt_z_bin2_sd};
data_indiv.min_pt = pt_s(pt_m);
data_indiv.Dir_e = Dir_e;
data_indiv.Dir_a = Dir_a;
data_indiv.ViewTime = Data.ViewTime;
data_indiv.types = {type0, type1, type2, type3, type4};
data_indiv.kinematics = {Kin_x, Kin_y};
data_indiv.kin_summary = {kin_PV, kin_MT, kin_VAR};
data_indiv.targ_len = TARG_LEN;

varargout{1} = error_queue;