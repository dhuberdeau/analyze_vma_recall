function table_out = plot_behavior_over_pt(view_time, behavior_var, trial_types, varargin)
% function to plot the desired behavioral variable across different
% preparation times ('view_time').
%
% INPUTS:
%  - view_time          (n_trials x n_subs) - PT
%  - behavior_var       (n_trials x n_subs) - relevant behavioral variable
%  - trial_types        (n_trials x n_subs) - trial types
%  - n_bins (optional)  scalar              - how many bins of PT desired
%
% OUTPUTS:
%  - table_out 
%
% David Huberdeau

%% define global variables:
global EARLIEST_VALID_PT LATEST_VALID_PT MIN_N_PT_FOR_MEASURE RED_COLOR...
    GREEN_COLOR BLUE_COLOR

%% handle inputs:
if nargin > 3
    n_bins = varargin{1};
else
    n_bins = 8;
end

N_SUB = size(view_time,2);
N_TYPES = 5;

%% define desired output structures:
pt_ind = nan(N_SUB, n_bins);
behavior_pt = nan(N_SUB, n_bins, N_TYPES);
n_samples_pt = nan(N_SUB, n_bins, N_TYPES);
var_pt = nan(N_SUB, n_bins, N_TYPES);

type_codes = [0, 1, 2, 3, 4];
%%
for i_sub = 1:N_SUB
    %% Split trials into types / exclude invalid trials:
    valid_vt = view_time(:, i_sub) > EARLIEST_VALID_PT & ...
        view_time(:, i_sub) < LATEST_VALID_PT;
    vt = view_time(valid_vt, i_sub);
    
    hist_bins = linspace(EARLIEST_VALID_PT, LATEST_VALID_PT, n_bins+1);
    [n_pt_all, edge_pt_all] = histcounts(vt, hist_bins);
    pt_ind(i_sub, :)  = edge_pt_all(2:end) - diff(edge_pt_all)/2;
    
    pt_ = cell(1, N_TYPES);
    beh_ = cell(1, N_TYPES);
    for i_type = 1:N_TYPES
        pt_{i_type} = view_time(valid_vt & trial_types(:, i_sub) == type_codes(i_type), i_sub);
        beh_{i_type} = behavior_var(valid_vt & trial_types(:, i_sub) == type_codes(i_type), i_sub);
        
        for i_edge = 2:length(edge_pt_all)
            this_pt_bin = pt_{i_type} > edge_pt_all(i_edge - 1) & ...
                pt_{i_type} <= edge_pt_all(i_edge);
            
            if sum(this_pt_bin) >= MIN_N_PT_FOR_MEASURE
                behavior_pt(i_sub, i_edge - 1, i_type) = nanmean(beh_{i_type}(this_pt_bin));
                var_pt(i_sub, i_edge - 1, i_type) = nanstd(beh_{i_type}(this_pt_bin));
                n_samples_pt(i_sub, i_edge - 1, i_type) = sum(this_pt_bin);
            end
        end
    end
    
end

%% plot output
figure; hold on;
errorbar(nanmean(pt_ind) + 0.02, nanmean(behavior_pt(:, :, 1)),...
sqrt(nanvar(behavior_pt(:, :, 1))./sum(~isnan(behavior_pt(:, :, 1)))), ...
'.-', 'LineWidth', 2, 'MarkerSize', 20, 'Color', RED_COLOR);

errorbar(nanmean(pt_ind) - 0.02, nanmean(behavior_pt(:, :, 2)),...
sqrt(nanvar(behavior_pt(:, :, 2))./sum(~isnan(behavior_pt(:, :, 2)))), ...
'.-', 'LineWidth', 2, 'MarkerSize', 20, 'Color', GREEN_COLOR);

errorbar(nanmean(pt_ind), nanmean(behavior_pt(:, :, 3)),...
sqrt(nanvar(behavior_pt(:, :, 3))./sum(~isnan(behavior_pt(:, :, 3)))), ...
'.-',  'LineWidth', 2, 'MarkerSize', 20, 'Color', BLUE_COLOR);

% errorbar(nanmean(pt_ind), nanmean(behavior_pt(:, :, 4)),...
% sqrt(nanvar(behavior_pt(:, :, 4))./sum(~isnan(behavior_pt(:, :, 4)))), ...
% 'o-',  'LineWidth', 2, 'MarkerSize', 12, 'Color', GREEN_COLOR);
% 
% errorbar(nanmean(pt_ind), nanmean(behavior_pt(:, :, 5)),...
% sqrt(nanvar(behavior_pt(:, :, 5))./sum(~isnan(behavior_pt(:, :, 5)))), ...
% 'o-',  'LineWidth', 2, 'MarkerSize', 12, 'Color', BLUE_COLOR);

axis([EARLIEST_VALID_PT, LATEST_VALID_PT, min(min(min(behavior_pt))), max(max(max(behavior_pt)))])
set(gca,'fontsize',18)
legend('None', 'Direct', 'Symbolic', 'location', 'bestoutside')

%% assign output
table_out = {pt_ind, behavior_pt, var_pt, n_samples_pt};