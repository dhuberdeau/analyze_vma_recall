function varargout = compute_variability_by_pt(view_time, dir_error, min_pt)
% function varargout = compute_var_across_pt(varargin)

reach_var_persub_0 = nan(size(view_time,2), 2); %subject x low-PT or high-PT
reach_var_persub_1 = nan(size(view_time,2), 2); %subject x low-PT or high-PT
reach_var_persub_2 = nan(size(view_time,2), 2); %subject x low-PT or high-PT

pt_low_bound = -.5;
pt_hgh_bound = .5;
N_min_samples = 3;
for i_sub = 1:size(view_time,2)
    % select out this subject's view time and directional error for each type
    this_view_time_0 = view_time(:, i_sub, 1) - min_pt(i_sub);
    this_direrr_0 = dir_error(:, i_sub, 1);
    
    this_view_time_1 = view_time(:, i_sub, 2) - min_pt(i_sub);
    this_direrr_1 = dir_error(:, i_sub, 2);
    
    this_view_time_2 = view_time(:, i_sub, 3) - min_pt(i_sub);
    this_direrr_2 = dir_error(:, i_sub, 3);
    
 % compute the variability (st dev) for subset of PT's in prespecified range
    this_inds_0 = this_view_time_0 > pt_low_bound & this_view_time_0 < 0 & abs(this_direrr_0) < 30;
    if sum(this_inds_0) > N_min_samples
        reach_var_persub_0(i_sub, 1) = sqrt(nanvar(this_direrr_0(this_inds_0)));
    end
    this_inds_0 = this_view_time_0 > 0 & this_view_time_0 < pt_hgh_bound;
    if sum(this_inds_0) > N_min_samples
        reach_var_persub_0(i_sub, 2) = sqrt(nanvar(this_direrr_0(this_inds_0)));
    end
    
    this_inds_1 = this_view_time_1 > pt_low_bound & this_view_time_1 < 0 & abs(this_direrr_1) < 30;
    if sum(this_inds_1) > N_min_samples
        reach_var_persub_1(i_sub, 1) = sqrt(nanvar(this_direrr_1(this_inds_1)));
    end
    this_inds_1 = this_view_time_1 > 0 & this_view_time_1 < pt_hgh_bound;
    if sum(this_inds_1) > N_min_samples
        reach_var_persub_1(i_sub, 2) = sqrt(nanvar(this_direrr_1(this_inds_1)));
    end
    
    this_inds_2 = this_view_time_2 > pt_low_bound & this_view_time_2 < 0 & abs(this_direrr_2) < 30;
    if sum(this_inds_2) > N_min_samples
        reach_var_persub_2(i_sub, 1) = sqrt(nanvar(this_direrr_2(this_inds_2)));
    end
    this_inds_2 = this_view_time_2 > 0 & this_view_time_2 < pt_hgh_bound;
    if sum(this_inds_2) > N_min_samples
        reach_var_persub_2(i_sub, 2) = sqrt(nanvar(this_direrr_2(this_inds_2)));
    end
end
varargout = {reach_var_persub_0, reach_var_persub_1, reach_var_persub_2};