function dir_err = compute_relative_target_error(reach_dir, targets_plus)

% function dir_err = compute_relative_target_error(reach_dir)
%
% Compute the reach directional error relative to the apparently chosen
% target (the target that is closest to the reach direction) regardless of
% which target was presented. This will give a measure of reach deviance to
% a target direction even if the choice of target was incorrect.
%
% Inputs:
%   reach_dir = the absolute reach direction (N-trals x 1);
%   targets_plus = list of targets (with extra one at end for wrap around)
%
% Outputs:
%   dir_err - the relative directional error.
%
% David Huberdeau, NTB lab, July 13, 2019

dir_err = nan(length(reach_dir),1);
for i_tr = 1:length(reach_dir)
    target_errors = abs(reach_dir(i_tr) - targets_plus);
    [~, targ_min] = min(target_errors);
    
    dir_err(i_tr) = targets_plus(targ_min) - reach_dir(i_tr);
end