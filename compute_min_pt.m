function [min_pt, varargout] = compute_min_pt(pt, type, dir_err)
% function min_pt = compute_min_pt(pt, type, dir_err)
%
% Find the minimum preparation time for an individual given their
% directional error vs. preparation time function (for no pre-cue trials).
%
% David Huberdeau, Dec. 29, 2018

global EARLIEST_VALID_PT LATEST_VALID_PT SUCCESS_TH_ANGLE
func_tolerance = .05; %rather than taking the hard maximum of the function, 
 % take latest point at which the function is within 7% of the maximum.
 % Experience holds that this is more accurate and robust.

valid_vt = pt > EARLIEST_VALID_PT & pt < LATEST_VALID_PT;

pt0 = pt(valid_vt & type == 0);
de0 = dir_err(valid_vt & type == 0);
pe0 = de0 < SUCCESS_TH_ANGLE & de0 > -SUCCESS_TH_ANGLE;

pt_s = linspace(min(pt0), max(pt0), 30);
fun = nan(size(pt_s));
fun_corr = nan(size(pt_s));
fun_err = nan(size(pt_s));
for i_s = 1:length(pt_s)
    fun(i_s) = sum(pe0(pt0 > pt_s(i_s))) + sum(1 - pe0(pt0 <= pt_s(i_s)));
    fun_corr(i_s) = sum(1 - pe0(pt0 <= pt_s(i_s)));
    fun_err(i_s) = sum(pe0(pt0 > pt_s(i_s)));
end

[fun_m, pt_m] = max(fun);
closeness_to_maximum = fun_m*func_tolerance;
if sum(abs(fun - fun_m) < closeness_to_maximum) > 0 % there is a point nearly
                                                    % equal to the max
    % take the later point's index as the min-PT
    temp_inds = 1:length(fun);
    equal_max_inds = temp_inds(abs(fun - fun_m) < closeness_to_maximum);
    pt_m = max(equal_max_inds); %take the latest index among (near) equals
end

min_pt = pt_s(pt_m);
varargout = {pt_s, fun, pt_m, fun_corr, fun_err};