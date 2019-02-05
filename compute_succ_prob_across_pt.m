function [sp0, sp1, sp2] = compute_succ_prob_across_pt(pt, dir_err, type)
% function [sp0, sp1, sp2] = compute_succ_prob_across_pt(pt, dir_err, type)
%
% Compute the probability of a movement success within bins of pt for each
% trial type given the preparation times (pt), directional errors
% (dir_err), and trial types (type).
%
% David Huberdeau, 12/21/2018

global n_bins EARLIEST_VALID_PT LATEST_VALID_PT SUCCESS_TH_ANGLE

hist_bins = linspace(EARLIEST_VALID_PT, LATEST_VALID_PT, n_bins+1);

[n0, edge0, bins0] = histcounts(pt(type == 0), hist_bins);
[n1, edge1, bins1] = histcounts(pt(type == 1), hist_bins);
[n2, edge2, bins2] = histcounts(pt(type == 2), hist_bins);

de0 = dir_err(type == 0);
sp0 = nan(size(n0));
for i_bin = 1:(length(edge0) - 1)
    sp0(i_bin) = ...
        sum(abs(de0(bins0 == i_bin)) < SUCCESS_TH_ANGLE)/sum(bins0 == i_bin);
end

de1 = dir_err(type == 1);
sp1 = nan(size(n1));
for i_bin = 1:(length(edge1) - 1)
    sp1(i_bin) = ...
        sum(abs(de1(bins1 == i_bin)) < SUCCESS_TH_ANGLE)/sum(bins1 == i_bin);
end

de2 = dir_err(type == 2);
sp2 = nan(size(n2));
for i_bin = 1:(length(edge2) - 1)
    sp2(i_bin) = ...
        sum(abs(de2(bins2 == i_bin)) < SUCCESS_TH_ANGLE)/sum(bins2 == i_bin);
end
