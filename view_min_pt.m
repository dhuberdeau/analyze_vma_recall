% view all of the subject's (from a particular group) minimum pt and the
% function used to compute it. Mainly used to select example subject for
% supplemental figure.
%
% Must run this function after analyze_statLearn_v2.m
%
% David Huberdeau

grp = 1;

figure;
for i_sub = 1:length(subject_list_3T)
    this_pt = pt_all{grp}(type_all{grp}(:, i_sub) ==0, i_sub);
    this_type = type_all{grp}(type_all{grp}(:, i_sub) ==0, i_sub);
    this_de = dir_err_all{grp}(type_all{grp}(:, i_sub) ==0, i_sub);
    
%     [min_pt, pt_s, fun, pt_m] = compute_min_pt(this_pt, this_type, this_de);

    min_pt = min_pt_all(i_sub, grp);
    
    subplot(4,5,i_sub); hold on;
    plot(this_pt, this_de, 'r.');
    plot(pt_s, fun, 'g-');
    plot([min_pt, min_pt], [-200 200]);
    axis([-.2 .85 -200 200])
    
end
