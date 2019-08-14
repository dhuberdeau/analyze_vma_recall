% must have run analyze_retention_group_v1 first.


%re-assign catch trials as if they were no-cue trials:

% pt_catch = nan(size(pt_all));
% de_catch = nan(size(direrror_all));
% type_catch = nan(size(type_all));
min_pt_catch = nan(1, size(pt_all,2));

for i_sub = 1:size(pt_all,2)
    catch_trials = type_all(:, i_sub) == 3 | type_all(:, i_sub) == 4;
    
    temp_pt = pt_all(catch_trials, i_sub);
    temp_de = direrror_all(catch_trials, i_sub);
    temp_type = zeros(sum(catch_trials), 1);
    
    min_pt_catch(i_sub) = compute_min_pt(temp_pt,temp_type,temp_de);
    
%     pt_catch(1:length(temp_pt), i_sub) = temp_pt;
%     de_catch(1:length(temp_de), i_sub) = temp_de;
%     type_catch(1:length(temp_type), i_sub) = temp_type;
end

% min_pt_catch = compute_min_pt(pt_catch,type_catch,de_catch);