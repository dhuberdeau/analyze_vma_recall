function [recall_prob_lowPT_2, recall_prob_highPT_2] = ...
    compute_recall_probability_appearance_order(...
    mov_class, target, type, pt, min_pt, varargin)
% recall_prob = compute_recall_probability_appearance_order(...
%               mov_class, target, type, pt, min_pt, [trial_type_to_display])
%
% Computes the probability of recall of each symbol (in target input
% variable) as a function of the appearance of that symbol.
%
% Inputs:
%   mov_class - the tareget moved towards on each trial (trs x subjs)
%   target - the actual symbol of movement
%   type - the type of each trial
%   pt - the actual pt of each trial
%   the minimum pt measured from each person
%   trial_type_to_display (optional) - symbolic (2), 
%     direct (1), or no precue (0) [default = 2, symbolic];
%
% David Huberdeau, Dec. 10, 2018

% target_ = nan(
% for i_sub = 1:size(target,2)
%     sub_type = type(:, i_sub);
%     sub_mov_class = mov_class(sub_type == 2, i_sub);
%     sub_target = target(sub_type == 2, i_sub);
%     sub_target_num = target_num(sub_type == 2, i_sub);
%     sub_pt = pt(sub_type == 2, i_sub);
%     sub_mpt = min_pt(i_sub);
% end

if nargin > 5
    search_type = varargin{1};
else
    search_type = 2;
end

N_TYPE = max(sum(type == search_type,1));
mov_class_ = nan(N_TYPE, size(target,2));
target_ = nan(N_TYPE, size(target,2));
pt_ = nan(N_TYPE, size(target,2));
for i_sub = 1:size(target,2)
    sub_type = type(:, i_sub);
    sub_mov_class = mov_class(sub_type == search_type, i_sub);
    sub_target = target(sub_type == search_type, i_sub);
    sub_pt = pt(sub_type == search_type, i_sub);
    
    mov_class_(1:length(sub_mov_class), i_sub) = sub_mov_class;
    target_(1:length(sub_target), i_sub) = sub_target;
    pt_(1:length(sub_pt), i_sub) = sub_pt;
end
mov_class = mov_class_;
target = target_;
pt = pt_;

unique_syms = unique(target(~isnan(target(:,1)), 1));
unique_targs = unique(mov_class(~isnan(mov_class(:,1)), 1));
N_SYMS = length(unique_syms);
targ_appearances = nan(1, N_SYMS);
for i_targ = 1:N_SYMS
    targ_appearances(i_targ) = ...
        max(sum(target == unique_syms(i_targ)));
end
N_APPEAR = max(targ_appearances);
N_SUBS = size(mov_class, 2);
N_TARGS = length(unique_targs);
target_num = target;
for i_sub = 1:size(target,2)
    if N_TARGS == 3
        target_num(target(:, i_sub) == 4, i_sub) = 1;
        target_num(target(:, i_sub) == 5, i_sub) = 2;
        target_num(target(:, i_sub) == 6, i_sub) = 3;
    elseif N_TARGS == 4
        target_num(target(:, i_sub) == 5, i_sub) = 1;
        target_num(target(:, i_sub) == 6, i_sub) = 2;
        target_num(target(:, i_sub) == 7, i_sub) = 3;
        target_num(target(:, i_sub) == 8, i_sub) = 4;
        target_num(target(:, i_sub) == 9, i_sub) = 1;
        target_num(target(:, i_sub) == 10, i_sub) = 2;
        target_num(target(:, i_sub) == 11, i_sub) = 3;
        target_num(target(:, i_sub) == 12, i_sub) = 4;
    elseif N_TARGS == 6
        target_num(target(:, i_sub) == 7, i_sub) = 1;
        target_num(target(:, i_sub) == 8, i_sub) = 2;
        target_num(target(:, i_sub) == 9, i_sub) = 3;
        target_num(target(:, i_sub) == 10, i_sub) = 4;
        target_num(target(:, i_sub) == 11, i_sub) = 5;
        target_num(target(:, i_sub) == 12, i_sub) = 6;
    else
       error('Invalid number of targets specified.') 
    end
end

% for type 0 (no pre-cue trials) only targets 1,.. N were used (as opposed
% to 2*N or 3*N when multiple symbols corresponded to each target, thus,
% also re-assign targets greater than N)
if search_type == 0
    for i_sub = 1:size(target,2)
        if N_TARGS == 3
            target(target(:, i_sub) == 4, i_sub) = 1;
            target(target(:, i_sub) == 5, i_sub) = 2;
            target(target(:, i_sub) == 6, i_sub) = 3;
        elseif N_TARGS == 4
            target(target(:, i_sub) == 5, i_sub) = 1;
            target(target(:, i_sub) == 6, i_sub) = 2;
            target(target(:, i_sub) == 7, i_sub) = 3;
            target(target(:, i_sub) == 8, i_sub) = 4;
            target(target(:, i_sub) == 9, i_sub) = 1;
            target(target(:, i_sub) == 10, i_sub) = 2;
            target(target(:, i_sub) == 11, i_sub) = 3;
            target(target(:, i_sub) == 12, i_sub) = 4;
        elseif N_TARGS == 6
            target(target(:, i_sub) == 7, i_sub) = 1;
            target(target(:, i_sub) == 8, i_sub) = 2;
            target(target(:, i_sub) == 9, i_sub) = 3;
            target(target(:, i_sub) == 10, i_sub) = 4;
            target(target(:, i_sub) == 11, i_sub) = 5;
            target(target(:, i_sub) == 12, i_sub) = 6;
        else
           error('Invalid number of targets specified.') 
        end
    end
end

recall_prob_lowPT_2 = nan(N_SYMS, N_APPEAR, N_SUBS);
recall_prob_highPT_2 = nan(N_SYMS, N_APPEAR, N_SUBS);

for i_sub = 1:N_SUBS
%     sub_type = type(:, i_sub);
%     sub_mov_class = mov_class(sub_type == 2, i_sub);
%     sub_target = target(sub_type == 2, i_sub);
%     sub_target_num = target_num(sub_type == 2, i_sub);
%     sub_pt = pt(sub_type == 2, i_sub);
%     sub_mpt = min_pt(i_sub);

%     sub_type = type(:, i_sub);
    sub_mov_class = mov_class(:, i_sub);
    sub_target = target(:, i_sub);
    sub_target_num = target_num(:, i_sub);
    sub_pt = pt(:, i_sub);
    sub_mpt = min_pt(i_sub);
    
    % appropriately truncate data:
%     sub_type = sub_type(~isnan(sub_target));
    sub_mov_class = sub_mov_class(~isnan(sub_target));
    sub_target_num = sub_target_num(~isnan(sub_target));
    sub_pt = sub_pt(~isnan(sub_target));
    sub_target = sub_target(~isnan(sub_target));
    
    symbol_counter = ones(N_SYMS, 1);
    
    for i_tr = 1:length(sub_target_num)
        this_target_num = sub_target_num(i_tr);
        this_mov_class = sub_mov_class(i_tr);
        this_symbol = sub_target(i_tr);
        this_pt = sub_pt(i_tr);
%         this_type = sub_type(i_tr);
        
        if ~isnan(this_mov_class)
            move_bool = this_target_num == this_mov_class;
        else
            move_bool = nan;
        end
        try
        if ~isnan(this_symbol)
            if this_pt < sub_mpt %&& this_type == 2
                % lowPT trial
                recall_prob_lowPT_2(this_symbol, symbol_counter(this_symbol), i_sub) = ...
                    move_bool;
%                 symbol_counter(this_symbol) = symbol_counter(this_symbol) + 1;
            elseif this_pt >= sub_mpt %&& this_type == 2
                % highPT trial
                recall_prob_highPT_2(this_symbol, symbol_counter(this_symbol), i_sub) = ...
                    move_bool;
%                 symbol_counter(this_symbol) = symbol_counter(this_symbol) + 1;
            else
                % must be a different trial type.. don't record these.
            end
            symbol_counter(this_symbol) = symbol_counter(this_symbol) + 1;
        else
            error('A symbol is inappropriately NaN');
        end
        catch err
            rethrow(err);
        end
    end
end

% for i_sub = 1:size(mov_class, 2)
%     sub_mov_class_ = mov_class(:, i_sub);
%     sub_target_ = target(:, i_sub);
%     sub_target_num_ = target_num(:, i_sub);
% %     sub_trial_ = 1:length(sub_target_);
%     
%     sub_pt = pt(:, i_sub);
%     sub_mpt = min_pt(i_sub);
%     
%     sub_mov_class = sub_mov_class_(sub_pt < sub_mpt);
%     sub_target = sub_target_(sub_pt < sub_mpt);
%     sub_target_num = sub_target_num_(sub_pt < sub_mpt);
% %     sub_trial = sub_trial_(sub_pt < sub_mpt);
%     sub_trial = 1:length(sub_target);
%     
%     for i_targ = 1:N_SYMS
%         targ_trial = sub_trial(sub_target == i_targ);
%         for i_tr = 1:length(targ_trial)
%             try
%             recall_prob_lowPT(i_targ, i_tr, i_sub) = ...
%                 sub_mov_class(targ_trial(i_tr)) == sub_target_num(targ_trial(i_tr));
%             catch
%                 pause;
%             end
%         end
%     end
% end

