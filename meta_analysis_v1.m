%% compare minimum PT across experiments (which had differing numbers of 
% targets) to see if the power-law holds / Hick's law (i.e. the tendency for choice RT
% to increase as a function of the log-number of options). 

% Collect Min-PT of no-precue trials:
min_pt_comb = nan(...
    max([size(min_pt_E1, 1), size(min_pt_E2,1)]),...
    size(min_pt_E2,2) + size(min_pt_E1,2));

% assign group min_pt into *comb matrix in order of target number:
min_pt_comb(1:size(min_pt_E2,1), 1:2) = min_pt_E2(:, [1,4]); %3-targets
min_pt_comb(1:size(min_pt_E1,1), 3) = min_pt_E1; %4-targets
min_pt_comb(1:size(min_pt_E2,1), 4) = min_pt_E2(:, 2); %4-targets
min_pt_comb(1:size(min_pt_E2,1), 5) = min_pt_E2(:, 3); %6-targets

targ_nums = [3*ones(size(min_pt_comb,1), 2), ...
    4*ones(size(min_pt_comb,1), 2), ...
    6*ones(size(min_pt_comb,1), 1)];

% plot:
f_mpt = figure; hold on;
min_pt_corr = min_pt_comb - .06; % correct the minimum PT.
errorbar([2.95 3.05 3.95 4.05 6], ...
    nanmean(min_pt_corr,1), ...
    sqrt(nanvar(min_pt_corr)./sum(~isnan(min_pt_corr))),...
    'ks', 'LineWidth', 2);
plot(targ_nums, min_pt_corr, 'k.');
axis([2.5 6.5 .15 .3])
ylabel('Min PT (sec)');
xlabel('Number of targets');
set(f_mpt, 'Position', [560, 700, 487 242])
saveas(f_mpt, 'MinPT_by_TargetNum.pdf')

min_pt_dat_mat = [min_pt_corr(:), targ_nums(:)];
csvwrite('min_pt_data', min_pt_dat_mat);