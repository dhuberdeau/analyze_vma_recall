function table_out = make_data_table(table_in, varargin)

% function table_out = make_data_table(table_in, varargin)
%
% Function to organize behavioral data into a matrix of dimensions n_bins 
%
% INPUT: table_in - a cell array with elements {pt_ind, behavior_pt, var_pt, n_samples_pt};
%        behavior_index (optional) - 2 = mean, 3 = standard deviation
%
% OUTPUT: formated table for linear mixed effects analysis.


if nargin > 1
    beh_index = varargin{1};
else
    beh_index = 2; %index corresponds to the behavioral variable mean, not standard deviation.
    % use 3 for the standard deviation.
end
table_out = [...
    repmat(table_in{1}(:), 3, 1),... %PT
    reshape(table_in{beh_index}(:,:,1:3), numel(table_in{1})*3, 1), ... % behavior
    [zeros(numel(table_in{1}), 1); ...
    ones(numel(table_in{1}), 1);...
    2*ones(numel(table_in{1}), 1)],...% trial type
    repmat((1:size(table_in{1},1))', 3*size(table_in{1},2), 1)]; % participant