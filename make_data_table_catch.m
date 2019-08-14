function table_out = make_data_table_catch(table_in, varargin)

% function table_out = make_data_table(table_in, varargin)
%
% Function to organize behavioral data into a matrix of dimensions n_bins 

if nargin > 1
    beh_index = varargin{1};
else
    beh_index = 2; %index corresponds to the behavioral variable mean, not standard deviation.
    % use 3 for the standard deviation.
end
table_out = [...
    repmat(table_in{1}(:), 5, 1),... %PT
    reshape(table_in{beh_index}(:,:,:), numel(table_in{1})*5, 1), ... % behavior
    [zeros(numel(table_in{1}), 1); ...
    ones(numel(table_in{1}), 1);...
    2*ones(numel(table_in{1}), 1);...
    3*ones(numel(table_in{1}), 1);...
    4*ones(numel(table_in{1}), 1)],...% trial type
    repmat((1:size(table_in{1},1))', 5*size(table_in{1},2), 1)]; % participant