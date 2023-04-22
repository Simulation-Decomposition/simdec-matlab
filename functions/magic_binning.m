function [bin_idx, n_bins_out] = magic_binning(X, n_bins_default)
% bin_idx bins an array X into number of bins equal or smaller than n_bins and returns the corresponding indices. 
%
% The following binning rules apply:
% 1. The size of each bin should be no less than length(X)/n_bins.
% 2. The same values of X should be allocated to a single bin. 
% 3. Nan values should get the 0 bin index. 
%
%
% * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
% INPUTS
%   X               - an array of values
%   n_bins_default  - default number of bins
% 
% OUTPUT
%   bin_idx         - an array with bin indices of the same length as X
%   n_bins_out      - resulting number of bins (might be less than default if many NaNs or same values) 
%
% * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
% EXAMPLE
%   X = [NaN NaN 6 6 6 6 6 6 100 200 2 1]; % length(X) = 12
%   n_bins = 6;
%   bin_idx(X, n_bins) 
%              = [0 0 3 3 3 3 3 3 2 2 1 1]
%
%
% * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
% Written by Mariia Kozlova & ChatGPT, last updated 27.3.2023
% Many thanks for the grant #220178 from Finnish Foundation for Economic
% Education (lsr.fi)
% * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *


% Sort X
[X_sorted, idx] = sort(X);

% Compute minimum bin size
min_bin_size = floor(length(X_sorted)/n_bins_default);

% Initialize bin index array
bin_idx = zeros(size(X_sorted));

% If number of unique values is less than or equal to the number of bins,
% then use them to form the bins
if numel(unique(X_sorted)) <= n_bins_default
    [~, ~, bin_idx] = unique(X_sorted);
else

    % initialize remaining size & edge index
    not_nan_size = length(X) - sum(isnan(X)); % fixed
    remaining_size = length(X) - sum(isnan(X)); % decreases in every iteration of for loop
    current_edge_idx = min_bin_size; % increases in every iteration of for loop

    % Iterate over bins to produce a corresponding bin_idx
    for b = 1:n_bins_default
        current_bin_size = min_bin_size; % offsets in every iteration of for loop

        % While the edge is between same values, move one element further
        while current_edge_idx < not_nan_size && X_sorted(current_edge_idx + 1) == X_sorted(current_edge_idx)
            current_edge_idx = current_edge_idx + 1;
            current_bin_size = current_bin_size + 1;
        end

        % Assign bin indices to elements in the bin
        bin_idx(current_edge_idx - current_bin_size + 1 : current_edge_idx) = b;

        % Update remaining size and bin edge index for next iteration
        remaining_size = remaining_size - current_bin_size;

            % If not enough elements left for two bins, assign all to the
            % last bin and brerak the for loop.
            if remaining_size < min_bin_size * 2
                bin_idx(current_edge_idx + 1 : end) = b + 1;
                break
            end
        current_edge_idx = current_edge_idx + min_bin_size;
    end
end

% Map bin indices back to original order
bin_idx(idx) = bin_idx;

% Assign 0 bin index to NaN values
bin_idx(isnan(X)) = 0;

% The resulting number of bins
n_bins_out = max(bin_idx);

end
