function [bin_avg, bin_count] = bin_data_1D(X, Y, n_bins_default)
% bin_data computes averages and counts of Y in bins of X 
%
% Uses function magic_binning.
%
% * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
% INPUTS
%   X               - a 1D array of X values
%   Y               - a 1D array of Y values
%   n_bins_default  - default number of bins
% 
% OUTPUT
%   bin_avg         - averages of Y values in bins of X
%   bin_count       - number of Y values in bins of X
%
% * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
% Written by Mariia Kozlova, last updated 16.4.2023
% Many thanks for the grant #220178 from Finnish Foundation for Economic
% Education (lsr.fi) and the grant #6713/31/2021 from Business Finland.
% * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *


    % Getting the bins indices
    [bin_idx, n_bins_X] = magic_binning(X, n_bins_default);

    % Initializing matrices for recording averages and counts of Y 
    bin_avg = NaN(n_bins_X, 1);
    bin_count = NaN(n_bins_X, 1);

    % Computing averages and counts of Y     
    for b = 1:n_bins_X
        bin_avg(b) = mean(Y(bin_idx==b));
        bin_count(b) = length(Y(bin_idx==b));
    end

end