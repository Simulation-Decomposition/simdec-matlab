function [bin_avg_ij_string, bin_count_ij_string, bin_avg_i, bin_count_i, bin_avg_j, bin_count_j] = bin_data_2D(Xi,Xj,Y,n_bins_default)
% bin_data_2D computes averages and counts of Y in bins of XiXj, Xi, and Xj 
%
% Uses function magic_binning.
%
% * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
% INPUTS
%   Xi              - a 1D array of Xi values
%   Xj              - a 1D array of Xj values
%   Y               - an array of Y values
%   n_bins_default  - default number of bins
% 
% OUTPUT
%   bin_avg_ij_string    - an array of average Y values in 2D bins of Xi and Xj     
%   bin_count_ij_string  - an array of corresponding count of Y values     
%   bin_avg_i            - an array of average Y values in 1D bins of Xi 
%   bin_count_i          - an array of corresponding count of Y values  
%   bin_avg_j            - an array of average Y values in 1D bins of Xj 
%   bin_count_j          - an array of corresponding count of Y values  
%
% * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
% Written by Mariia Kozlova with help of Leo Tran, last updated 16.4.2023
% Many thanks for the grant #220178 from Finnish Foundation for Economic
% Education (lsr.fi) and the grant #6713/31/2021 from Business Finland.
% * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *


    % Getting the bins indices
    [bin_idx_i, n_bins_i] = magic_binning(Xi, n_bins_default);
    [bin_idx_j, n_bins_j] = magic_binning(Xj, n_bins_default);


    % Initializing matrices for recording averages and counts of Y 

    bin_avg_ij = NaN(n_bins_i, n_bins_j);         
    bin_count_ij = NaN(n_bins_i, n_bins_j);           
    bin_avg_i = NaN(n_bins_i, 1);
    bin_count_i = NaN(n_bins_i, 1);
    bin_avg_j = NaN(n_bins_j, 1);
    bin_count_j = NaN(n_bins_j, 1);


    % Computing averages and counts of Y 

        % for Xi
        
        for n = 1 : n_bins_i
            bin_avg_i(n) = mean(Y(bin_idx_i==n));
            bin_count_i(n) = length(Y(bin_idx_i==n));
        end


        % for Xj

        for m = 1 : n_bins_j
            bin_avg_j(m) = mean(Y(bin_idx_j==m));
            bin_count_j(m) = length(Y(bin_idx_j==m));
        end

        
        % for XiXj

        for n = 1 : n_bins_i
            for m = 1 : n_bins_j
                mask_i = find(bin_idx_i==n);
                mask_j = find(bin_idx_j==m);
                [bin_indices,~] = intersect(mask_i,mask_j,'stable');
                bin_avg_ij(n,m) = mean(Y(bin_indices));         
                bin_count_ij(n,m) = length(Y(bin_indices));    
            end
        end


    % Converting 2D bin_avg and bin_count into 1D with no NaNs (empty bins) 
    bin_avg_ij_string = reshape(bin_avg_ij.',1,[]);
    ind_exists = find(~isnan(bin_avg_ij_string));
    bin_avg_ij_string = bin_avg_ij_string(ind_exists)';
    
    bin_count_ij_string = reshape(bin_count_ij.',1,[]);
    bin_count_ij_string = bin_count_ij_string(ind_exists)';

end