
function [SI, FOE, SOE] = significance (output, inputs)
% significance calculates how much variability of the output is explained by inputs 
% 
% Uses functions number_of_bins, bin_data, magic_binning. 
%
% * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
% INPUTS
%   output - Target variable (Y), size [N_runs, 1]
%   inputs - Input variables (Xs), size [N_runs, N_factors]
%
% OUTPUTS
%   FOE - first-order effects (also called 'main' or 'individual'), size [N_factors, 1]  
%   SOE - second-order effects (also called 'interaction'), size [N_factors, N_factors]
%   SI  - significance index, combined effect of each input, size [N_factors, 1]  
%
% * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
% Written by Mariia Kozlova & ChatGPT, last updated 16.4.2023
% Many thanks for the grant #220178 from Finnish Foundation for Economic
% Education (lsr.fi) and the grant #6713/31/2021 from Business Finland.
% * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *


% Get data size

N_runs = size(output,1);
N_factors = size(inputs,2);


% Define number of bins
[N_bins_foe, N_bins_soe] = number_of_bins(N_runs, N_factors); 

% Initialize effects
FOE = NaN(N_factors,1);
SOE = zeros(N_factors);
SI = NaN(N_factors,1); % Combined effect

% Overall variance of the output
VarY = var(output);

% For each input variable 
for i = 1 : N_factors

    % Calculate its first-order effect
    [bin_avg, bin_count] = bin_data_1D(inputs(:,i), output, N_bins_foe); % compute average Y and its count in every bin
    FOE(i) = var(bin_avg, bin_count) / VarY; % calculate weighted variance of those and divide by the overal variance of the ouptut


        % Calculate second-order effect 
        for j = 1 : N_factors
            if i ~= j && j > i % for each unique pair of input variables
                
                [bin_avg_ij, bin_count_ij, bin_avg_i, bin_count_i, bin_avg_j, bin_count_j] ...
                    = bin_data_2D(inputs(:,i),inputs(:,j),output,N_bins_soe);
                var_ij = var(bin_avg_ij, bin_count_ij); % Var(E(Y|XiXj), calculation over 2D bins
                var_i = var(bin_avg_i, bin_count_i); % Var(E(Y|Xi) calculated with n_bins_soe 
                var_j = var(bin_avg_j, bin_count_j); % Var(E(Y|Xj) calculated with n_bins_soe 

                SOE(i,j) = (var_ij - var_i - var_j) / VarY;

            end
        end
end

    
% Calcualte combined effect (FOE plus halves of all SOE)
SOE_full = SOE + SOE.';
for k = 1 : N_factors
    SI(k) = FOE(k) + sum(SOE_full(:,k))/2;
end

end


