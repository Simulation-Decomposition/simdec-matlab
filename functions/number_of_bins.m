function [N_bins_foe, N_bins_soe] = number_of_bins(N_runs, N_factors)

% number_of_bins defines the optimal number of bins for first-order &
% second-order significance indices calculation. 
%
% * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
% INPUTS
%   N_runs     - the length of data / number of observations / number of
%              simulation runs
%   N_factors  - number of inputs variables
%
% OUTPUTS
%   N_bins_foe - number of bins for first-order effect calucation
%   N_bins_soe - number of bins for second-order effect calculation
%
% * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
% Written by Mariia Kozlova, last updated 5.3.2023
% Many thanks for the grant #220178 from Finnish Foundation for Economic
% Education (lsr.fi) and the grant #6713/31/2021 from Business Finland.
% * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *


    % for first-order effects (foe)
    N_bins_foe = ceil(36 - 2.7* N_factors + (0.0017 - 0.00008 * N_factors) * N_runs); % linear approximation of the experimental results from (Marzban & Lahmer, 2016)
        if N_bins_foe <= 30
            N_bins_foe = 10; % setting a limit to fit the experimental results
        end

% 
%         while rem(N_runs, N_bins_foe) ~= 0
%             N_bins_foe = N_bins_foe + 1;
%         end  

N_bins_foe

    % for second-order effects (soe)

    N_bins_soe = max(4,round(sqrt(N_bins_foe)));

end