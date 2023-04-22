
function [scenarios,scen_legend,thresholds_out,var_names_dec] = decomposition (output, inputs, SI, dec_limit, manual_vars, manual_states, manual_thresholds, threshold_type, var_names)
% decomposition created the scenarios and maps them onto output values.
%
% * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
% INPUTS
%
%   output            - target variable (Y), size [N_runs, 1]
%   inputs            - input variables (Xs), size [N_runs, N_factors]
%   SI                - sensitivity indices
%   dec_limit         - threshold of cumulative significance for selection of 
%                       input variables for decomposition 
%   manual_vars       - [optional] for custom decomposition specify the order
%                       of variables for decomposition, use zero to exclude.
%                       For example, if 4 input variables, third and second 
%                       are desired for decomposition, then 
%                       manual_vars = [0 2 1 0].
%   manual_states     - [optional] the number of states for each input
%                       variable, i.e. [0 3 2 0]
%   manual_thresholds - [optional] maximums (numeric thresholds) 
%                       of every state, leave the rest as NaN, e.g. 
%                                          [NaN   3  -1    NaN;
%                                           NaN   5   0    NaN;
%                                           NaN   7   NaN  NaN]
%   threshold_type    - 1 for 'precentile-based' (same amount of observations in each state),   
%                       2 for 'median-based' (equaly-spaced ranges of states)
%   var_names         - cell array of input variables names
%
% OUTPUTS
%
%   scenarios      - an array of the same size as Y with scenario indices 
%                    for every simulation run
%   scen_legend    - a scenario table that shows which states of which
%                    variables compose different sceanrios
%   thresholds_out - thresholds of states of input variables
%   var_names_dec  - a cell array with sorted input variables' names
%
% * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
% Written by Mariia Kozlova, last updated 27.4.2023
% Many thanks for the grant #220178 from Finnish Foundation for Economic
% Education (lsr.fi) and the grant #6713/31/2021 from Business Finland.
% * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *


%% 1. Variables for decomposition

% Deciding on number and order of variables for decomposition

if isempty(manual_vars) % automatic selection of input variables
    [SI_sorted,var_order]=sort(SI,1,"descend"); % sorted SI and indices of most influential variables
    N_var_dec = find(cumsum(SI_sorted) > dec_limit, 1 ); % how many variables for decomposition
    var_order (N_var_dec+1:end) = 0;
else
    var_order = zeros(1, length(manual_vars));
    for i = 1:length(manual_vars) % modification of the variables selection array into a convenient form
            if manual_vars(i) > 0
                var_order(manual_vars(i)) = i;
            end
    end
    N_var_dec = nnz(var_order);
end

% Sorting var_names
var_names_dec = {};
    for f = 1 : N_var_dec
        var_names_dec(f) = var_names(var_order(f));
    end


%% 2. States formation
N_var=size(inputs,2);

 
    if isempty(manual_states) && isempty(manual_thresholds) % automatic definition of states
        states = zeros(1,N_var);
        
        if N_var_dec < 3      % assigning default number of states
            states(1:N_var) = 3;
        else
            states(1:N_var) = 2;
        end
        
        % checking for categorical

        for f = 1 : N_var
            n_unique = numel(unique(inputs(:,f)));
                if n_unique < 5
                    states(f) = n_unique; % number of states == number of categories
                end
        end

    else
        states = manual_states;
    end

    % zeroing out not participating in decomposition 
    for f = 1 : N_var
        if ismember(f,var_order(var_order~=0))
        else
            states(f) = 0;
        end
    end
    
    


%% 3. Numeric thresholds
N_runs = size(output,1);

if isempty(manual_thresholds) % if automatic definition of thresholds
    thresholds = NaN (max(states),N_var); 
    if threshold_type == 1 % percentile-based thresholds
            for f = 1 : N_var
                if states(f) == 0 % skip variables that are not for decomposition
                    f = f + 1;
                else
                    x = inputs(:,f); 
                    x_sorted = zeros (size(x,1)+1,1);
                    x_sorted(1:end-1) = sort (x); x_sorted(end)=x_sorted(end-1)+1;
                    min_threshold = x_sorted(1);
                    state_size = round(N_runs/states(f));
                    for s = 1 : states(f)-1
                        thresholds(s,f) = x_sorted (state_size*s+1);
                    end
                    thresholds(states(f),f) = max(x) + 1;
                end
            end
    
    else % median-based thresholds

            for f = 1 : N_var
                if states(f) ~= 0
                f_min = min (inputs(:,f));
                f_max = max (inputs(:,f));
                n_states = states(f);
                step = (f_max - f_min) / n_states;
                    for s = 1 : n_states
                        thresholds(s,f) = f_min + step * s;
                    end
                    thresholds(s,f) = max(inputs(:,f)) + 1; % shifting max of the last state upward to include all datapoints
                end
            end
    end
else
    thresholds = manual_thresholds(2:end,:);
end


%% 5. Scenario matrix
 N_scen = prod(nonzeros(states)); % number of scenarios
  
 Scen_matrix = ones(N_scen,N_var_dec+1);
 Scen_matrix(:,1) = [1:N_scen];

 j = N_scen;
 
    for v = 1 : N_var_dec
        ind_inf = var_order(v); % index of the most influential variable
        j = j / states(ind_inf);
        s = 1;
        
        for row = 1 : N_scen 
            Scen_matrix(row,v+1) = s;
            
            if rem(row,j) == 0 % index update
                s = s + 1;
                if s > states(ind_inf)
                    s = 1;
                end
             end
         end
     end
            
%% 6. Scenario matching

states_matching = zeros(N_runs,N_var_dec);
scenario_matching = zeros(N_runs,1);

for i = 1: N_runs
    
    for v = 1 : N_var_dec
        
        ind_inf = var_order(v); % picking variables in order of influence
        
        for s = 1 : states(ind_inf)

            if inputs(i,ind_inf) <  thresholds (s,ind_inf)
                states_matching(i,v) = s;
                break
            end 
        end
    end

    [q, scenario_matching(i)] = ismember(states_matching(i,:),Scen_matrix(:,2:end),'rows');
end


scenarios = scenario_matching;

% Adding min values to thresholds for export

    if isempty(manual_thresholds) % if automatic definition of thresholds
        
        thresholds_out = NaN(max(states)+1,N_var);
            
        for f = 1 : N_var_dec
            thresholds_out(1,var_order(f)) = min(inputs(:,var_order(f)));
            th = thresholds(states(var_order(f)),var_order(f));
            thresholds(states(var_order(f)),var_order(f)) = th - 1;
        end
       
        thresholds_out(2:end,:) = thresholds;
    else
        thresholds_out = thresholds;
    end


%% 7. Filling scenario legend

scen_legend = [Scen_matrix, NaN(N_scen,4)];

for sc = 1 : N_scen
    y = output(scenarios==sc);

    if isempty(y)
        
    else
        scen_legend(sc,end-3) = min(y);
        scen_legend(sc,end-2) = mean(y);
        scen_legend(sc,end-1) = max(y);
        scen_legend(sc,end) = size(y,1)/N_runs;
    end

end

end
