function [scenarios, scen_legend, boundaries_out, stacked_histogram] = simdec_visualization(output, inputs, SI, varargin)
% builds SimDec visualization using data decomposition
%
% * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
% INPUTS
%
%   REQUIRED
%   output               - target variable (Y), size [N_runs, 1]
%   inputs               - input variables (Xs), size [N_runs, N_factors]
%   SI                   - sensitivity indices
%
%   OPTIONAL
%   DecompositionLimit   - threshold of cumulative importance (sum(SI)) 
%                        for selection of input variables for decomposition 
%   OrderOfVariables     - for custom decomposition specify the order
%                        of variables for decomposition, use zero to exclude.
%                        For example, if 4 input variables, third and second 
%                        are desired for decomposition, then specify
%                        OrderOfVariables as [0 2 1 0].
%   NumberOfStates       - the number of states for each input
%                        variable, i.e. [0 3 2 0].
%   BoundaryType        - defines how the numerical boundaries between 
%                        the states of inputs are defined:
%                        'precentile-based' for same amount of observations 
%                        in each state,   
%                        'median-based' for equaly-spaced ranges of states.
%   StateBoundaries      - maximums (numeric boundaries) 
%                        of every state, leave the rest as NaN, e.g. 
%                                          [NaN   3  -1    NaN;
%                                           NaN   5   0    NaN;
%                                           NaN   7   NaN  NaN].
%   MainColors           - a cell array with HEX numbers of the main colors
%                        for decomposition (should correspond to the number
%                        of states of the first for decomposition input variable).  
%   OutputName           - name of the output variable.
%   InputNames           - names of the input variables in the order of
%                        their appearance in the original dataset. 
%
% 
% OUTPUTS
%
%   scenarios          - an array of the same size as Y with scenario indices 
%                      for every simulation run
%   scen_legend        - a scenario table that shows which states of which
%                      variables compose different sceanrios
%   boundaries_out     - numeric boundaries of states of input variables
%   stacked_histogram  - object that returns the visualization
%
% * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
% Written by Mariia Kozlova, last updated 5.10.2023
% Many thanks for the grant #220178 from Finnish Foundation for Economic
% Education (lsr.fi) and the grant #6713/31/2021 from Business Finland.
% * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *



%% Initializing inputs

    default_decomposition_limit = 0.8*sum(SI);  
    default_var_order = []; % default order of variables is dictated by sensitivity indices
    expected_size = [1, size(inputs, 2)];
    default_number_states = []; % default number of states is computed automatically
    default_boundaries = []; % default boundaries between states are computed automatically
    default_boundary_type = 'percentile-based'; % same amount of observations in each state
    expected_boundary_types = {'percentile-based','median-based'}; % 'median-based' for equaly-spaced ranges
    default_colors = {'#DC267F'; '#E8EA2F'; '#26DCD1';'#C552E4';'#3F45D0'};
    defult_output_name = 'Y';
    default_input_names = cell(1, size(inputs,2));
        for i = 1:size(inputs,2)
            default_input_names{i} = ['X' num2str(i)];
        end


   p = inputParser;
   validScalarPosNum = @(x) isnumeric(x) && isscalar(x) && (x > 0);
   addRequired(p,'output');
   addRequired(p,'inputs');
   addRequired(p,'SI');
   addOptional(p,'DecompositionLimit',default_decomposition_limit,validScalarPosNum);
   addOptional(p,'OrderOfVariables',default_var_order, @(x) isequal(size(x), expected_size));
   addOptional(p,'NumberOfStates',default_number_states, @(x) isequal(size(x), expected_size));
   addOptional(p,'BoundaryType',default_boundary_type,@(x) any(validatestring(x,expected_boundary_types)));
   addOptional(p,'StateBoundaries',default_boundaries);
   addParameter(p,'MainColors',default_colors,@iscellstr);
   addParameter(p,'OutputName',defult_output_name);
   addParameter(p,'InputNames',default_input_names,@iscellstr);
   
   
   parse(p,output,inputs,SI,varargin{:});
   


%% 1. Variables for decomposition

N_var=size(inputs,2);

% Deciding on number and order of variables for decomposition

if isempty(p.Results.OrderOfVariables) % automatic selection of input variables
    [SI_sorted,var_order]=sort(SI,1,"descend"); % sorted SI and indices of most influential variables
    N_var_dec = find(cumsum(SI_sorted) > p.Results.DecompositionLimit, 1 ); % how many variables for decomposition
    var_order (N_var_dec+1:end) = 0;
else % modification of the variables selection array into a convenient form
    var_order = zeros(1, N_var);
    for i = 1:N_var 
            if p.Results.OrderOfVariables(i) > 0
                var_order(p.Results.OrderOfVariables(i)) = i;
            end
    end
    N_var_dec = nnz(var_order);
end

% Sorting var_names
var_names_dec = {};
    for f = 1 : N_var_dec
        var_names_dec(f) = p.Results.InputNames(var_order(f));
    end


%% 2. States formation
 
    if isempty(p.Results.NumberOfStates) && isempty(p.Results.StateBoundaries) % automatic definition of states
        states = zeros(1,N_var);
        
        if N_var_dec < 3      % assigning default number of states
            states(1:N_var) = 3;
        else
            states(1:N_var) = 2;
        end
        
        % checking for categorical

        for f = 1 : N_var
            n_unique = numel(unique(inputs(:,f)));
                if n_unique <= 5
                    states(f) = n_unique; % number of states == number of categories
                end
        end

    elseif isempty(p.Results.NumberOfStates) % boundaries are given but number of states is not specified
        states=sum(~isnan(p.Results.StateBoundaries(2:end,:)),1);
    else
        states = p.Results.NumberOfStates;
    end

    % zeroing out not participating in decomposition 
    for f = 1 : N_var
        if ismember(f,var_order(var_order~=0))
        else
            states(f) = 0;
        end
    end
    
    


%% 3. Numeric boundaries
N_runs = size(output,1);

if isempty(p.Results.StateBoundaries) % if automatic definition of boundaries
    boundaries = NaN (max(states),N_var); 
    if matches (p.Results.BoundaryType , 'percentile-based' )
            for f = 1 : N_var
                if states(f) == 0 % skip variables that are not for decomposition
                    f = f + 1;
                else
                    x = inputs(:,f); 
                    if numel(unique(x)) <= 5
                        uniq = unique(x);
                        margin = min(uniq(2:end) - uniq(1:end-1))*0.1; % finding a margin that wouldn't overlap with other states
                            for s = 1 : states(f)-1
                            boundaries(s,f) = uniq(s) + margin;
                            end
                    else
                        x_sorted = zeros (size(x,1)+1,1);
                        x_sorted(1:end-1) = sort (x); x_sorted(end)=x_sorted(end-1)+1;
                        state_size = round(N_runs/states(f));
                        for s = 1 : states(f)-1
                            boundaries(s,f) = x_sorted (state_size*s+1);
                        end
                    end
                    boundaries(states(f),f) = max(x) + 1;
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
                        boundaries(s,f) = f_min + step * s;
                    end
                    boundaries(s,f) = max(inputs(:,f)) + 1; % shifting max of the last state upward to include all datapoints
                end
            end
    end
else
    boundaries = p.Results.StateBoundaries(2:end,:);
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

            if inputs(i,ind_inf) <  boundaries (s,ind_inf)
                states_matching(i,v) = s;
                break
            end 
        end
    end

    [q, scenario_matching(i)] = ismember(states_matching(i,:),Scen_matrix(:,2:end),'rows');
end


scenarios = scenario_matching;

% Adding min values to boundaries for export

    if isempty(p.Results.StateBoundaries) % if automatic definition of boundaries
        
        boundaries_out = NaN(max(states)+1,N_var);
            
        for f = 1 : N_var_dec
            boundaries_out(1,var_order(f)) = min(inputs(:,var_order(f)));
            th = boundaries(states(var_order(f)),var_order(f));
            boundaries(states(var_order(f)),var_order(f)) = th - 1;
        end
       
        boundaries_out(2:end,:) = boundaries;
    else
        boundaries_out = boundaries;
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

%% 8. Visualization

% separating result by scenario
result_dec=cell(1,max(scenarios));

for i=1:max(scenarios)
    s=size(output(scenarios==i),1);
    result_dec(i)=mat2cell(output(scenarios==i),s,1);
end

% defining edges of bins
l=min(output);
h=max(output);
bins=100;
edges=[l:((h-l)/bins):h];

% define frequency of each scenario NPV for each bin
f=zeros(max(scenarios),100);
for i=1:max(scenarios)
    f(i,:)=histcounts(cell2mat(result_dec(i)),edges);
end

% building stacked histogram
figure
stacked_histogram=bar(1:bins,permute(f,[2 1]),1,'stack');


% axis titles
xlabel(p.Results.OutputName);
ylabel('Probability');

% colours
    
    N_main_colors = max(scen_legend(:,2));
    N_shades = scen_legend(end,1) / N_main_colors;
    color = NaN(scen_legend(end,1),3);

        % shades
        sc = 1;
        if N_shades >= 3
            range = 0.8; 
        else
            range = 0.5;
        end
        step = range * 2 / (N_shades - 1);
        beta = [-range : step : range];
        
        for m = 1 : N_main_colors
            current_main_color = hex2rgb(p.Results.MainColors(m));
            
            for sh = 1 : N_shades
                color(sc,:) = brighten (current_main_color, beta(sh));
                sc = sc+1;
            end
        end

        %filling

        for i=1:max(scenarios)
            stacked_histogram(i).FaceColor = color(i,:);
        end
        
        
        for i=1:max(scenarios)
            stacked_histogram(i).EdgeColor = [0.5 0.5 0.5];
        end
        hold on

% ticks
    % real axis ticks
    figure
    h=histogram(output,bins);
    xticks=get(gca,'xtick');

    % to fix matlab R2022b update issue, when the first/last X tick is not
    %   necesseraly at the very edge
        xtick_step = xticks(2)-xticks(1);

        if max(h.BinEdges) > max(xticks)
            xticks = [xticks, xticks(end) + xtick_step]; 
        end

        if min(h.BinEdges) < min(xticks)
            xticks = [xticks(1) - xtick_step, xticks]; 
        end

yticks=get(gca,'ytick');
close
% X axis

% getting corresponding min and max on the artificial axis
min_art=(xticks(1)-min(output))*bins/(max(output)-min(output));
max_art=(xticks(end)-max(output))*bins/(max(output)-min(output))+bins;
distance_art=max_art-min_art;

% corresponding ticks on the artificial axis
xticks_art=(xticks-xticks(1)).*distance_art./(xticks(end)-xticks(1))+min_art;

% Y axis
total=size(output,1);
yticks_art=round(yticks/total,3);
a=[cellstr(num2str(yticks_art'*100))]; % converting values into percentage
pct = char(ones(size(a,1),1)*'%'); % creating vector of % signs
new_yticks = [char(a),pct]; % add the '%' signs after the percentage values

% getting new labels to the graph
set(gca,'XLim',[min(xticks_art), max(xticks_art)],'XTick',xticks_art,'XTickLabel',xticks,'YLim',[0 yticks(end)],'YTickLabel',new_yticks);
hold off


%% showing the legend

    % rename states
    
var1 = {'Color'};
var2 = var_names_dec;
var3 = {'Ymin','Ymean','Ymax', 'Probability'};
var = [var1  var2  var3];


t = array2table(scen_legend,'VariableNames',var);
fig = uifigure;
uit = uitable(fig,'Data', t);

    for i = 1 : scen_legend(end,1)
        s = uistyle('BackgroundColor', color(i,:));
        addStyle(uit,s,'cell',[i,1]);
    end


end