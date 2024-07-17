function [scenarios, scen_legend, boundaries_out] = simdec_visualization(output, inputs, SI, varargin)
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
%                        'percentile-based' for the same amount of observations 
%                        in each state,   
%                        'interval-based' for equally-spaced ranges of states.
%   NumberOfBins         - number of bins for histogram
%   XLim                - Minimum and maximum values for x axis
%                         [xmin xmax].
%   YLim                - Minimum and maximum values for y axis [ymin ymax]. 
%                         Fox boxplot, YLim is ignored
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
%   GraphType           - stacked histogram by default, boxplot as an alterantive
%   Output2             - a second output variable, that will be displayed
%                         with a scatterhit
%   Output2Name         - name for the second output variable
%   ScatterFraction     - the portion of data / points displayed on the
%                         scatterplot. The default value is 1 - the entire 
%                         dataset is displayed. A value of e.g. 0.5 will show
%                         every second point. 
%   XLim2                - Minimum and maximum values for x axis of the
%                          second histogram (output2 values)[xmin xmax].
%                          This scales the scatteplot (ylim) accordingly.
%   YLim2                - Minimum and maximum values for y axis of the second histogram 
%                          (% values of Probability) [ymin ymax]. 
%
% 
% OUTPUTS
%
%   scenarios          - an array of the same size as Y with scenario indices 
%                      for every simulation run
%   scen_legend        - a scenario table that shows which states of which
%                      variables compose different scenarios
%   boundaries_out     - numeric boundaries of states of input variables
%  
% * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
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
    expected_boundary_types = {'percentile-based','interval-based'}; % 'interval-based' for equaly-spaced ranges
    default_colors = {'#DC267F'; '#E8EA2F'; '#26DCD1';'#C552E4';'#3F45D0'};
    defult_output_name = 'Y';
    default_input_names = cell(1, size(inputs,2));
        for i = 1:size(inputs,2)
            default_input_names{i} = ['X' num2str(i)];
        end
    default_graph_type = 'stacked_histogram';
    expected_graph_types = {'stacked_histogram','boxplot'};
    default_output2 = [];
    default_output2name = 'Y2';
    default_scatterfraction = 1;


   p = inputParser;
   validScalarPosNum = @(x) isnumeric(x) && isscalar(x) && (x > 0);
   addRequired(p,'output');
   addRequired(p,'inputs');
   addRequired(p,'SI');
   addOptional(p,'DecompositionLimit',default_decomposition_limit,validScalarPosNum);
   addOptional(p,'OrderOfVariables',default_var_order, @(x) isequal(size(x), expected_size));
   addOptional(p,'NumberOfStates',default_number_states, @(x) isequal(size(x), expected_size));
   addOptional(p,'BoundaryType',default_boundary_type,@(x) any(validatestring(x,expected_boundary_types)));
   addOptional(p,'NumberOfBins',100,validScalarPosNum);
   addOptional(p,'XLim',[]); 
   addOptional(p,'YLim',[]); 
   addOptional(p,'StateBoundaries',default_boundaries);
   addParameter(p,'MainColors',default_colors,@iscellstr);
   addParameter(p,'OutputName',defult_output_name);
   addParameter(p,'InputNames',default_input_names,@iscellstr);
   addOptional(p,'GraphType',default_graph_type,@(x) any(validatestring(x,expected_graph_types)));
   addOptional(p,'Output2',default_output2); 
   addOptional(p,'Output2Name',default_output2name); 
   addOptional(p,'ScatterFraction',default_scatterfraction); 
   addOptional(p,'XLim2',[]); 
   addOptional(p,'YLim2',[]); 

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
    
    else % interval-based thresholds

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
        boundaries_out = p.Results.StateBoundaries;
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


% 8.1 Single output visualization
        
figure

if isempty(p.Results.Output2) 

    if matches (p.Results.GraphType, 'stacked_histogram') % if stacked histogram
        
        % setting the axes limits
        hold on;
        rotate90 = 0; 
        hide_axes = 0;     
        stacked_histogram(output, p.Results.XLim, p.Results.YLim, p.Results.OutputName, scenarios, color, p.Results.NumberOfBins, rotate90, hide_axes)  

    
    else % if boxplot
        ind_0=find(scenarios~=0);
        b = boxplot(output(ind_0), scenarios(ind_0),'Orientation','horizontal','FactorDirection','data','Symbol','ok','Color','k'); % 'FactorDirection','list' to switch the order of boxes
        b = findobj(gca,'Tag','Box');
        for j=1:length(b) % coloring
            patch(get(b(j),'XData'),get(b(j),'YData'),color(length(b)-j+1,:),'FaceAlpha',.8);
        end
    
    
        xlabel(p.Results.OutputName);
        ylabel('Scenario');
    
        % set x-axis
        if ~isempty(p.Results.XLim)
            xlim(p.Results.XLim)
        end    
    end





% 8.2 Two-output visualization

else % if two outputs


    % Scatter plot
    subplot(2,2,3);
    X = output(1:round(1/p.Results.ScatterFraction):end);
    Y = p.Results.Output2(1:round(1/p.Results.ScatterFraction):end);
    scenarios_scatter = scenarios(1:round(1/p.Results.ScatterFraction):end);
    scatter(X, Y, 20, color(scenarios_scatter,:), 'filled');
    xlabel(p.Results.OutputName);
    ylabel(p.Results.Output2Name);
        % set / record axes
        scatter_xlim = get(gca, 'XLim');
        scatter_ylim = get(gca, 'YLim');
        
        if ~isempty(p.Results.XLim)
            scatter_xlim = p.Results.XLim;
            xlim(p.Results.XLim);
            if isempty(p.Results.XLim2) % updating ylim if not specified
                scatter_ylim = get(gca, 'YLim'); 
            end
        end

        if ~isempty(p.Results.XLim2)
            scatter_ylim = p.Results.XLim2;
            ylim(p.Results.XLim2);
            if isempty(p.Results.XLim) % updating xlim if not specified 
                scatter_xlim = get(gca, 'XLim'); 
            end
        end

    hold on;


    % Top histogram for output1
    subplot(2,2,1);
    hold on;
    rotate90 = 0; % the upper histogram is vertical
    hide_axes = 1; % axes are hidden for subplots
    stacked_histogram(output, scatter_xlim, p.Results.YLim, p.Results.OutputName, scenarios, color, p.Results.NumberOfBins, rotate90, hide_axes)     
    hold off;
    
    
    % Right histogram for output2
    subplot(2,2,4);
    hold on;
    rotate90_2 = 1; % rotated 90 degrees
    hide_axes_2 = 1; % axes are hidden for subplots
    stacked_histogram(p.Results.Output2, scatter_ylim, p.Results.YLim2, p.Results.Output2Name, scenarios, color, p.Results.NumberOfBins, rotate90_2, hide_axes_2)     
    hold off;
    
    
    % Bring histograms closer
    s1 = subplot(2,2,3); 
    set(s1, 'Units', 'normalized');
    set(s1, 'Position', [0.1300 0.1100 0.44 0.44]);
    
    h1 = subplot(2,2,1); 
    set(h1, 'Units', 'normalized');
    set(h1, 'Position', [0.1300 0.549 0.44 0.4]);
    
    h2 = subplot(2,2,4); 
    set(h2, 'Units', 'normalized');
    set(h2, 'Position', [0.569 0.1100 0.4 0.44]);



end

%% showing the legend

    % rename states
    
var1 = {'Color'};
var2 = var_names_dec;
var3 = {'Ymin','Ymean','Ymax', 'Probability'};
var = [var1  var2  var3];


t = array2table(scen_legend,'VariableNames',var);
fig = uifigure('HandleVisibility', 'on');
uit = uitable(fig,'Data', t);

    for i = 1 : scen_legend(end,1)
        s = uistyle('BackgroundColor', color(i,:));
        addStyle(uit,s,'cell',[i,1]);
    end


end