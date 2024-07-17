
%% DATA
             
    Matrix = xlsread ("example_data.xlsx");
    
    output = Matrix(:,1);
    inputs = Matrix(:,2:end);      

%% AUTOMATIC SIMDEC
    
    [SI, FOE, SOE]  = sensitivity_indices (output, inputs) % Calculate sensitivity indices
    sum(SI) % shows what portion of the variance of the output is explained 
            % by the combined first- and second-order effects of the inputs
            
    [scenarios, scen_legend, boundaries] = simdec_visualization (output, inputs, SI); % Build visualization
            

%% BOXPLOT
    
    [scenarios, scen_legend, boundaries] = simdec_visualization (output, inputs, SI,'GraphType','boxplot'); % Boxplot alternative visualization


%% POLISHING: colors, names of variables

    output_name = 'Output';
    input_names = {'Input1','Input2','Input3','Input4'};
    colors = {'#3F45D0','#DC267F','26DCD1'}; % HEX codes for the main colors
    
    [scenarios, scen_legend, boundaries] = simdec_visualization (output, inputs, SI,'OutputName',output_name,'InputNames',input_names,'MainColors',colors);
    

%% TUNING: more variables, different state formation
    
    dec_limit = 0.9; % minimum overall importance [sum(SI)] of chosen for decomposition input variables
    boundary_type = 'interval-based'; % divides input ranges into states using equaly-spaced intervals rather than default same amount of observations in each state 
                                    % This would not make any difference for independent uniformly distributed inputs
                                   
    [scenarios, scen_legend, boundaries] = simdec_visualization (output, inputs, SI,'DecompositionLimit',dec_limit,'BoundaryType',boundary_type);


%% CUSTOMIZATION: own order of variables, own states


    manual_vars = [0 2 1 0]; % specify the order of variables for decomposition, use 0 to exclude, size (1, N_inputs)
                             % in this example we set that the third input variable is used first, and then the second one.  
    manual_states = [0 3 2 0]; % specify the number of states for each variable, size (1, N_inputs), the position corresponds to the original order of inputs
    manual_boundaries =  [ NaN    min(inputs(:,2))     min(inputs(:,3))    NaN
                           NaN         100                  657.5          NaN
                           NaN         650             max(inputs(:,3))    NaN      
                           NaN    max(inputs(:,2))          NaN            NaN]; % specify numeric boundaries for every state, size(max(manual_states)+1, N_inputs).


    [scenarios, scen_legend, boundaries] = simdec_visualization (output, inputs, SI,'OrderOfVariables',manual_vars,'NumberOfStates',manual_states,'StateBoundaries',manual_boundaries);

   
%% Two-output visualization 

    Matrix = xlsread ("example_data2.xlsx");
    
    output = Matrix(:,1);
    output2 = Matrix(:,2);
    inputs = Matrix(:,3:end); 
    
    [SI, FOE, SOE]  = sensitivity_indices (output, inputs) 
    
    outputname = 'Output1';
    output2name = 'Output2';

    share_of_data_shown = 0.005;
    n_bins = 40;

    [scenarios, scen_legend, boundaries] = simdec_visualization (output, inputs, SI,'OutputName',outputname...
        ,'Output2',output2,'Output2Name',output2name,'ScatterFraction',share_of_data_shown,'NumberOfBins',n_bins);

    % Changing XY limits for the first histogram (scatter plot and the second histogram scale accordingly)
    xlim_values = [1000 3000]; % affects x-axis of the top histogram and x-axis of the scatterplot.
    ylim_values = [0 4]; % value in % for Probability y-axis of top histogram, does not affect the scatterplot.
    [scenarios, scen_legend, boundaries] = simdec_visualization (output, inputs, SI,'OutputName',outputname...
        ,'Output2',output2,'Output2Name',output2name,'ScatterFraction',share_of_data_shown,'NumberOfBins',n_bins...
        ,'XLim',xlim_values,'YLim',ylim_values);
    
    % Enforcing a desired combination of XY limits for the scatterplot
    xlim_values = [1000 3000]; % affects x-axis of the top histogram and x-axis of the scatterplot.
    xlim_values_2 = [0 1000]; % affects (rotated) x-axis of the right histogram and y-axis of the scatterplot.
    [scenarios, scen_legend, boundaries] = simdec_visualization (output, inputs, SI,'OutputName',outputname...
        ,'Output2',output2,'Output2Name',output2name,'ScatterFraction',share_of_data_shown,'NumberOfBins',n_bins...
        ,'XLim',xlim_values,'XLim2',xlim_values_2);