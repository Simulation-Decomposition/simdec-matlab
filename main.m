
%% DATA
             
    Matrix = xlsread ("example_data.xlsx");
    
    output = Matrix(:,1);
    inputs = Matrix(:,2:end);      

%% AUTOMATIC SIMDEC
    
    [SI, FOE, SOE]  = sensitivity_indices (output, inputs) % Calculate sensitivity indices

    [scenarios, scen_legend, boundaries] = simdec_visualization (output, inputs, SI); % Build visualization


%% POLISHING: colors, names of variables

    output_name = 'Output';
    input_names = {'Input1','Input2','Input3','Input4'};
    colors = {'#3F45D0','#DC267F','26DCD1'}; % HEX codes for the main colors
    
    [scenarios, scen_legend, boundaries] = simdec_visualization (output, inputs, SI,'OutputName',output_name,'InputNames',input_names,'MainColors',colors);
    

%% TUNING: more variables, different state formation
    
    dec_limit = 0.9; % minimum overall importance [sum(SI)] of chosen for decomposition input variables
    boundary_type = 'median-based'; % divides input ranges into states using equaly-spaced intervals rather than default same amount of observations in each state 
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
