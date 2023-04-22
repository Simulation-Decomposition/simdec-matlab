
%% GET DATA
             
    [Matrix, text_data] = xlsread ("example_data.xlsx");
    
    output = Matrix(:,1);
    inputs = Matrix(:,2:end); % reads only numeric variables           
    [inputs, cats] = append_with_cat(output,inputs,text_data(2:end,:)); % if you have categorical
            % variables, append inputs with string data (does not harm if no categorical variables)



%% RUN AUTOMATIC SIMDEC (the most influential variables are used for decomposition)

    % Calculate significance
    [SI, FOE, SOE]  = significance (output, inputs)
    sum(SI) % shows what portion of the variance of the output is explained by the combined first- and second-order effects of the inputs


    % Initialize decomposition
    dec_limit = 0.8; % cummulative significance threshold (used to decide how many variables to take for decomposition)
    
    output_name = text_data(1,1);
    input_names = text_data(1,2:end);

    threshold_type = 2; % 1 for 'percentile-based' (same amount of observations in each state),
                        % 2 for 'median-based' (equaly-spaced ranges)

    [scenarios, scen_legend, thresholds, var_names_dec] = decomposition (output, inputs, SI, dec_limit, [], [], [], threshold_type, input_names);
    
    % Build the graph
    sm = build_simdec_chart(output,scenarios,scen_legend,[],output_name,var_names_dec);


%% CREATE CUSTOM DECOMPOSITION (choose inputs & alter states) 

    % Provide specifications  
    manual_vars = [0 2 1 0]; % specify the order of variables for decomposition, use 0 to exclude, size (1, N_inputs)
                             % in this example we set that the third input variable is used first, and then the second one.  
    manual_states = [0 3 2 0]; % specify the number of states for each variable, size (1, N_inputs), the position corresponds to the original order of inputs
    manual_thresholds = [ NaN    min(inputs(:,2))     min(inputs(:,3))    NaN
                          NaN         100                  657.5          NaN
                          NaN         650             max(inputs(:,3))    NaN      
                          NaN    max(inputs(:,2))          NaN            NaN]; % specify numeric thresholds for every state, size (max(manual_states)+1, N_inputs).


    [scenarios, scen_legend, thresholds, var_names_dec] = decomposition (output, inputs, SI, [], manual_vars, manual_states, manual_thresholds, [], input_names);
        

    % Build the graph
    main_colors = {'#8c5eff'; '#ffe252'; '#0dd189'}; % Number of main colors should be equal to the number of states of the first for decomposition variable. 
    sm = build_simdec_chart(output,scenarios,scen_legend,main_colors,output_name,var_names_dec);
