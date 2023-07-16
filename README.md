> **Warning**
> This library is under active development and things can change at anytime! Suggestions and help are greatly appreciated.

![image](https://user-images.githubusercontent.com/37065157/233836694-5312496e-4ada-47cb-bc09-3bf8c00be135.png)

**Simulation decomposition** or **SimDec** is an uncertainty and sensitivity analysis method, which is based on Monte Carlo simulation. SimDec consists of three major parts:
1. computing significance indices, 
2. creating multi-variable scenarios and mapping the output values to them, and 
3. visualizing the scenarios on the output distribution by color-coding its segments. 

SimDec reveals the nature of causalities and interaction effects in the model.  
See our [publications](https://www.simdec.fi/publications) and join our [discord community](https://discord.gg/54SFcNsZS4).



## Starting to work
Simply download the matlab [functions](functions) and add them to your path when working with matlab.

## Example
The following procedure is saved in the [main.m](main.m) and uses [example_data.xlsx](example_data.xlsx).

### Load data 
First the simulated `inputs` and the `output` need to be specified. They can result from a Monte Carlo simulation arranged directly in matlab, or conducted elsewhere and then loaded through a file, like in this example. The `append_with_cat` function reads any string variables and converts them into numeric ones for further processing. 

```matlab
%% GET DATA
             
    [Matrix, text_data] = xlsread ("example_data.xlsx");
    
    output = Matrix(:,1);
    inputs = Matrix(:,2:end); % reads only numeric variables           
    [inputs, cats] = append_with_cat(output,inputs,text_data(2:end,:)); % if you have categorical
            % variables, append inputs with string data (does not harm if no categorical variables)
```


### Compute significance indices
Function `significance` computes first-order effects `FOE` (main individual effect of every input variable), second-order effects `SOE` (interaction effects between pairs of variables and combined sensitivity indices `SI`. 

```matlab
    [SI, FOE, SOE]  = significance (output, inputs)
    sum(SI) % shows what portion of the variance of the output is explained 
            % by the combined first- and second-order effects of the inputs
```
Here is the result it returns:

SI =

    0.0409
    0.5155
    0.0955
    0.3506

FOE =

    0.0367
    0.4910
    0.1069
    0.2777

SOE =

         0    0.0034    0.0015    0.0035
         0         0   -0.0605    0.1059
         0         0         0    0.0363
         0         0         0         0

sum(SI) = 

    1.0024

Each value shows what portion of the variance of the output is explained (negative SOE values indicate correlation). In this example, SI shows that the most significant inputs are X2 (52%) and X4 (35%). SOE points out that there is interaction between X2 and X3 (11%) and correlation between X2 and X3 (-6%). In total, 101% of the output variance is explained (1% can be attributed to noise).


### Run decomposition
Function `decomposition` chooses the most important input variables, breaks them down into states, forms scenarios out of all combinations of those states and maps the scenarios onto the output values.

```matlab
    % Initialize decomposition
    dec_limit = 0.8; % cummulative significance threshold 
                     % (used to decide how many variables to take for decomposition)
    
    output_name = text_data(1,1);
    input_names = text_data(1,2:end);

    threshold_type = 2; % 1 for 'percentile-based' (same amount of observations in each state),
                        % 2 for 'median-based' (equaly-spaced ranges)

    [scenarios, scen_legend, thresholds, var_names_dec] = ...
            decomposition (output, inputs, SI, dec_limit, [], [], [], threshold_type, input_names);
```


### Visualize
The SimDec graph and the corresponding legend is created with the function `build_simdec_chart`. 

```matlab
    % Build the graph
    sm = build_simdec_chart(output,scenarios,scen_legend,[],output_name,var_names_dec);
```
![image](https://user-images.githubusercontent.com/37065157/233839183-7852735f-4c80-47da-9f26-d4ac8d962dd3.png)

Feel free to go an extra step, - name the states and merge the cells of the legend with the same state. 


### Customize
There are a number of ways to customize the visuals. One can choose different input variables for decomposition, predefine the number of states and specific numeric threshold, and most importantly, change the colors. Here is an example of all of those.

```matlab
%% CREATE CUSTOM DECOMPOSITION (choose inputs & alter states) 

    % Provide specifications  
    manual_vars = [0 2 1 0]; % specify the order of variables for decomposition, use 0 to exclude, size (1, N_inputs)
                             % in this example we set that the third input variable is used first, and then the second one.  
    manual_states = [0 3 2 0]; % specify the number of states for each variable, size (1, N_inputs), the position corresponds 
                               % to the original order of inputs
    manual_thresholds = [ NaN    min(inputs(:,2))     min(inputs(:,3))    NaN
                          NaN         100                  657.5          NaN
                          NaN         650             max(inputs(:,3))    NaN      
                          NaN    max(inputs(:,2))          NaN            NaN]; % specify numeric thresholds for every state, 
                                                                                % size (max(manual_states)+1, N_inputs).


    [scenarios, scen_legend, thresholds, var_names_dec] = ...
             decomposition (output, inputs, SI, [], manual_vars, manual_states, manual_thresholds, [], input_names);
        

    % Build the graph
    main_colors = {'#8c5eff'; '#ffe252'; '#0dd189'}; % Number of main colors should be equal to the number 
                                                     % of states of the first for decomposition variable. 
    sm = build_simdec_chart(output,scenarios,scen_legend,main_colors,output_name,var_names_dec);
```
![image](https://user-images.githubusercontent.com/37065157/233839368-f3cc430b-f974-4de2-9d80-27179eaa1458.png)


## Code structure
Each block in Figure below is a matlab function. The green ones are higher-level functions that are called in the main script (i.e. example above).


![scheme](https://user-images.githubusercontent.com/37065157/234074889-719ea46b-f542-4ef5-8709-542747fc17c1.png)



## Citations

The algorithms and visualizations used in this package came primarily out of research at LUT University, Lappeenranta, Finland, and Stanford University, California, U.S., supported with grants from Business Finland, Wihuri Foundation, and Finnish Foundation for Economic Education. If you use SimDec in your research we would appreciate a citation to the following publications:

- Kozlova, M., & Yeomans, J. S. (2022). Monte Carlo Enhancement via Simulation Decomposition: A “Must-Have” Inclusion for Many Disciplines. _INFORMS Transactions on Education, 22_(3), 147-159. [Available here](https://pubsonline.informs.org/doi/10.1287/ited.2019.0240).
- Kozlova, M., Moss, R. J., Yeomans, J. S., & Caers, J. (forthcoming). Uncovering Heterogeneous Effects in Computational Models for Sustainable Decision-making. _Environmental Modelling & Software_. 
- Kozlova, M., Moss, R. J., Roy, P., Alam, A., & Yeomans, J. S. (forthcoming). SimDec algorithm. In M. Kozlova & J. S. Yeomans (Eds.), _Sensitivity Analysis for Business, Technology, and Policymaking Made Easy with Simulation Decomposition_. Routledge.
