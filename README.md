> **Warning**
> This library is under active development and things can change at anytime! Suggestions and help are greatly appreciated.

![image](https://user-images.githubusercontent.com/37065157/233836694-5312496e-4ada-47cb-bc09-3bf8c00be135.png)

**Simulation decomposition** or **SimDec** is an uncertainty and sensitivity analysis method, which is based on Monte Carlo simulation. 
SimDec consists of two major parts (and two corresponding funciton):
1. Computing sensitivity indices, 
2. Visualizing the effects by creating multi-variable scenarios, mapping the output values to them, and color-coding the corresponding segments of the output distribution. 

SimDec reveals the nature of causalities and interaction effects in the model.  



## Starting to work
Simply download the matlab [functions](functions) and add them to your path when working with matlab. 

Use [sensitivity_indices.m](sensitivity_indices.m) for computing the indices and [simdec_visualization.m](simdec_visualization.m) for building the graphics.

## Example
The following procedure is saved in the [main.m](main.m) and uses [example_data.xlsx](example_data.xlsx).

### 1. Load data 
First the simulated `inputs` and the `output` need to be specified. They can result from a Monte Carlo simulation arranged directly in matlab, or conducted elsewhere and then loaded through a file, like in this example. 

```matlab   
    Matrix = xlsread ("example_data.xlsx");
    
    output = Matrix(:,1);
    inputs = Matrix(:,2:end);         
```


### 2. Compute sensitivity indices
Function [sensitivity_indices.m](sensitivity_indices.m) computes first-order effects `FOE` (main individual effect of every input variable), second-order effects `SOE` (interaction effects between pairs of variables) and combined sensitivity indices `SI`. 

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

Each value shows what portion of the variance of the output is explained (negative SOE values indicate correlation). In this example, SI shows that the most important inputs are X2 (52%) and X4 (35%). SOE points out that there is interaction between X2 and X3 (11%) and correlation between X2 and X3 (-6%). In total, 100% of the output variance is explained (a few percentage points can be attributed to noise).


### 3. Visualize
Function [simdec_visualization.m](simdec_visualization.m) 
1. chooses the most important input variables, 
2. breaks them down into states, 
3. forms scenarios out of all combinations of those states, 
4. maps the scenarios onto the output values, and 
5. visualizes these scenarios by color-coding the distribution of the output.

```matlab
[scenarios, scen_legend, boundaries] = simdec_visualization (output, inputs, SI);
```


The SimDec graph and the corresponding legend is generated automatically when running this function. 

![image](https://github.com/Simulation-Decomposition/simdec-matlab/assets/37065157/2304f44a-05b1-4c9d-88c8-d4862ca54258)


Feel free to go an extra step, - name the states (i.e., *low*, *medium*, *high*) and merge the cells of the legend with the same state. 

That's it, you SimDec analysis is completed. Unless you want to customize it furhter.

### 4. Customize

The [simdec_visualization.m](simdec_visualization.m) function has numerious optional arguments that can be used to polish the outlook of the results, tune and play with the decomposition set-up.

#### 4.1. Polishing: colors, names of variables

Use `'OutputName'`, `'InputNames'` to add names of the variables. 
Specify `'MainColors'` as a cell array of HEX codes of the desired main colors (the shades for subdivisions are produced automatically from the main colors). 

```matlab
    output_name = 'Output';
    input_names = {'Input1','Input2','Input3','Input4'};
    colors = {'#3F45D0','#DC267F','26DCD1'}; % HEX codes for the main colors
    
    [scenarios, scen_legend, boundaries] = simdec_visualization (output, inputs,...
        SI,'OutputName',output_name,'InputNames',input_names,'MainColors',colors);
```

![image](https://github.com/Simulation-Decomposition/simdec-matlab/assets/37065157/0d9be2e4-7f8f-462d-a3cc-2e34a19f0ea4)



#### 4.2. Tuning decomposition: more variables, different state formation

Deafult number of variables for decomposition is defined beased on the threshold `0.8*sum(SI)`. The threshold can be changed by using `'DecompositionLimit'` argument.

The `'BoundaryType'` argument defines how the numeric range of input variables is broken down into states. The deafult value `'percentile-based'` forms the states ensuring the same amount of observations in each state. The alternative `'median-based'` approach divides input ranges into equaly-spaced intervals. Changing this argument would not make any difference for independent uniformly distributed inputs.


```matlab
    dec_limit = 0.9; % minimum overall importance [sum(SI)] of chosen for decomposition input variables
    boundary_type = 'median-based'; 
                                   
    [scenarios, scen_legend, boundaries] = simdec_visualization (output, inputs,...
        SI,'DecompositionLimit',dec_limit,'BoundaryType',boundary_type);
```

![image](https://github.com/Simulation-Decomposition/simdec-matlab/assets/37065157/305b660f-df12-4a78-bc34-aa61d004d760)



#### 4.3. Tuning decomposition: more variables, different state formation
The decomposition can be further fully customized by altering the amount and the order of inputs for decomposition(`'OrderOfVariables'`), defining custom amount of states (`'NumberOfStates'`) and their boundaries (`'StateBoundaries'`). 

```matlab
    manual_vars = [0 2 1 0]; % specify the order of variables for decomposition, use 0 to exclude, size (1, N_inputs). In this example we set that the third input variable is used first, and then the second one.  
    manual_states = [0 3 2 0]; % specify the number of states for each variable, size (1, N_inputs), the position corresponds to the original order of inputs. Three states for the second input variable and two states for the third. 
    manual_boundaries =  [ NaN    min(inputs(:,2))     min(inputs(:,3))    NaN
                           NaN         100                  657.5          NaN
                           NaN         650             max(inputs(:,3))    NaN      
                           NaN    max(inputs(:,2))          NaN            NaN]; % specify numeric boundaries for every state, size(max(manual_states)+1, N_inputs).


    [scenarios, scen_legend, boundaries] = simdec_visualization (output, inputs,...
        SI,'OrderOfVariables',manual_vars,'NumberOfStates',manual_states,'StateBoundaries',manual_boundaries);
```

![image](https://github.com/Simulation-Decomposition/simdec-matlab/assets/37065157/ea3b79e1-c969-467d-a817-b23c55a01402)

The optional function attributes can be used in any combination. 


## Links
- See [how to read SimDec on wikipedia](https://en.wikipedia.org/wiki/SimDec)
- [People behind SimDec](https://www.simdec.fi/team)
- Join our [Sensitivity Analysis discord community](https://discord.gg/54SFcNsZS4)


## Citations

The algorithms and visualizations used in this package came primarily out of research at **LUT University**, Lappeenranta, Finland, and **Stanford University**, California, U.S., supported with grants from ***Business Finland***, ***Wihuri Foundation***, ***Foundation for Economic Education***, and ***Natural Sciences and
Engineering Research Council***. If you use SimDec in your research we would appreciate a citation to the following publications:

- Kozlova, M., & Yeomans, J. S. (2022). Monte Carlo Enhancement via Simulation Decomposition: A “Must-Have” Inclusion for Many Disciplines. _INFORMS Transactions on Education, 22_(3), 147-159. [Available here](https://pubsonline.informs.org/doi/10.1287/ited.2019.0240).
- Kozlova, M., Moss, R. J., Yeomans, J. S., & Caers, J. (forthcoming). Uncovering Heterogeneous Effects in Computational Models for Sustainable Decision-making. _Environmental Modelling & Software_. [Available here](https://papers.ssrn.com/sol3/papers.cfm?abstract_id=4550911)
- Kozlova, M., Moss, R. J., Roy, P., Alam, A., & Yeomans, J. S. (forthcoming). SimDec algorithm. In M. Kozlova & J. S. Yeomans (Eds.), _Sensitivity Analysis for Business, Technology, and Policymaking. Made Easy with Simulation Decomposition_. Routledge.