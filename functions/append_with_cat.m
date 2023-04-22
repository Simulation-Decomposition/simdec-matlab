function [all_inputs, sorted_X] = append_with_cat(output,numeric_inputs,string_data)
% Appends numeric inputs with categorical variables which are transformed
% into numeric ones. 
%
% Uses function cat_transform.
%
% * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
% INPUTS
%   output            - a 1D array of Y values
%   numeric_inputs    - a 2D array of X values, all inputs that a numeric
%   string_data       - all non-numeric values in the dataset
% 
% OUTPUT
%   all_inputs      - a 2D array containing all numeric and all transformed
%                     input variables
%   sorted_X        - the order of categories in increasing averag Y for
%                     each non-numeric input variables
%
% * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
% Written by Mariia Kozlova, last updated 16.4.2023
% Many thanks for the grant #220178 from Finnish Foundation for Economic
% Education (lsr.fi) and the grant #6713/31/2021 from Business Finland.
% * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *


% get rid of empty columns
string_data(:, all(cellfun(@isempty, string_data), 1)) = [];

% get number of categorical variables
n_cats = size(string_data,2); 

% initialize matrix for transofrmed categorical variables & cell array for sorted categories
X_cat_num = NaN(length(output),n_cats); 
sorted_X = cell(1, n_cats);

% check Y averages to wisely index the categories
        for c = 1 : n_cats
            [sorted_X{c}, X_cat_num(:,c)] = cat_transform(output, string_data(:,c));
        end
        
% append inputs 
all_inputs = [numeric_inputs X_cat_num];


end