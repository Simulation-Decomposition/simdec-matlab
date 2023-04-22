function [sortedX, indexMap] = cat_transform(Y, X)
% Transforms categorical variable into numeric, by assigning indices in
% accordance with averag Y values in each category.
%
% * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
% INPUTS
%   Y               - a 1D array of Y values
%   X               - a 1D array of X values
% 
% OUTPUT
%   sortedX         - categories of X sorted by increasing averag Y
%   indexMap        - an array of the same size as X with assigned numeric
%                     indices correspodning to original categorical values
%
% * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
% Written by Mariia Kozlova, last updated 16.4.2023
% Many thanks for the grant #220178 from Finnish Foundation for Economic
% Education (lsr.fi) and the grant #6713/31/2021 from Business Finland.
% * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *


% Find all unique values of X
uniqueX = unique(X);

% Compute average of Ys for each unique value of X
averages = zeros(size(uniqueX));
for i = 1:length(uniqueX)
    idx = strcmp(X, uniqueX(i));
    averages(i) = mean(Y(idx));
end

% Sort the unique values of X in accordance by average Y
[~, sortedIndices] = sort(averages);
sortedX = uniqueX(sortedIndices);

% Assign indices to sorted unique values of X
indexMap = zeros(size(X));
for i = 1:length(sortedX)
    idx = strcmp(X, sortedX(i));
    indexMap(idx) = i;
end


end