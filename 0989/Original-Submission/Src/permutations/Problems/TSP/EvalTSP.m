function[val] =  EvalTSP(vector)
% [val] =  EvalTSP(vector)
% EvalTSP: Using a distance matrix computes the cost of a given permutation
% Input:
% vector: Permutation vector
% val:  Cost of the permutation according to global TSPInstance distance matrix variable
%
% NOTE: In order to use the function correctly you must declare a global
% boolean variable depending on these cases:
%       - true: The permutation will be inversed before the evaluation.
%       - false: The permutation will NOT be inversed.
%
%   By default the value is true
%
% Created version 12/15/2013. Roberto Santana (roberto.santana@ehu.es)  
%
% Last version 05/13/2014. Josian Santamaria (jasantamaria003@ikasle.ehu.es)

global TSPInstance;
global inverse;

%if inverse is not defined by the user, put the default value.
if isempty(inverse)
    inverse=true;
end

PermutatioSize = size(vector,2);
val = 0;

[distance, n] = TSPInstance{:}; % Get the variables from the processed instance

if inverse
	genes=Invert(vector);
else
	genes=vector;
end

%Calculate the cost of the given permutation, visiting all the cities
%returning to the origin city.
for i=1:PermutatioSize-1 %Visit from the first city to the last city
    val = val + distance(genes(i),genes(i+1));
end
%Return to the origin city
val = val + distance(genes(PermutatioSize),genes(1));

val = - val;

