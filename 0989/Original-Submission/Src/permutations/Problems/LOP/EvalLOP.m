function[val] =  EvalLOP(vector)
% [val] =  EvalLOP(vector)
% EvalLOP: Using a nxn matrix computes the cost of a given permutation
% Input:
% vector: Permutation vector
% val:  Cost of the permutation according to global LOPInstance processingtimes variable
%
% NOTE: In order to use the function correctly you must declare a global
% boolean variable depending on these cases:
%       - true: The permutation will be inversed before the evaluation.
%       - false: The permutation will NOT be inversed.
%
%   By default the value is true
%
% Created version 02/21/2014. Josian Santamaria (jasantamaria003@ikasle.ehu.es) 
%
% Last version 03/19/2014. Josian Santamaria (jasantamaria003@ikasle.ehu.es)

global LOPInstance;
global inverse;

%if inverse is not defined by the user, put the default value.
if isempty(inverse)
    inverse=true;
end

PermutatioSize = size(vector,2);
val = 0;

[matrix, n] = LOPInstance{:}; % Get the variables from the processed instance

if inverse
	genes=Invert(vector);
else
	genes=vector;
end
C=matrix; % nxn matrix of LOPInstance
m = C(genes,genes);

uppertri = triu(m,1); % superdiagonal
val = sum(uppertri(:)); % sum of the superdiagonal