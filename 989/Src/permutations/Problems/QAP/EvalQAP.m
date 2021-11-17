function[val] =  EvalQAP(vector)
% [val] =  EvalQAP(vector)
% EvalQAP: Using two nxn matrix computes the cost of a given permutation
% Input:
% vector: Permutation vector
% val:  Cost of the permutation according to global QAPInstance processingtimes variable
%
% NOTE: In order to use the function correctly you must declare a global
% boolean variable depending on these cases:
%       - true: The permutation will be inversed before the evaluation.
%       - false: The permutation will NOT be inversed.
%
%   By default the value is true
%
% Created version 03/03/2014. Josian Santamaria (jasantamaria003@ikasle.ehu.es) 
%
% Last version 03/19/2014. Josian Santamaria (jasantamaria003@ikasle.ehu.es) 

global QAPInstance;
global inverse;

%if inverse is not defined by the user, put the default value.
if isempty(inverse)
    inverse=true;
end

PermutatioSize = size(vector,2);
val = 0;

[distance, flow, n] = QAPInstance{:}; % Get the variables from the processed instance

if inverse
	genes=Invert(vector);
else
	genes=vector;
end

% index variable
i=1:PermutatioSize;

% For use in matlab, we use less variables
m = flow(genes,genes).*distance(i,i); %in the papper m = flow(i,j).*distance(genes(i),genes(j));


val = -sum(m(:)); % sum of the matrix, negative in order to maximize the fitness