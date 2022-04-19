function C = uminus(A)
% C = uminus(A)
% Computes -A for cells
% Input: 
%   A           cell arrays
%
% E.g.: num2cell(randn(3))-num2cell(randn(3))
%       {1 2;3 4}-{5 7;9 11}
%       {1 2;3 4;5 6}-[7 9]
% 
% See also: cell/plus, minus
%
% Written by: tommsch, 2018

C=cellfun(@(x) -x,A,'UniformOutput',false);  %# Apply plus cell-wise

end