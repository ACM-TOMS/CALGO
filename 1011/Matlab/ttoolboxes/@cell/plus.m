function C = plus(A,B)
% C = plus(A,B)
% Computes A+B for cells
% Input: 
%   A and B are cell arrays                     Computes A plus B, Same Behaviour as minus for arrays
%   A or B is a cell array, the other is not    Compute Cell plus Element for all elements in the cell
%
% E.g.: num2cell(randn(3))+num2cell(randn(3))
%       {1 2;3 4;5 6}+[1 2]
%       {1 2;3 4;5 6}+3
% 
% See also: cell/minus, plus
%
% Written by: tommsch, 2018

if(iscell(A) && ~iscell(B))
    C=cellfun(@(x) x+B,A,'UniformOutput',false);  %# Apply plus cell-wise
elseif(~iscell(A) && iscell(B))
    C=cellfun(@(x) A+x,B,'UniformOutput',false);  %# Apply plus cell-wise
elseif(iscell(A) && iscell(B))
  C = cellfun(@plus,A,B,'UniformOutput',false);
end
end