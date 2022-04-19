function C = mrdivide(A,B)
% C = mrdivide(A,B)
% Elementwise /-operator for cells.
%
% Input: 
%   A and B are cell arrays                     Computes a/b, where a and b are the corresponding cells in A.
%                                               In other words, computes A/B elementwise
%   A is a cell, B is not                       Computes a/B, where a are the cells of A
%   B is a cell, A is not                       Computes A/b, where b are the cells of B
%
% E.g.: vdisp({[1 2; 1 0], [1 1]}\{[1 0;0 1],[1 0]})
%       vdisp({[1 2; 1 0],[1 0; 0 1]}\[1 0; 0 2])
%       vdisp([1 2; 1 0]\{[1 0; 0 1],[1 0; 0 2]})
% 
% See also: cell/ldivide
%
% Written by: tommsch, 2018

error("mrdivide for cells not yet implemented.");

if(iscell(A) && ~iscell(B))
    C=cellfun(@(x) x/B,A,'UniformOutput',false);  %# Apply / cell-wise
elseif(~iscell(A) && iscell(B))
    C=cellfun(@(x) A/x,B,'UniformOutput',false);  %# Apply / cell-wise
elseif(iscell(A) && iscell(B))
  C = cellfun(@rdivide,A,B,'UniformOutput',false);
end
end