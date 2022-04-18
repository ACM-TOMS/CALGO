function C = ldivide(A,B)
% times(A,B)
% Scalar * Cell : mldivide elementwise
% Cell * Scalar : mldivide elementwise
% Cell * Cell   : ldivide elementwise
%
% See also: cell/mtimes
%
% Written by: tommsch, 2018

if(iscell(A) && ~iscell(B))
    C=cellfun(@(x) x\B,A,'UniformOutput',false);  %# Apply mtimes cell-wise
elseif(~iscell(A) && iscell(B))
    C=cellfun(@(x) A\x,B,'UniformOutput',false);  %# Apply mtimes cell-wise
elseif(iscell(A) && iscell(B))
  C = cellfun(@(x,y) x.\y,A,B,'UniformOutput',false);
end

end