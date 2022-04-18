function C = minus(A,B)
% C = minus(A,B)
% Computes A-B for cells
% Input: 
%   A and B are cell arrays                     Computes A minus B, Same Behaviour as minus for arrays
%   A or B is a cell array, the other is not    Compute Cell minus Element, for all elements in the cell
%

% See also: cell/plus, minus
%
% Written by: tommsch, 2018

C=A+(-B);

end