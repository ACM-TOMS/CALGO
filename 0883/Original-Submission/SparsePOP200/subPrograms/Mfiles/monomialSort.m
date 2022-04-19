% Kojima, 02/15/005
% Monomial_Sort.m ---> monomialSort.m
function [BMat, I] = monomialSort(AMat, subspaceIdxSet, order)

%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This file is a component of SparsePOP 
% Copyright (C) 2007 SparsePOP Project
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307 USA
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
dimSubspace = length(subspaceIdxSet); 
rowSize = size(AMat,1);
if rowSize <= 1
  BMat = AMat;
  I = 1; 
  return;
elseif dimSubspace == 1 
  [tempCol,I] = sort(AMat(:,subspaceIdxSet(1))); 
  BMat = AMat(I,:);
  return;
elseif   nargin < 3 || strcmp(order, 'lex') 
  [BMat,I] = lexicoSort(AMat,subspaceIdxSet);
elseif strcmp(order, 'grevlex')
  [BMat,I] = GraededReverseLexSort(AMat,subspaceIdxSet);
elseif strcmp(order, 'grlex')
  [BMat,I] = GraededLexSort(AMat,subspaceIdxSet);
end
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [BMat,I] = lexicoSort(AMat,subspaceIdxSet)
% Sort the m row vectors of the matrix AMat using row indices startIndex,\ldots,n, 
% where [m,n] = size(AMat). 
% For example, if AMat = [1,4,3; 1,2,5; 1,5,1] and startIndex = 1 or 2 then 
% BMat = [1,2,5;1,4,3;1,5,1]. If startIndex = 3 then BMat  = [1,5,1;1,4,3;1,2,5]

[BMat,I] = sortrows(AMat,subspaceIdxSet);
return;

function [BMat,I] = GraededReverseLexSort(AMat,subspaceIdxSet)
%%Graded Reverse Lex Order (increase case)

degree = sum(AMat(:,subspaceIdxSet),2);
CMat = [degree,-AMat];
[DMat,I] = lexicoSort(CMat,[1,subspaceIdxSet+1]);
BMat = -DMat(:,2:end);
return ;

function [BMat,I] = GraededLexSort(AMat,subspaceIdxSet)
%%Graded Lex Order (increase case)

degree = sum(AMat(:,subspaceIdxSet),2);
CMat = [degree,AMat];
[DMat,I] = lexicoSort(CMat,[1,subspaceIdxSet+1]);
BMat = DMat(:,2:end);
return ;


% $Header: /home/waki9/CVS_DB/SparsePOPdev/subPrograms/monomialSort.m,v 1.1.1.1 2007/01/11 11:31:50 waki9 Exp $
