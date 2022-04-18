%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This module contains 
%	genBasisSupports;
% 	genSimplexSupport;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [basisSupports] = genBasisSupports(objPoly,ineqPolySys,param,...
    basisIndices)
%========================================================================
% If ineqPolySys{i}.typeCone = 1 then u_{g_i}u_{g_i}^T will be multiplied
% to the constraint polynomial ineqPolySys{i},
% where g_i = basisSupports{i}. 
%
% If ineqPolySys{i}.typeCone = -1 then u_{g_i} will be multiplied to the 
% constraint polynomial ineqPolySys{i}, where g_i =  basisSupports{i}. 
%
% Each  basisSupports{i} is determined by sosDim and basisIndices{i},
% where basisIndices{i} denotes the index set of coordinates along 
% which basisSupports{i} can be positive and sosDim denotes the degree 
% of the support for u_{g_i}
%
% relaxOrder means the degree of the support u_{g_j} u_{g_j}^T 
% j=m+1,...m', where m' = the number of cliques + m. 
% Relations: 
%		param.relaxOrder >= ceil((the max degree of f_j)/2) 
%       (j=0,1,2,\ldots,m). 
%		param.relaxOrder-1 <= 2*sosDim + deg f_j <= relaxOrder
%       (j=1,2,\ldots,m).
%========================================================================

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
nDim = objPoly.dimVar;
if isfield(param,'relaxOrder')
    rO = param.relaxOrder;
else
    error('Please define param.relaxOrder.');
end
m = size(ineqPolySys,2);
ell = size(basisIndices,2);
basisSupports = cell(1,ell);
for i=1:m
    sosDim = rO - ceil(ineqPolySys{i}.degree/2);
    if sosDim < 0
        error('## Increase relax Order##\n');
    end
    if (ineqPolySys{i}.typeCone == -1)
        %if (mod(ineqPolySys{i}.degree,2) == 0)
        basisSupports{i} = genSimplexSupport(1,nDim,2*sosDim,...
                basisIndices{i},'grevlex');
        %else
        %basisSupports{i} = genSimplexSupport(1,nDim,2*sosDim+1,...
        %        basisIndices{i},'grevlex');
        %end
    else
        basisSupports{i} = genSimplexSupport(1,nDim,sosDim,...
            basisIndices{i},'grevlex');
    end
end

for i=(m+1):ell
	basisSupports{i} = genSimplexSupport(1,nDim,rO,...
        basisIndices{i},'grevlex');
end

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function supSet = genSimplexSupport(supType,nDim,r,subspaceIdxSet,order)
% supType = 0 or 1
% supspaceIdxSet \subset 1:nDim
%
% If supType == 0, then supSet consists of the elements 
% \{ v \in Z^n_+ : \sum_{i=1}^n v_i = r, v_j = 0 
% (j \not\in supSpaceIdxSet \}
% in the lexico graphical order. 
%
% If supType == 1 then supSet consists of the elements 
% \{ v \in Z^n_+ : \sum_{i=1}^n v_i \le r, v_j = 0 
% (j \not\in supSpaceIdxSet \}
% in the lexico graphical order. 
% 
if nargin < 5
  order = 'grevlex';
end
dimSubspace = length(subspaceIdxSet); 
if nDim <= 0
  error('!!! nDim <= 0 !!!');
elseif r < 0
  error('!!! r < 0 !!!');	
elseif dimSubspace == nDim
  if subspaceIdxSet(nDim) == nDim
    if supType == 0
      supSet = flatSimpSup(nDim,r);
    elseif supType == 1
      supSet = fullSimpSup(nDim,r);
    else
      error('!!! supType is neither 0 nor 1 !!!');
    end
  else
    error('!!! dimSubspace = nDim but subspaceIdxSet(nDim) not= nDim !!!');
  end
else
  supSet = restSimpSup(supType,nDim,r,subspaceIdxSet);
end

% Kojima, 02/15/2005
% [supSet,I] = Monomial_Sort(supSet, subspaceIdxSet, order);
[supSet,I] = monomialSort(supSet, subspaceIdxSet, order);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function supSet = flatSimpSup(nDim,r)
% \{ v \in Z^n_+ : \sum_{i=1}^n v_i = r \}

if nDim == 0
  supSet = 0; 
elseif r == 0
  supSet = sparse(1,nDim);
elseif nDim == 1 
  supSet = r; 
else
  NumElem = nchoosek(nDim+r-1,r);
  supSet = sparse(NumElem,nDim);
  index = 0;
  for i=0:r
    aSupSet = flatSimpSup(nDim-1,i);
    [m,n] = size(aSupSet);
    Idx = index + (1:m);
    supSet(Idx,1) = repmat(r-i,m,1);
    supSet(Idx,2:n+1) = aSupSet;
    index = index + m;
  end
end
return; 		

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function supSet = fullSimpSup(nDim,r)
% \{ v \in Z^n_+ : \sum_{i=1}^nDim v_i \leq r \}

if nDim == 1;
  supSet = (0:r)';
  supSet = sparse(supSet);
elseif r == 0
  supSet = sparse(1,nDim);
else
  NumElem = nchoosek(nDim+r,r);
  supSet = sparse(NumElem,nDim);
  index = 0;
  for i=0:r
    aSupSet = flatSimpSup(nDim,i);
    m = size(aSupSet,1);
    Idx = index + (1:m);
    supSet(Idx,:) = aSupSet;
    index = index + m;
  end
end
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function supSet = restSimpSup(supType,nDim,r,subspaceIdxSet)

dimSubspace = length(subspaceIdxSet); 
%
% error handlng
%
if (nDim < dimSubspace) || (nDim < subspaceIdxSet(dimSubspace)) 
  error('!!! nDim < dimSubspace !!!');
end
if nDim == 0;
  error('!!! nDim = 0 !!!'); 
end
if dimSubspace == 0;
  error('!!! dimSupspace = 0 !!!');
end
if nDim < subspaceIdxSet(dimSubspace)
  error('!!! nDim < subspaceIdxSet(dimSubspace) !!!'); 
end
%
% end of error handling
%

if supType == 1
  %\{ v \in Z^n_+ : \sum_{i=1}^n v_i \le r, v_j = 0
  % (j \not\in supSpaceIdxSet) \}
  aSupSet = fullSimpSup(dimSubspace,r);
elseif supType ==0
  %\{ v \in Z^n_+ : \sum_{i=1}^n v_i = r, v_j = 0 
  % (j \not\in supSpaceIdxSet) \}
  aSupSet = flatSimpSup(dimSubspace,r); 
else
  error('!!! You should choose 1 or 0 as supType !!!'); 
end    
m = size(aSupSet,1); 
supSet = sparse(m,nDim);
supSet(:,subspaceIdxSet) = aSupSet;

return;

% $Header: /home/waki9/CVS_DB/SparsePOPdev/subPrograms/genBasisSupports.m,v 1.2 2007/01/31 00:57:20 waki9 Exp $
