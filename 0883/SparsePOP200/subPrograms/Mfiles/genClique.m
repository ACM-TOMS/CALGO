function clique = genClique(objPoly,inEqPolySys,sparseSW)

% FUNCTION genClique
% This function outputs the information about maximal cliques of
% the correlative sparrsity pattern graph for given POP.
%
% Procedure to find the maximal cliques:
% Step 1: generate correlative sparsity pattern matrix
% Step 2: apply chordal extension to the csp matrix
% Step 3: find the maximal cliques of the extension
%
% <Input argument>
% objPoly    : objective function
% inEqPOlySys: constraints
% sparseSW   : a kind of SDP relaxation
%              if sparseSW = 1, this finds cliques.
%
% <Output argument>
% clique: information on the maximal cliques, this variable is the
%         structure with the following fields.
%
% clique.Set:  The maximal clieques, which are represented by 0-1.
%              e.g. The number of variable is 3 and the maximal 
%                   cliques are C1={1,2} and C2={2,3}.
%                   
%                   Then, clique.Set(1,:) = [1,1,0];(= C1)
%                         clique.Set(2,:) = [0,1,1];(= C2)
% clique.NoC:  The number of cliques
% clique.maxC: The maximum size in the maximum cliques
% clique.minC: The minimum size in the maximum cliques
%

%
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
if sparseSW == 1
  %%
  %% Step 1:
  %% generate csp matrix by gathering all the supports.
  %%
  tSup = objPoly.supports;
  tSup = spones(tSup);
  tempSup = unique(tSup,'rows');% to simplify support
  m = size(tempSup,1);
  
  %% construct csp matrix
  rmat = speye(nDim);
  for i=1:m
    %%
    %% objPoly: if x_p & x_q in monomial, rmat_pq = rand
    %%
    nonIdx = find(tempSup(i,:));
    rmat(nonIdx,nonIdx) = rand(1,1);
  end
  clear tempSup;
 
  noOfInequalities = size(inEqPolySys,2);
  for i=1:noOfInequalities
    if inEqPolySys{i}.dimVar == nDim
      %%
      %%inEqPolySys{i}: if x_p & x_q in support, rmat_pq = rand
      %%
      tSup = any(inEqPolySys{i}.supports,1);
      nonIdx = find(tSup);
      rmat(nonIdx,nonIdx) = rand(1,1);
    else
      %% check the dimension
      error('Warning: Please set the same as dimension of objective func. to inEqPolySys{%d}.',i);
    end
  end  
  %% make rmat diagonally dominant
  rmat = rmat + 5*nDim*speye(nDim);

  clear I J;
  %%
  %% Step 2
  %% Chordal Extension by Cholesky decomposition
  %%
  %% minimum degree ordering     
  I = symamd(rmat);
  %% cholesky decomposition
  [R,p] = chol(rmat(I,I));	
  if (p > 0) 
    error('Correlative sparsity matrix is not positive definite.');
  end
  
  %%
  %% Step3
  %% Finding the maxmal clieques
  %%
  
  %% put 1 for nonzero element ofR
  Cliques = spones(R);
  [value,orig_idx] = sort(I);
  remainIdx = 1;
  for i=2:nDim
    checkSet = Cliques(i,i:nDim);
    one = find(checkSet);
    noOfone = length(one);
    cliqueResult = Cliques(1:i-1,i:nDim)*checkSet';
    yesno = find(cliqueResult == noOfone);
    %%
    %% Remove the set included into other set.
    %%
    if ~any(yesno)
      remainIdx = [remainIdx;i];
    end
  end
  clique.Set = Cliques(remainIdx,orig_idx);
else
  clique.Set = ones(1,nDim);
end
%%
%%Clique Information
%%
clique.NoC  = size(clique.Set,1);
sumClique = sum(clique.Set,2);
clique.maxC = max(sumClique);
clique.minC = min(sumClique);

return;

% $Header: /home/waki9/CVS_DB/SparsePOPdev/subPrograms/genClique.m,v 1.1.1.1 2007/01/11 11:31:50 waki9 Exp $
