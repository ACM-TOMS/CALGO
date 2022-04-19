function basisIndices = genBasisIndices(objPoly,...
    ineqPolySys,cliqueSet,param)

%
% function genBasisIndices
% This function outputs index sets of variables that consists of each
% localizing and moment matrices induced from a maximal clique sets.
%
% Input:
% objPoly --- objective function. Here, we use only this in order to
% dimension of POP.
% ineqPolySys --- constraints of POP.
% cliqueSet --- a maximal clique set found by csp matrix.
% param --- Here, use only param.multiCliquesFactor. This decides the rate
% of mixing some cliques.
%
% Output:
% basisIndices --- Index sets of variables of consisisting of localizing and
% moment matrices.
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

% the number of variables in POP
nDim = objPoly.dimVar;
rowSize = size(ineqPolySys,2);
%the number of maximal cliques of the csp graph induced from POP.
nClique = size(cliqueSet,1);

% the dense case
if param.sparseSW == 0
    idx = 1:nDim;
    for i=1:rowSize
        basisIndices{i} = idx;
    end
    basisIndices{rowSize+1} = idx;
else% the sparse case
    % the maximum size of the maximaz cliques.
    basisIndices = cell(1,rowSize+nClique);
    maxSizeClique = max(sum(cliqueSet,2));
    for i=1:rowSize
        % the n-dimensional row vector whose positive elements indicate
        % nonzeros.
        nzIndicator = any(ineqPolySys{i}.supports,1);
        IntersectClique = cliqueSet*nzIndicator';
        candidates = find(IntersectClique == nnz(nzIndicator));
        % the cliques each of which covers nzIndicator =
        % the candidate of cliques whose union replaces nzIndicator.
        noCandidates = length(candidates);

        nzIndicator = cliqueSet(candidates(1),:);
        noOfNz = nnz(nzIndicator);
        maxSize = maxSizeClique * param.multiCliquesFactor;
        j = 2;
        % expanding nzIndicator until its size does not exceed maxSize.
        while (j <= noCandidates) && (noOfNz < maxSize)
            newNzIndicator = nzIndicator + cliqueSet(candidates(j),:);
            newNoOfNz = nnz(newNzIndicator);
            if newNoOfNz <= maxSize
                nzIndicator = newNzIndicator;
                noOfNz = newNoOfNz;
            end
            j = j+1;
        end
%         for k=1:nDim
%             if(nzIndicator(k) ~= 0)
%                 fprintf('%d ',k-1);
%             end
%         end
%         fprintf('\n');
%         if i == rowSize -2
%             fprintf('\n noCandidates = %d\n',noCandidates);
%         end
        basisIndices{i} = find(nzIndicator);
    end
    for i=(1:nClique)+rowSize
        basisIndices{i} = find(cliqueSet(i-rowSize,:));
    end
end


return

% $Header: /home/waki9/CVS_DB/SparsePOPdev/subPrograms/genBasisIndices.m,v 1.2 2007/01/31 00:58:58 waki9 Exp $
