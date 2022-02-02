%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Modified by Kojima, 02/06/2005
% This module includes
% 	subStituteEq(basisSupports,ineqPolySys,ConstraintInfo);
% M. Kojima, 02/15/2005
% SubStituteEq ---> substituteEq
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
function [basisSupports,momentSup,ineqBasis] = substituteEq(basisSupports,ineqBasis,ineqPolySys,CompSup,param)

% This function inserts simple equations into basisSupports,
% not Moment matrix.
% 
% In this version, this does only for complementarities


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

noOfBasisSupports = size(basisSupports,2);
[m,n] = size(CompSup);
if m*n == 0
   param.complementaritySW = 0; 
end

%%
%% Substitute complementarities into all basisSupports.
%%
if param.complementaritySW == 1
    for i=1:noOfBasisSupports
        SupSet = basisSupports{i};
        SupSet = spones(SupSet);
        P = SupSet*CompSup';
        P = max(P,[],2);
        I = find(P == 2);
        if ~isempty(I)
            basisSupports{i}(I,:) = [];
        end
    end
end
%%
%% Substitute into monomials appeared in localizing mat.
%%
if param.boundSW > 0 && param.complementaritySW == 1
    SupSet = ineqBasis;
    SupSet = spones(SupSet);
    P = SupSet*CompSup';
    P = max(P,[],2);
    I = find(P == 2);
    if ~isempty(I)
       ineqBasis(I,:) = []; 
    end
end

%%
%% Gather all monomials appeared in Moment Mat. and 
%% if monomials can be removed by complimentarity, we can do it!
%%
if param.boundSW > 0
    t = size(basisSupports,2);
    s = size(ineqPolySys,2);
    momentSup = [];
    for j = s+1:t
        SupSet = makeMoment(basisSupports{j});
        if param.complementaritySW == 1
            tempSet = spones(SupSet);
            P = tempSet*CompSup';
            P = max(P,[],2);
            I = find(P == 2);
            if ~isempty(I)
                SupSet(I,:) = [];
            end
        end
        momentSup = [momentSup;SupSet];
    end
else
    momentSup = [];
end


return

function Sup = makeMoment(basisSupports)

% This function returns monomials appeared in Moment matrix and
% their row and column indeces.
%
% Sup is monomial set.

[m1,nDim] = size(basisSupports);
[Col,Row] = find(tril(ones(m1)));
Sup = sparse(length(Row),nDim);
UsedVar = find(any(basisSupports,1));
Sup(:,UsedVar) = basisSupports(Row,UsedVar) + basisSupports(Col,UsedVar);
return


% $Header: /home/waki9/CVS_DB/SparsePOPdev/subPrograms/substituteEq.m,v 1.1.1.1 2007/01/11 11:31:50 waki9 Exp $
