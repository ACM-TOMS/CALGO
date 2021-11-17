function [ineqPolySys,lbdIdx,ubdIdx] = boundToIneqPolySys(ineqPolySys,lbd,ubd) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Incorprate lower and upper bounds into ineqPolySys to strengthen 
% the relaxation
%
% lbdIdx is cell. lbdIdx{i} indicates the number of ineqPolysys that
% represents constraint "x_i >= lbd(i)".
%
% ubdIdx is cell. ubdIdx{i} indicates the number of ineqPolysys that
% represents constraint "x_i <= ubd(i)".
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
pointer = size(ineqPolySys,2); 
noOfVariables = length(lbd); 
[ineqPolySys,lbdIdx,pointer] = lbdToInEqPolySys(ineqPolySys,pointer,...
    noOfVariables,lbd);
[ineqPolySys,ubdIdx,pointer] = ubdToInEqPolySys(ineqPolySys,pointer,...
    noOfVariables,ubd);
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ineqPolySys,lbdIdx,pointer] = lbdToInEqPolySys(ineqPolySys,pointer,...
    noOfVariables,lbd)
epsilon = 1.0e-8;
lbdIdx = cell(1,noOfVariables);
for i=1:noOfVariables
    if lbd(i) > -1.0e8
        pointer = pointer + 1;
        ineqPolySys{pointer}.typeCone = 1;
        ineqPolySys{pointer}.sizeCone = 1;
        ineqPolySys{pointer}.degree = 1;
        ineqPolySys{pointer}.dimVar = noOfVariables;
        lbdIdx{i} = pointer;
        if abs(lbd(i)) >= epsilon
            ineqPolySys{pointer}.noTerms = 2;
            ineqPolySys{pointer}.supports = sparse(2,...
                ineqPolySys{pointer}.dimVar);
            ineqPolySys{pointer}.supports(2,i) = 1;
            ineqPolySys{pointer}.coef = [-lbd(i);1];
        else
            ineqPolySys{pointer}.noTerms = 1;
            ineqPolySys{pointer}.supports = sparse(1,...
                ineqPolySys{pointer}.dimVar);
            ineqPolySys{pointer}.supports(1,i) = 1;
            ineqPolySys{pointer}.coef = 1;
        end
    end
end

return

function [ineqPolySys,ubdIdx,pointer] = ubdToInEqPolySys(ineqPolySys,...
    pointer,noOfVariables,ubd)
ubdIdx = cell(1,noOfVariables);
for i=1:noOfVariables
    if ubd(i) < 1.0e8
        pointer = pointer + 1;
        ineqPolySys{pointer}.typeCone = 1;
        ineqPolySys{pointer}.sizeCone = 1;
        ineqPolySys{pointer}.degree = 1;
        ineqPolySys{pointer}.dimVar = noOfVariables;
        ineqPolySys{pointer}.noTerms = 2;
        ineqPolySys{pointer}.supports = sparse(2,...
            ineqPolySys{pointer}.dimVar);
        ineqPolySys{pointer}.supports(2,i) = 1;
        ineqPolySys{pointer}.coef = [ubd(i);-1];
        ubdIdx{i} = pointer;
    end
end
return


% $Header: /home/waki9/CVS_DB/SparsePOPdev/subPrograms/boundToIneqPolySys.m,v 1.2 2007/01/12 09:45:15 waki9 Exp $
