function [CompSup,ConstraintInfo] = separateSpecialMonomial(inEqPolySys,param)
% Kojima, 02/15/2005
% Separate_Special_Monomial ---> separateSpecialMonomial
% function [ConstraintInfo] = ...
%    Separate_Special_Monomial(inEqPolySys);
%
%substituteEqList --- list of equality that has only 2 or 1 monomial

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
ConstraintInfo.Term1.comp = [];
ConstraintInfo.Term2.combzero = [];
ConstraintInfo.Term2.combsquare = [];
ConstraintInfo.Term2.equiv = [];
ConstraintInfo.Term2.constant = [];
ConstraintInfo.Term2.others = [];
ConstraintInfo.Others = [];

EPS = 1.0e-8;
noOfInequality = size(inEqPolySys,2);
if param.complementaritySW == 1
    for i=1:noOfInequality
        typeCone = inEqPolySys{i}.typeCone;
        if typeCone == -1
            for j=1:inEqPolySys{i}.sizeCone
                %I = find(abs(inEqPolySys{i}.coef(:,j)) > EPS);
                %noTerms = length(I);
                noTerms = inEqPolySys{i}.noTerms;
				if noTerms == 1%% Complementality Constraint
                    ConstraintInfo.Term1.comp = [ConstraintInfo.Term1.comp,i];
	                %[i,noTerms]
					%full([inEqPolySys{i}.supports, inEqPolySys{i}.coef])
                elseif noTerms == 2
                    a1 = inEqPolySys{i}.coef(1,j);
                    a2 = inEqPolySys{i}.coef(2,j);
                    if a1*a2 == 0%% Complementality Constraint
                        ConstraintInfo.Term1.comp = [ConstraintInfo.Term1.comp,i];
                		%[i,noTerms]
						%full([inEqPolySys{i}.supports, inEqPolySys{i}.coef])
                        %elseif a1*a2 == -1%% Combinatorial Constraint
                        %substituteEqList.Terms2  = [substituteEqList.Terms2;i];
                    else
                        ConstraintInfo.Others = [ConstraintInfo.Others, i];
                    end
                else
                    ConstraintInfo.Others = [ConstraintInfo.Others, i];
                end
            end
        end
    end
end
if (~isempty(inEqPolySys)) && (param.complementaritySW == 1)
    CompSup = SupTerm1(inEqPolySys,ConstraintInfo);
else
    CompSup = [];
end
return

function CompSup = SupTerm1(inEqPolySys,ConstraintInfo)

EPS = 1.0e-10;
LenTerm1 = length(ConstraintInfo.Term1.comp);
nDim = inEqPolySys{1}.dimVar;
CompSup = sparse(LenTerm1,nDim);
t = 0;
for j=1:LenTerm1
    i= ConstraintInfo.Term1.comp(j);
    for k=1:inEqPolySys{i}.sizeCone
        I = find(abs(inEqPolySys{i}.coef(:,k)) > EPS);
        if length(I) == 1
            t = t + 1;
            CompSup(t,:) = inEqPolySys{i}.supports;
        end
    end
end
return


% $Header: /home/waki9/CVS_DB/SparsePOPdev/subPrograms/separateSpecialMonomial.m,v 1.1.1.1 2007/01/11 11:31:50 waki9 Exp $
