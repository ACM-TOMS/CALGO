function [objPoly,ineqPolySys,lbd,ubd,trans] ... 
	= convert2(objPoly0,ineqPolySys0,origLbd,origUbd)
% 
% Kojima September 11, 2005
%
% Symbolic Math Tool is necessary to execute this function. 
%
% Convert a given polynomial to the one with positive variables of lower 
% bound 0 and free variables. 
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
%
% Variable type for each variable is usesd for identifying the conversion needed.
% Variable type 
%        0: free varable with no lower and upper bounds
%	 	-1: variable with a lower bound but no upper bound
%		+1: variable with no lower bound but upper bound
%		+2: variable with lower and upper bounds

nDim = objPoly0.dimVar;

%===========================================================================
% Kojima 09/11/2005 ----> 
% Linear transformation from the new variable xVecNew to the old variable xVecOld
% xVecOld = trans.Amat * xVecNew + trans.bVect;

% The variable indices are classified into 
%   freeVar \cup zeroLowerBdd \cup varToConvert
% varToConvert is the union of lowBddVar and onlyUpBddVar. 
% Bur 
%   freeVar = find(origLbd < -0.999e10 & +0.999e+10 <= origUbd); 
% is not used. 
zeroLowBddVar = find(-1.0e-10  < origLbd &  origLbd < +1.0e-10);
varToConvert = find((-1.0e+8 < origLbd | origUbd < 1.0e+8)... 
    & (origLbd <= -1.0e-10 | +1.0e-10 <= origLbd));
lowBddVar   = find((-1.0e+8 < origLbd) & (origLbd <= -1.0e-10 | +1.0e-10 <= origLbd));
onlyUpBddVar = find(origLbd <=  -1.0e+8 & origUbd < 1.0e+8);
% trans.Amat and trans.bVect
trans.Amat = speye(nDim); 
trans.bVect = zeros(nDim,1);

for i=onlyUpBddVar
	trans.Amat(i,i) = -1.0;
end
trans.bVect(lowBddVar,1) = origLbd(lowBddVar)';
trans.bVect(onlyUpBddVar,1) = origUbd(onlyUpBddVar)';
% lower and upper bounds 
lbd = -1.0e10*ones(1,nDim);
ubd = 1.0e10*ones(1,nDim);
if length(zeroLowBddVar) > 0 
    lbd(zeroLowBddVar) = zeros(1,length(zeroLowBddVar));
    ubd(zeroLowBddVar) = origUbd(zeroLowBddVar);
end
for i=lowBddVar
	if origUbd(i) >= 1.0e+8
		lbd(i) = 0;
	end
end
if length(varToConvert) > 0 
    lbd(varToConvert) = zeros(1,length(varToConvert));
    if length(lowBddVar) > 0
		for j=lowBddVar
			if origUbd(j) < 1.0e+8 
        		ubd(j) = origUbd(j) - origLbd(j);
			else
				ubd(j) = 1.0e+10;
			end
    	end
	end
    if length(onlyUpBddVar) > 0 
        %ubd(onlyUpBddVar) = origUbd(onlyUpBddVar);
        ubd(onlyUpBddVar) = 1.0e+10;
    end
end
%{
[origLbd;origUbd]
zeroLowBddVar
varToConvert
lowBddVar
onlyUpBddVar
[lbd;ubd]
%}

tempSW = 1;
syms y real

% objPoly
objPoly0 = simplifyPolynomial(objPoly0);
objPoly.typeCone = objPoly0.typeCone; 
objPoly.sizeCone = objPoly0.sizeCone; 
objPoly.dimVar = objPoly0.dimVar; 
objPoly.degree = objPoly0.degree;
if isempty(varToConvert)
    objPoly.noTerms = objPoly0.noTerms; 
    objPoly.supports = objPoly0.supports; 
    objPoly.coef = objPoly0.coef;
else 
    for j=varToConvert 
        objPoly.noTerms = 0; 
        objPoly.supports = []; 
        objPoly.coef = [];
        for i=1:objPoly0.noTerms          
            pp = objPoly0.supports(i,j); 
            if pp >= 1 
                cc = sym2poly((trans.Amat(j,j)*y+trans.bVect(j,1))^pp)';
                supSet = objPoly0.supports(i*ones(pp+1,1),:);
                supSet(:,j) = (pp:-1:0)';
%                
                objPoly.noTerms = objPoly.noTerms + size(supSet,1); 
                objPoly.supports = [objPoly.supports; supSet];
                objPoly.coef = [objPoly.coef; cc*objPoly0.coef(i,:)];
            else
                objPoly.noTerms = objPoly.noTerms + 1; 
                objPoly.supports = [objPoly.supports; objPoly0.supports(i,:)];
                objPoly.coef = [objPoly.coef; objPoly0.coef(i,:)];
            end
        end
        if tempSW == 1 
            objPoly0 = simplifyPolynomial(objPoly);
        else 
            objPoly0.noTerms = objPoly.noTerms; 
            objPoly0.supports = objPoly.supports; 
            objPoly0.coef = objPoly.coef;
        end
    end
    if tempSW == 1 
        objPoly.noTerms = objPoly0.noTerms; 
        objPoly.supports = objPoly0.supports; 
        objPoly.coef = objPoly0.coef;
    else 
        objPoly =simplifyPolynomial(objPoly);
    end
end

%full([objPoly.supports, objPoly.coef])

% ineqPolySys; 
if ~isempty(ineqPolySys0)
    m = size(ineqPolySys0,2);
    ineqPolySys = cell(1,m);
    for k=1:m 
        ineqPolySys0{k} = simplifyPolynomial(ineqPolySys0{k});
        ineqPolySys{k}.typeCone = ineqPolySys0{k}.typeCone; 
        ineqPolySys{k}.sizeCone = ineqPolySys0{k}.sizeCone; 
        ineqPolySys{k}.dimVar = ineqPolySys0{k}.dimVar; 
        ineqPolySys{k}.degree = ineqPolySys0{k}.degree;
        if isempty(varToConvert)
            ineqPolySys{k}.noTerms = ineqPolySys0{k}.noTerms; 
            ineqPolySys{k}.supports = ineqPolySys0{k}.supports; 
            ineqPolySys{k}.coef = ineqPolySys0{k}.coef;
        else 
            for j=varToConvert 
                ineqPolySys{k}.noTerms = 0; 
                ineqPolySys{k}.supports = []; 
                ineqPolySys{k}.coef = [];
                for i=1:ineqPolySys0{k}.noTerms
                    pp = ineqPolySys0{k}.supports(i,j); 
                    if pp >= 1 
                        cc = sym2poly((trans.Amat(j,j)*y+trans.bVect(j,1))^pp)';
                        supSet = ineqPolySys0{k}.supports(i*ones(pp+1,1),:);
                        supSet(:,j) = (pp:-1:0)';
%                
                        ineqPolySys{k}.noTerms = ineqPolySys{k}.noTerms + size(supSet,1); 
                        ineqPolySys{k}.supports = [ineqPolySys{k}.supports; supSet];
                        ineqPolySys{k}.coef = [ineqPolySys{k}.coef; cc*ineqPolySys0{k}.coef(i,:)];
                    else
                        ineqPolySys{k}.noTerms = ineqPolySys{k}.noTerms + 1; 
                        ineqPolySys{k}.supports = [ineqPolySys{k}.supports; ineqPolySys0{k}.supports(i,:)];
                        ineqPolySys{k}.coef = [ineqPolySys{k}.coef; ineqPolySys0{k}.coef(i,:)];
                    end
                end
                if tempSW == 1
                    ineqPolySys0{k} = simplifyPolynomial(ineqPolySys{k});
                else 
                    ineqPolySys0{k}.noTerms = ineqPolySys{k}.noTerms; 
                    ineqPolySys0{k}.supports = ineqPolySys{k}.supports; 
                    ineqPolySys0{k}.coef = ineqPolySys{k}.coef;
                end
            end
            if tempSW == 1
                ineqPolySys{k}.noTerms = ineqPolySys0{k}.noTerms; 
                ineqPolySys{k}.supports = ineqPolySys0{k}.supports; 
                ineqPolySys{k}.coef = ineqPolySys0{k}.coef;
            else
                ineqPolySys{k} = simplifyPolynomial(ineqPolySys{k});
            end
        end
    end
else
    ineqPolySys = [];
end

clear y; 
  
return;



% $Header: /home/waki9/CVS_DB/SparsePOPdev/subPrograms/convert2.m,v 1.2 2007/01/12 09:34:38 waki9 Exp $
