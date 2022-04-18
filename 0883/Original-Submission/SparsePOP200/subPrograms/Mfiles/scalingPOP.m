function [nObjPoly,nInEqPolySys,nLbd,nUbd,objValScale,inEqValScale,variableScale] = ...
	scalingPOP(objPoly,inEqPolySys,lbd,ubd)
% 
% A given problem with objPoly,inEqPolySys,lbd and ubd 
% is scaled so that 
%	max{ |lbd(j)|, |ubd(j)| } = 1
% and that the maximum absolute value of coefficients of all monomials in each 
% polynomial function (objective or constraint polynomial function) is 1. 
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
lbdForAbsBounds = 1.0; 
ubdForAbsBounds = 1.0;
lbdForMaxCoef = 1.0; 
ubdForMaxCoef = 1.0;

% Determining variableScale and scaleSW
nDim = length(lbd);
variableScale = ones(1,nDim);
nLbd = lbd;
nUbd = ubd;
for j=1:nDim
    if (-1.0e8 < lbd(j)) && (ubd(j) < 1.0e+8)
        absBound = max(abs(lbd(j)),abs(ubd(j)));
        if absBound > 1.0e-8
            if ubdForAbsBounds < absBound
                variableScale(j) = absBound / ubdForAbsBounds;
                nLbd(j) = lbd(j) / variableScale(j);
                nUbd(j) = ubd(j) / variableScale(j);
            elseif absBound < lbdForAbsBounds
                variableScale(j) = absBound / lbdForAbsBounds;
                nLbd(j) = lbd(j) / variableScale(j);
                nUbd(j) = ubd(j) / variableScale(j);
            end
        end
    end
end
% "Old x(j)" = variableScale * "New x(j)" if scaleSW(j) == 1

funcValScaleSW = 1; 
[nObjPoly,objValScale] = scaleSinglePolynomial(funcValScaleSW,...
    objPoly,variableScale,lbdForMaxCoef,ubdForMaxCoef);

funcValScaleSW = 1; 
mDim = size(inEqPolySys,2); 
inEqValScale = ones(1,mDim);
if mDim == 0
  nInEqPolySys = [];
  inEqValScale = [];
else  
  for i=1:mDim
	[nInEqPolySys{i},inEqValScale(i)] = scaleSinglePolynomial(funcValScaleSW,inEqPolySys{i},...
		variableScale,lbdForMaxCoef,ubdForMaxCoef);
  end
end

debug = 0;
if debug == 1 
	checkScalingPOP(objPoly,inEqPolySys,lbd,ubd, ...  
		nObjPoly,nInEqPolySys,nLbd,nUbd,objValScale,inEqValScale,variableScale);
end
return 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [polyNew,functValScale] = scaleSinglePolynomial(funcValScaleSW,...
    polyOld,variableScale,lbdForMaxCoef,ubdForMaxCoef)

if isempty(polyOld)
  polyNew = [];
  functValScale = [];
  return;
end

noTerms = polyOld.noTerms;

polyNew.typeCone = polyOld.typeCone;
polyNew.sizeCone = polyOld.sizeCone;
polyNew.dimVar = polyOld.dimVar;
polyNew.degree = polyOld.degree;
polyNew.noTerms = polyOld.noTerms;
polyNew.supports = polyOld.supports;
polyNew.coef = polyOld.coef;
maxCoef = 0; 
for i=1:noTerms
	I = find(polyNew.supports(i,:));
	if ~isempty(I)
        ratio = prod(power(variableScale(I),polyNew.supports(i,I)));
	polyNew.coef(i,:) = ratio * polyOld.coef(i,:); 
	end
	maxCoef = max(maxCoef,norm(polyNew.coef(i,:),inf));
end

if funcValScaleSW == 1
	if (maxCoef < lbdForMaxCoef) && abs(maxCoef) > 1.0e-8 
		polyNew.coef = (lbdForMaxCoef / maxCoef) * polyNew.coef;
		functValScale = maxCoef / lbdForMaxCoef; 
	elseif (ubdForMaxCoef < maxCoef) && abs(maxCoef) > 1.0e-8 
		polyNew.coef = (ubdForMaxCoef / maxCoef) * polyNew.coef;
		functValScale = maxCoef / ubdForMaxCoef;  
	else
		functValScale = 1.0;
	end
else
	functValScale = 1.0;
end
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function checkScalingPOP(objPoly,inEqPolySys,lbd,ubd, ... 
	nObjPoly,nInEqPolySys,nLbd,nUbd,objValScale,inEqValScale,variableScale)

% Check lbd, ubd, nLbd and nUbd
[I] = find(lbd > -1.0e9 & lbd < 1.0e9); 
difference = norm(lbd(I) - variableScale(I) .* nLbd(I),inf); 
fprintf('\nCheck lbd and nLbd : %+6.2e\n',difference); 
difference = norm(ubd(I) - variableScale(I) .* nUbd(I),inf); 
fprintf('Check ubd and nUbd : %+6.2e\n',difference); 

% Check objPoly and nObjPoly
nDim = nObjPoly.dimVar;
nxVect = rand(nDim,1);
xVect = variableScale' .* nxVect;
nObjValue = evalPolynomials(nObjPoly,nxVect);
objValue = evalPolynomials(objPoly,xVect);
difference = objValScale * nObjValue -objValue; 
fprintf('Check objPoly and nObjPoy : %+6.2e\n',difference); 

% Check inEqPolySys and nInEqPolySys
mDim = size(inEqPolySys,2);
difference = 0.0; 
for i=1:mDim
	nInEqVal = evalPolynomials(nInEqPolySys{i},nxVect);
	inEqVal = evalPolynomials(inEqPolySys{i},xVect); 
	difference = difference + abs(inEqValScale(i) * nInEqVal -inEqVal);
end
fprintf('Check inEqPolySys and nInEqPolysys : %+6.2e\n', difference); 
return

% $Header: /home/waki9/CVS_DB/SparsePOPdev/subPrograms/scalingPOP.m,v 1.1.1.1 2007/01/11 11:31:50 waki9 Exp $
