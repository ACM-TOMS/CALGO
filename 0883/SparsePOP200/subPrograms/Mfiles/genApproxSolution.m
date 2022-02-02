%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% M. Kojima, 02/06/2005
% Revised M. Kojima, 04/24/2005
% This module includes:
%	genApproxSolution(y,linearterms,SDPobjValue,objPoly,objPoly0,ineqPolySys0,param,trans,fixedVar);
% 	infeasibility(inEqPolySys,xVect);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
function [POP,SDPobjValue] = ...
    genApproxSolution(y,linearterms,SDPobjValue,objPoly,objPoly0,...
    ineqPolySys0,lbd0,ubd0,param,trans,fixedVar)

nDim = objPoly0.dimVar;
removeIdx = [];
remainIdx = (1:nDim);
if ~isempty(fixedVar)
    removeIdx = fixedVar(:,1);
    remainIdx = setdiff((1:nDim),fixedVar(:,1));
end

if length(linearterms) ~= objPoly.dimVar;
    POP.xVect = [];
    SDPobjValue = SDPobjValue * trans.objValScale+ trans.objConstant;
    POP.absError = [];
    POP.scaledError = [];
	return;
else
    POP.xVect = zeros(nDim,1);
    POP.xVect(remainIdx) = y(linearterms);
    if ~isempty(fixedVar)
        POP.xVect(removeIdx) = fixedVar(:,2);
    end
end

if objPoly.noTerms == 0
	POP.objValueScaled = 0;
else
	[POP.objValueScaled,temp] = evalPolynomials(objPoly,POP.xVect(remainIdx));
end

if param.scalingSW == 1 
    SDPobjValue = SDPobjValue * trans.objValScale;
	POP.objValueScaled = POP.objValueScaled*trans.objValScale;
    if ~isempty(POP.xVect)
        POP.xVect(remainIdx) = trans.Amat * POP.xVect(remainIdx);
        POP.xVect(remainIdx) = POP.xVect(remainIdx) + trans.bVect;
    end
end

if isfield(param,'sdpaDataFile') && ischar(param.sdpaDataFile) && ~isempty(param.sdpaDataFile) 
	idx = findstr(param.sdpaDataFile,'.dat-s');
	if ~isempty(idx)
		infoFile = strcat(param.sdpaDataFile(1:idx),'info');
		fid = fopen(infoFile,'w+');
		%%fprintf(fid,'%10.8f\t# scaling coefficient for objective function\n',trans.objValScale);
		%%fprintf(fid,'%10.8f\t# constant coefficient for objective function\n',trans.objConstant);
		DiagMat = zeros(nDim,1);	
		Order   = -1*ones(nDim,1);
		Vect    = zeros(nDim,1);
		if ~isempty(fixedVar)	
			DiagMat(remainIdx) = diag(trans.Amat);
			Order(remainIdx) = linearterms;
			Vect(removeIdx) = fixedVar(:,2);
			Vect(remainIdx) = trans.bVect;
		else
			DiagMat = diag(trans.Amat);
			Order = linearterms;
			Vect = trans.bVect;
		end
		
		for i=1:nDim
			fprintf(fid,'%2d\t%10.8f\t%10.8f\n',Order(i),DiagMat(i),Vect(i));
		end
		fclose(fid);
	end
end

SDPobjValue = SDPobjValue + trans.objConstant;
POP.objValueScaled = POP.objValueScaled + trans.objConstant;
% the objective value of the original POP
%[POP.objValue,maxAbsMonomial] = evalPolynomials(objPoly0,POP.xVect);
POP.objValue = POP.objValueScaled;

[ineqPolySys0,lbdIdx,ubdIdx] = boundToIneqPolySys(ineqPolySys0,lbd0,ubd0);
% the infeasiblity of the original POP

if size(ineqPolySys0,2) == 0
    POP.absError = 0.0;
    POP.scaledError = 0.0;
elseif isempty(POP.xVect)
    POP.absError = [];
    POP.scaledError = [];
else
    [POP.absError,POP.scaledError] = infeasibility(ineqPolySys0,POP.xVect);
end
return;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [infeasError,scaledError] = infeasibility(inEqPolySys,xVect)
%%%%%%
% Check the feasibility of xVect.
% If infeasError >= 0 and scaledError >= 0
% then xVect is feasible in given region.
%%%%%%
% Modified by M. Kojima, 01/06/05.
%
% The modified version outputs a relative infeasibility error,
% scaledError.
%
% Let maxAbsMonomial(f,x) be the maximum over the absolute values of all
% monimials evaluated at x if the maximum is greater than 1 or
%  maxAbsMonomial(f,x) = 1 otherweise;
% 	maxAbsMonomial(f,x)
%	= max { | c_{\alpha} x^{\alpha} | (\alpha \in \FC), 1 },
% where we assume f(x) = \sum_{\alpha \in \FC} c_{\alpha} x^{\alpha}.
% If inequalities f_i(x) \geq 0 are given as input, then
% 	infeasError = min_i { min{f_i(x),0} },
% 	scaledError = min_i { min{f_i(x),0} / max{1,maxAbsMonomial(f,x)} }.
% If inequalities f_i(x) = 0 are given as input, then
%       infeasError = min_i { - |f_i(x)| },
%       scaledError = min_i { - |f_i(x)| / max{1, maxAbsMonomial(f,x)} }.
%
%%%%%%
noOfinEqPolySys = size(inEqPolySys,2);
infeasError = 0;
scaledError = 0;
for i=1:noOfinEqPolySys
    if inEqPolySys{i}.typeCone == 1
        [value,maxAbsMonomial] = evalPolynomials(inEqPolySys{i},xVect);
        maxAbsMonomial = max(1, maxAbsMonomial);
        %fprintf('%d --- %f\n',i,value);
        scaledError = min(scaledError,value/maxAbsMonomial);
        infeasError = min(infeasError,value);
    elseif inEqPolySys{i}.typeCone == 2
        vecSize = inEqPolySys{i}.sizeCone;
        [funcValues,maxAbsMonomial] = evalPolynomials(inEqPolySys{i},xVect);
        maxAbsMonomial = max(1, maxAbsMonomial);
        value = funcValues(1)-sqrt(sum(funcValues(2:vecSize).^(2* ...
            ones(1,vecSize-1))));
        %fprintf('%d --- %f\n',i,value);
        scaledError = min(scaledError,value/maxAbsMonomial);
        infeasError = min(infeasError,value);
    elseif inEqPolySys{i}.typeCone == 3
        [funcValues,maxAbsMonomial] = evalPolynomials(inEqPolySys{i},xVect);
        maxAbsMonomial = max(1, maxAbsMonomial);
        matSize = inEqPolySys{i}.sizeCone;
        sMat = reshape(funcValues,matSize,matSize);
        d = eig(sMat);
        scaledError = min(scaledError,min(d)/maxAbsMonomial);
        infeasError = min(infeasError,min(d));
    elseif inEqPolySys{i}.typeCone == -1
        [value,maxAbsMonomial] = evalPolynomials(inEqPolySys{i},xVect);
        maxAbsMonomial = max(1, maxAbsMonomial);
        %fprintf('%d --- %f\n',i,value);
        scaledError = min(scaledError,-abs(value)/maxAbsMonomial);
        infeasError = min(infeasError,-abs(value));
        %fprintf('%d --- %f\n',i,scaledError);
    end
end
return

% $Header: /home/waki9/CVS_DB/SparsePOPdev/subPrograms/genApproxSolution.m,v 1.1.1.1 2007/01/11 11:31:50 waki9 Exp $
