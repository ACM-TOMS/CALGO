function [POP, SDPobjValue, cpuTime, SeDuMiInfo,infeasibleSW] ... 
    = nDimZero2(objPoly, objPoly0, ineqPolySys0, lbd0, ubd0,POP, ... 
    SDPobjValue, cpuTime, startingTime1, SeDuMiInfo, fixedVar)

% If all variables vanish from the original POP by "deleteVar.m", 
% SDP relaxations are not applied into the POP.

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
cpuTime.conversion = cputime - startingTime1;
cpuTime.total = cpuTime.conversion;
SeDuMiInfo.numerr = 0;
SeDuMiInfo.pinf   = 0;
SeDuMiInfo.dinf   = 0;
POP.xVect = zeros(objPoly0.dimVar,1);
POP.xVect(fixedVar(:,1)) = fixedVar(:,2);
POP.objValue = evalPolynomials(objPoly0,POP.xVect);
SDPobjValue = POP.objValue;
[ineqPolySys,lbdIdx,ubdIdx] = boundToIneqPolySys(ineqPolySys0,lbd0,ubd0); 
[POP.absError, POP.scaledError] = infeasibility(ineqPolySys,POP.xVect); 
if abs(POP.absError) > 1.0e-6 || abs(POP.scaledError) > 1.0e-6
%		fprintf('## POP is infeasible.                               ##\n\n');
    fprintf('## POP is infeasible.\n');
	infeasibleSW = 2;
else
    fprintf('## POP has a unique feasible solution.\n');
	infeasibleSW = -2;
end
return

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
