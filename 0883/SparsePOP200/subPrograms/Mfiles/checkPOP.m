function continueSW = checkPOP(objPoly,ineqPolySys,lbd,ubd,param)

% This function checks the input data of POP described in terms
% of objPoly,ineqPolySys,lbd, ubd, and their consistency with param.
% Specifically this function checks the followings:
% 1. A common dimension nDim of the variable x over all
%	polynomials involved in objPoly and ineqPolySys.
% 2. noTerm = #Support = #coef?
% 3. degree of each polynomial.
% 4. consistensy of typeCone with sizeCone.
% 5. supports has multiple identical row vectors.


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
mDim = size(ineqPolySys,2);

continueSW = 1;

% Check objPoly
[continueSW] = checkPolynomial(objPoly,nDim,0,continueSW);
for i=1:mDim
    % Check ineqPolySys{i}
    [continueSW] = checkPolynomial(ineqPolySys{i},nDim,i,continueSW);
end
% Check lbd
lenLbd = length(lbd);
if lenLbd ~= nDim
    fprintf('## the length of lbd=%d ~= objPoly.dimVar=%d ##\n',lenLbd,nDim);
    continueSW = 0;
end
% Check ubd
lenUbd = length(ubd);
if lenUbd ~= nDim
    fprintf('## the length of ubd=%d ~= objPoly.dimVar=%d ##\n',lenUbd,nDim);
    continueSW = 0;
end
% Check param
if isfield(param,'boundSW')
    if (param.boundSW ~= 0) && (param.boundSW ~= 1) && (param.boundSW ~= 2)
        fprintf('## param.boundSW=%d needs to be 0, 1 or 2 ##\n',param.boundSW);
        continueSW = 0;
    end
end
if isfield(param,'complementaritySW')
    if (param.complementaritySW ~= 0) && (param.complementaritySW ~= 1)
        fprintf('## param.complementaritySW=%d needs to be 0 or 1 ##\n',param.complementaritySW);
        continueSW = 0;
    end
end
if isfield(param,'eqTolerance')
    if (0 < param.eqTolerance) && (param.eqTolerance < 0.9e-10 )
        fprintf('## param.eqTolerance=%6.2e needs to be either 0.0 or not less than 1.0e-10 ##\n',...
            param.eqTolerance);
        continueSW = 0;
    end
    %%param.multiCliquesFactor = 1;
end
if isfield(param,'multiCliquesFactor')
    if (param.multiCliquesFactor < 0.0) || (nDim < param.multiCliquesFactor)
        fprintf('## 0 <= param.multiCliquesFactor=%6.2e <= objPoly.dimVar=%d ##\n',...
            param.multiCliquesFactor,nDim);
        continueSW = 0;
    end
end
if isfield(param,'perturbation')
    if (0.0 < abs(param.perturbation)) && (abs(param.perturbation) < 0.9e-10 )
        fprintf('## abs(param.perturbation)=%6.2e needs to be either 0.0 or not less than 1.0e-10 ##\n',...
            abs(param.perturbation));
        continueSW = 0;
    end
end
if isfield(param,'printLevel')
    tf = ismember(param.printLevel',[0;1;2],'rows');
    if ~all(tf)
        fprintf('## param.printLevel=%d needs to be 0, 1 or 2 ##\n',param.printLevel);
        continueSW = 0;
    end
end
if isfield(param,'reduceMomentMatSW')
    if (param.reduceMomentMatSW ~= 0) && (param.reduceMomentMatSW ~= 1)
        fprintf('## param.reduceMomentMatSW=%d needs to be 0 or 1 ##\n',param.reduceMomentMatSW);
        continueSW = 0;
    end
end
%{
if isfield(param,'relaxOrder')
    maxDegree = objPoly.degree;
    m = size(ineqPolySys,2);
    for i=1:m
        maxDegree = max(maxDegree,ineqPolySys{i}.degree);
    end
    if param.relaxOrder < ceil(maxDegree/2);
        fprintf('## param.relaxOrder=%d < ceil(maxDegree/2)=%d ##\n',param.relaxOrder,ceil(maxDegree/2));
        continueSW = 0;
    end
end
%}
if isfield(param,'scalingSW')
    if (param.scalingSW ~= 0) && (param.scalingSW ~= 1)
        fprintf('## param.scalingSW=%d needs to be 0 or 1 ##\n',param.scalingSW);
        continueSW = 0;
    end
end
if isfield(param,'sdpaOutputFile')
    if ischar(param.sdpaOutputFile) == 0;
        fprintf('## param.sdpaOutputFile needs to a file name ##\n');
        continueSW = 0;
    end
end
if isfield(param,'SeDuMiSW')
    if (param.SeDuMiSW ~= 1) && (param.SeDuMiSW ~= 0);
        fprintf('## param.SeDuMiSW=%d needs to be 0 or 1 ##\n',param.SeDuMiSW);
        continueSW = 0;
    end
end
if isfield(param,'sparseSW')
    if (param.sparseSW ~= 0) && (param.sparseSW ~= 1)
        fprintf('## param.sparseSW=%d needs to be 0 or 1 ##\n',param.sparseSW);
        continueSW = 0;
    end
end

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [continueSW] = checkPolynomial(poly,nDim,i,continueSW)
if (i == 0)
    % Check typeCone of objPoly
    if (abs(poly.typeCone) > 1)
        fprintf('## Only typeCone = -1, 0 or 1 is possible in objPoly ##\n');
        continueSW = 0;
    end
    % Check sizeCone of objPoly
    % Kojima, 02/15/05
    % Checking of sizeCone only when sizeCone is specified; --->
    if (i == 0) && (isfield(poly,'sizeCone')) && (poly.sizeCone ~= 1)
        fprintf('## Only sizeCone = 1 is possible in objPoly ##\n');
        continueSW = 0;
    end
    % <--- Checking of sizeCone only when sizeCone is specified;
else % i >= 1
    % Check typeCone and sizeCone of ineqPolysys{i}
    if (poly.typeCone == -1) || (poly.typeCone == 1) % real equality or inequality
        % Kojima, 02/15/05
        % Checking of sizeCone only when sizeCone is specified; --->
        if (isfield(poly,'sizeCone')) && (poly.sizeCone ~= 1)
            fprintf('## Only sizeCone = 1 is possible when sizeCone = -1 or 1 in ineqPolySys{%d} ##\n',i);
            continueSW = 0;
        end
        % <--- Checking of sizeCone only when sizeCone is specified;
    elseif poly.typeCone == 2 % SOCP cone
        if poly.sizeCone <= 1
            fprintf('## sizeCone=%d >=2 when typeCone = 2 in ineqPolySys{%d} ##\n',poly.sizeCone,i);
            continueSW = 0;
        end
    elseif poly.typeCone == 3 %SDP cone
        if poly.sizeCone <= 1
            fprintf('## sizeCone=%d >=2 when typeCone = 3 in ineqPolySys{%d} ##\n',poly.sizeCone,i);
            continueSW = 0;
        end
    else
        fprintf('## typeCone=%d needs to be either -1, 1, 2 or 3 in ineqPolySys{%d} ##\n',poly.tyoeCone,i);
        continueSW = 0;
    end
end
% Check dimVar.
if (0 < i) && (poly.dimVar ~= nDim)
    fprintf('## ineqPolySys{%d}.dimVar=%d ~= objPoly.dimVar=%d ##\n',i,poly.dimVar,nDim);
    continueSW = 0;
end
% Check the sizes of supports,
rowSizeS = size(poly.supports,1);
colSizeS = size(poly.supports,2);
if poly.noTerms ~= rowSizeS
    if i == 0
        fprintf('## objPoly.noTerms=%d ~= size(objPoly.supports,1)=%d ##\n',...
            poly.noTerms,size(poly.supports,1));
    else
        fprintf('## ineqPolySys{%d}.noTerms=%d ~= size(ineqPolySys{%d}.supports,1)=%d ##\n',...
            i,poly.noTerms,i,size(poly.supports,1));
    end
    continueSW = 0;
end
if poly.dimVar ~= colSizeS
    if i == 0
        fprintf('## objPoly.dimVar=%d ~= size(objPoly.supports,2)=%d ##\n',nDim,size(poly.supports,2));
    else
        fprintf('## ineqPolySys{%d}.dimVar=%d ~= size(ineqPolySys{%d}.supports,2)=%d ##\n',...
            i,poly.dimVar,i,size(poly.supports,2));
    end
    continueSW = 0;
end
% Check the sizes of coef.
rowSizeC = size(poly.coef,1);
colSizeC = size(poly.coef,2);
if poly.noTerms ~= rowSizeC
    if i == 0
        fprintf('## objPoly.noTerms=%d ~= size(objPoly.coef,1)=%d ##\n',...
            poly.noTerms,size(poly.coef,1));
    else
        fprintf('## ineqPolySys{%d}.noTerms=%d ~= size(ineqPolySys{%d}.coef,1)=%d ##\n',...
            i,poly.noTerms,i,size(poly.coef,1));
    end
    continueSW = 0;
end
if abs(poly.typeCone) <= 2
    if colSizeC ~= poly.sizeCone
        if i == 0
            fprintf('## objPoly.coef ~= noTerms x objPoly.sizeCone  ##\n');
        else
            fprintf('## ineqPolySys{%d}.coef ~= noTerms x ineqPolySys{%d}.sizeCone ##\n',i,i);
        end
        continueSW = 0;
    end
elseif (i > 0) && (poly.typeCone == 3)
    if colSizeC ~= poly.sizeCone * poly.sizeCone
        fprintf('## ineqPolySys{%d}.coef ~= noTerms x (ineqPolySys{%d}.sizeCone*ineqPolySys{%d}.sizeCone) ##\n',i,i,i);
        continueSW = 0;
    end
end
% Check whether the supports are nonnegative integers
indexSet = find(poly.supports < 0);
if ~isempty(indexSet)
    if i == 0
        fprintf('## Some element objPoly.supports is negative   ##\n');
    else
        fprintf('## Some element of ineqPolySys{%d}.supports is negative ##\n',i);
    end
    continueSW = 0;
end
indexSet = find(sum(poly.supports - ceil(poly.supports),2));
if ~isempty(indexSet)
    if i == 0
        fprintf('## Some element objPoly.supports is not integer ##\n');
    else
        fprintf('## Some element of ineqPolySys{%d}.supports is not integer ##\n',i);
    end
    continueSW = 0;
end
% Check degree.
degree = max(sum(poly.supports,2));
if degree ~= poly.degree
    if i == 0
        fprintf('## objPoly.degree is different from max(sum(objPoly.supports,2)) ##\n');
    else
        fprintf('## ineqPolySys{%d}.degree is different from max(sum(ineqPolySys{%d}.supports,2)) ##\n',i,i);
    end
    continueSW = 0;
end
%% Check more than one identical supports
if exist('mexconv3') == 3
	%disp(full(supmat));	
	[temp, M, N] = quickUnique(poly.supports);
	%disp(full(temp));
	%disp(length(M));
	%disp(length(N));	
else
	[temp, M, N] = unique(poly.supports,'rows');
end
if length(N) > length(M)
    if i == 0
        fprintf('## objPoly involves more than one idientical support row vectors in objPoly.supports ##\n');
    else
        fprintf('## ineqPolySys{%d} involves more than one idientical support row vectors in ineqPolySys{%d}.supports ##\n',i,i);
    end
	%full([poly.supports])
	%full([poly.coef])
    continueSW = 0;
end

return;

% $Header: /home/waki9/CVS_DB/SparsePOPdev/subPrograms/checkPOP.m,v 1.1.1.1 2007/01/11 11:31:50 waki9 Exp $
