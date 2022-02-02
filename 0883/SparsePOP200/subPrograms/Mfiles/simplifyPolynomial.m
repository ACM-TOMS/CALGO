function polyOut = simplifyPolynomial(polyIn, order);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Simplify a polynomial by combining terms with a common support %
% and eliminating terms with zero coefficients.                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   polynomial.typeCone = -1, 1, 2 0r 3
%   polynomial.sizeCone >= 1
% 	polynomial.dimVar --- a positive integer 
% 	polynomial.degree --- a positive integer 
% 	polynomial.noTerms --- a positive integer 
% 	polynomial.supports --- a noTerms \times nDim matrix
% 	polynomial.coef --- a noTerms \times noTerms vector
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
if nargin <2
  order = 'grevlex';
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Only order = 'grevlex' is available. %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
order = 'grevlex';

n = size(polyIn.supports,2);
col = size(polyIn.coef,2);

UsedVariable = find(any(polyIn.supports,1));
% [tmpSupSet,I] = monomialSort(polyIn.supports,UsedVariable,'grevlex');
[tmpSupSet,I] = monomialSort(polyIn.supports,UsedVariable,order);
tmpCoef = polyIn.coef(I,:);
tmpSupSet = [sum(tmpSupSet(:,UsedVariable),2),-tmpSupSet(:,UsedVariable)];
[tmpSupSet,M,N] = unique(tmpSupSet,'rows');
m = size(tmpSupSet,1);
SupSet = sparse(m,n);
SupSet(:,UsedVariable) = -tmpSupSet(:,2:end);
Coef = sparse(length(M),col);
shift_M = [0;M];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Add up coefficients of common supports. %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for j=1:length(M)
    Coef(j,:) = sum(tmpCoef(shift_M(j)+1:shift_M(j+1),:),1);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Remove support and coef such that abs(coef) = 0. % 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
CoefMat = (abs(Coef) > 1.0e-9);
NonZero = any(CoefMat,2);
NonZeroIndex = find(NonZero);
newCoef = Coef(NonZeroIndex,:);

% Kojima 02/17/2006
% For the case that NonZeroIndex is empty set ---> 
if isempty(NonZeroIndex)
    polyOut.typeCone = polyIn.typeCone;
	polyOut.sizeCone = polyIn.sizeCone;
    polyOut.dimVar   = polyIn.dimVar; 
    polyOut.degree   = 0; 
    polyOut.noTerms  = 1; 
    polyOut.supports = sparse(zeros(1,polyIn.dimVar));
	polyOut.coef     = sparse(zeros(1,col));
% <--- For the case that NonZeroIndex is empty set
else 
    newSupSet = sparse(length(NonZeroIndex),n);
    %
    % size(SupSet)
    % full(newSupSet)
    % NonZeroIndex
    % UsedVariable
    % size(newSupSet)
    % full(newSupSet)
    %
    % writePolynomials(1,polyIn); 
    %
    newSupSet(:,UsedVariable) = SupSet(NonZeroIndex,UsedVariable);
    maxDegree = max(sum(newSupSet(:,UsedVariable),2));
    maxDegree = full(maxDegree);
    polyOut.typeCone = polyIn.typeCone;
    polyOut.sizeCone = polyIn.sizeCone;
    polyOut.dimVar   = polyIn.dimVar; 
    if size(newCoef,1) >= 1 
        polyOut.degree   = maxDegree; 
        polyOut.noTerms  = size(newCoef,1); 
        polyOut.supports = newSupSet;
        polyOut.coef     = newCoef; 
    else 
        polyOut.noTerms = 1; 
        polyOut.degree = 0; 
        polyOut.supports = sparse(1,n); 
        polyOut.coef = sparse(1,col); 
    end
end

return;

% $Header: /home/waki9/CVS_DB/SparsePOPdev/subPrograms/simplifyPolynomial.m,v 1.1.1.1 2007/01/11 11:31:50 waki9 Exp $
