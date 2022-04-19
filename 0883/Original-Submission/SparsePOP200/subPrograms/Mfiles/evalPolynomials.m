function [value,maxAbsMonomial] = evalPolynomials(polyIn,xVect)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Modified by Kojima, 01/06/05
% (1) 	The modified version outputs the maximum absolute value,  
% 	maxAbsMonimial of monomials involved in the polynomial 
% 	evaluated at xVect. This gives information how many digits 
% 	are reliable in the evaluation of the polynomial.
% 	For example, if p(x_1,x_2) = 3*x_1^2 - 0.5*x_1*x_2^2 
% 	and (x_1,x_2) = (-2,3) then
% 		maxAbsMonomial = max{3*2^2, 0.5*2*3^2} = 12.  
% (2) 	Setting the negligiblly small number to be epsilon instead 
%	of a constant 1.0e-5 in the previous version.
%	epsilon = 1.0e-8; 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Let [dummy,noOfPolynomials] = size(polyIn); 
% If 1 = noOfPolynomials, then the polyIn is a single polynomial 
% described as
%	polyIn.typeCone = -1, 1, 2, or 3
%	polyIn.sizeCone >= 1
% 	polyIn.dimVar --- a positive integer 
% 	polyIn.degree --- a positive integer 
% 	polyIn.noTerms --- a positive integer 
% 	polyIn.supports --- a noTerms \times nDim matrix
% 	polyIn.coef --- a polyIn.noTerms \times size(polyIn.coef,2) matrix.
% In this case, evalPolynomials(polyIn,xVect) returns a real function 
% value or a row vector function value of the polynomial at xVect. 
% If 1 < noOfPolynomials, then the polyIn is a system of polynomials 
% described as 
%	polyIn{i}.typeCone = -1, 1, 2, or 3
%	polyIn{i}.sizeCone >= 1
% 	polyIn{i}.dimVar --- a positive integer 
% 	polyIn{i}.degree --- a positive integer 
% 	polyIn{i}.noTerms --- a positive integer 
% 	polyIn{i}.supports --- a noTerms \times nDim matrix
% 	polyIn{i}.coef --- a polyIn{i}.noTerms \times size(polyIn{i}.coef,2) matrix
%	(i=1,2,\ldots,noOfPolynomials). 
% In this case, evalPolynomials(polyIn,xVect) returns a column vector or a matrix 
% value of the polynomial system at xVect.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
epsilon = 0;%1.0e-8; 

nDim = size(xVect,1); 
if nDim ~= polyIn.dimVar 
	error('!!! nDim ~= polyIn.dimVar !!!');
end
colSize = size(polyIn.coef,2);
constant = sum(polyIn.supports,2);
I = find(constant == 0);
notI = find(constant > 0);
typeCone = polyIn.typeCone;
sizeCone = polyIn.sizeCone;

if ~isempty(I);
	value = sum(polyIn.coef(I,:),1);%1\times sizeCone vector 
	value = value';
	maxAbsMonomial = norm(value,inf); 
	supSet = polyIn.supports(notI,:);
	coef = polyIn.coef(notI,:);
	noTerms = length(notI); 
else
	value = zeros(colSize,1);
	maxAbsMonomial = 0.0; 
	supSet = polyIn.supports;
    coef = polyIn.coef;
    noTerms = polyIn.noTerms;
end
if noTerms > 0
    I = find(abs(xVect) > epsilon);
    notI = find(abs(xVect) <= epsilon);
    tempxVect = ones(noTerms,1)* xVect(I)';
    if ~isempty(I)
        monomial = power(tempxVect,supSet(:,I));
        monomial = prod(monomial,2);
    else
        t = size(supSet,1);
        monomial = zeros(t,1);
    end
    J = find(any(supSet(:,notI),2)>0);
    monomial(J) = 0.0;
    value = value + coef'* monomial;
else
    error('input polynomial is empty or constant.');
end

if typeCone == 3
	if sizeCone*sizeCone == colSize
        monomial = monomial(:,ones(1,colSize));
        AllElemVal = coef.*monomial;
        CandMaxVal = sum(abs(AllElemVal),2);
        maxAbsMonomial = max(maxAbsMonomial,norm(CandMaxVal,inf));
		value = reshape(value, sizeCone,sizeCone);
	else
		error('Warning: Please check sizeCone and coef size.');
	end
else
   	monomial = monomial(:,ones(1,sizeCone));
   	AllElemVal = coef.*monomial;
   	CandMaxVal = max(abs(AllElemVal),2);
   	maxAbsMonomial = max(maxAbsMonomial,norm(CandMaxVal,inf));
	if sizeCone > 1
		value = value';
	end
end

return

% $Header: /home/waki9/CVS_DB/SparsePOPdev/subPrograms/evalPolynomials.m,v 1.1.1.1 2007/01/11 11:31:50 waki9 Exp $
