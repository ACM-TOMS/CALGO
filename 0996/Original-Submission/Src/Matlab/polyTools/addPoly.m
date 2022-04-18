function newPoly = addPoly( poly1, poly2, sorted )
%ADDPOLY    Add two polynomials.
%  Usage:
%    newPoly = addPoly(poly1, poly2);
%    newPoly = addPoly(poly1, poly2, false); 
%
%    The latter is faster as it does not sort the supports.
%
%    newPoly, poly1, and poly2 are polynomials in the simplified sparsePOP format.
%
%    If poly1 (or poly2) is empty [], 
%    then it returns poly2 (or poly1, respectively).

if nargin == 2
    sorted = true;
end
if ~( (isempty(poly1)||checkPoly(poly1)) ...
        && (isempty(poly2)||checkPoly(poly2)) )
    error('Input format is not correct')
end
if isempty(poly1); newPoly = poly2; return; end;
if isempty(poly2); newPoly = poly1; return; end;

%[newPoly.supports,~,ic] = grevlexUnique([poly1.supports; poly2.supports]);
if sorted
    [newPoly.supports, ~, ic] = grevlexUnique([poly1.supports; poly2.supports]);
else
    [newPoly.supports, ~, ic] = fastUnique([poly1.supports; poly2.supports]);
end
newPoly.coef = accumarray(ic,[poly1.coef; poly2.coef]); 


end

