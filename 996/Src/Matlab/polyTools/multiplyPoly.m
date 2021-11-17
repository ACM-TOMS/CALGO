function newPoly = multiplyPoly(poly1, poly2, sorted)
%%MULTIPLYPOLY   Multiply two polynomials.
%  newPoly = multiplyPoly(poly1, poly2);
%  newPoly = multiplyPoly(poly1, poly2, false);
%
%  The latter is faster as it does not sort the supports.
%
%  newPoly, poly1, and poly2 are polynomials in the simplified sparsePOP format.
%
%  If poly1 (or poly2) is empty [], 
%  then it returns poly2 (or poly1, respectively).

if nargin == 2
    sorted = true;
end
if ~( (isempty(poly1)||checkPoly(poly1)) ...
        && (isempty(poly2)||checkPoly(poly2)) )
    error('Input format is not correct')
end
if isempty(poly1); newPoly = poly2; return; end;
if isempty(poly2); newPoly = poly1; return; end;

[idx1,idx2] = myMeshGrid(1:size(poly1.supports,1),1:size(poly2.supports,1));
newPoly.supports = poly1.supports(idx1,:) + poly2.supports(idx2,:);
newPoly.coef = poly1.coef(idx1) .* poly2.coef(idx2);

% [newPoly.supports,~,ic] = fastUnique(newPoly.supports); %grevlexUnique(newPoly.supports);
% newPoly.coef = accumarray(ic,newPoly.coef);

newPoly = simplifyPoly(newPoly, [], [], sorted);
end

