function poly = penalyzePoly(poly1, poly2, lambda)
%%PENALYZEPOLY Penalize a polynomial.
% Usage:
%   penalizedObjPoly = penalyzePoly(objPoly, poly, lambda);
% 
% The relation of the inputs and output is as follows:
%   penalizedObjPoly = objPoly + lambda * (poly)^2
    tmp = powPoly(poly2, 2);
    tmp.coef = tmp.coef * lambda;
    poly = addPoly(poly1, tmp);
end