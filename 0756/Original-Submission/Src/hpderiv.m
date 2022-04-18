function fprime = hpderiv(zp,x,beta)
%HPDERIV Derivative of the half-plane map.
%	HPDERIV(ZP,X,BETA) returns the derivative at the points of ZP of
%	the Schwarz-Christoffel half-plane map whose prevertices are X and
%	whose turning angles are BETA.
%
%       Don't forget the multiplicative constant in the SC map!
%	
%	See also HPPARAM, HPMAP.
%	
%	Written by Toby Driscoll.  Last updated 5/24/95.

x = x(:);
beta = beta(:);
zprow = zp(:).';
fprime = zeros(size(zp));
npts = length(zp(:));
terms = zprow(ones(length(beta),1),:) - x(:,ones(npts,1));
fprime(:) = exp(sum(log(terms).*beta(:,ones(npts,1))));
