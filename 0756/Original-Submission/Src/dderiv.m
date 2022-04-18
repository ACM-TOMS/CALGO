function fprime = dderiv(zp,z,beta)
%DDERIV Derivative of the disk map.
%	DDERIV(ZP,Z,BETA) returns the derivative at the points of ZP of
%	the Schwarz-Christoffel disk map whose prevertices are Z and
%	whose turning angles are BETA.
%
%       Don't forget the multiplicative constant in the SC map!
%
%	See also DPARAM, DMAP.
%	
%	Written by Toby Driscoll.  Last updated 5/24/95.

z = z(:);
beta = beta(:);
zprow = zp(:).';
fprime = zeros(size(zp));
npts = length(zp(:));
terms = 1 - zprow(ones(length(beta),1),:)./z(:,ones(npts,1));
fprime(:) = exp(sum(log(terms).*beta(:,ones(npts,1))));
