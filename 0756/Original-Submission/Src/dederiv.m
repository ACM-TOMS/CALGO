function fprime = dederiv(zp,z,beta)
%DEDERIV Derivative of the exterior map.
%	DEDERIV(ZP,Z,BETA) returns the derivative at the points of ZP of
%	the Schwarz-Christoffel exterior map whose prevertices are Z and
%	whose turning angles are BETA.
%
%       Don't forget the multiplicative constant in the SC map!
%	
%	See also DEPARAM, DEMAP.
%	
%	Written by Toby Driscoll.  Last updated 5/24/95.

z = z(:);
beta = [beta(:);-2];
zprow = zp(:).';
fprime = zeros(size(zp));
npts = length(zp(:));
terms = 1 - zprow(ones(length(z),1),:)./z(:,ones(npts,1));
terms(length(z)+1,:) = zprow;
fprime(:) = exp(sum(log(terms).*beta(:,ones(npts,1))));
