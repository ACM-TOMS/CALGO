function fprime = rderiv(zp,z,beta,L,zs)
%RDERIV Derivative of the rectangle map.
%	RDERIV(ZP,Z,BETA,L) returns the derivative at the points of ZP of
%	the Schwarz-Christoffel rectangle map whose prevertices are Z,
%	turning angles are BETA, and aspect ratio parameter is L.
%	
%	If a fifth argument is supplied, it is assumed to be the image
%	of Z on the intermediate strip; see R2STRIP.
%	
%	Don't forget the multiplicative constant in the SC map!
%	
%	See also RPARAM, RMAP, R2STRIP.
%
%	Written by Toby Driscoll.  Last updated 5/24/95.

n = length(z);

if nargin < 5
  % Find prevertices on the strip
  zs = r2strip(z,z,L);
  zs = real(zs) + i*round(imag(zs)); 	% put them *exactly* on edges
end

% First compute map and derivative from rectangle to strip
[F,dF] = r2strip(zp,z,L);

% Now compute derivative of map from strip to polygon
[tmp,j1] = min(zs);
renum = [j1:n,1:j1-1];
zs = zs(renum);
beta = beta(renum);
nb = sum(~imag(zs));
zs = zs(:);
zs = [-Inf; zs(1:nb); Inf; zs(nb+1:n)];
betas = [0; beta(1:nb); 0; beta(nb+1:n)];
dG = stderiv(F,zs,betas);

% Put it together
fprime = dF.*dG;

