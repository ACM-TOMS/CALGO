function wp = rmap(zp,w,beta,z,c,L,qdat)
%RMAP   Schwarz-Christoffel rectangle map.
%       RMAP(ZP,W,BETA,Z,C,L,QDAT) computes the values of the
%       Schwarz-Christoffel rectangle map at the points in vector ZP.
%       The remaining arguments are as in RPARAM.  RMAP returns a vector
%       the same size as ZP.
%	
%       RMAP(ZP,W,BETA,Z,C,L,TOL) uses quadrature data intended to give
%       an answer accurate to within TOL.
%	
%	RMAP(ZP,W,BETA,Z,C,L) uses a tolerance of 1e-8.
%
%	See also RPARAM, RPLOT, RINVMAP.
%
%	Written by Toby Driscoll.  Last updated 5/24/95.

n = length(w);
wp = z;
w = w(:);
beta = beta(:);
z = z(:);
[w,beta,z,corners] = rcorners(w,beta,z);

% Map prevertices to strip
K = max(real(z));
Kp = max(imag(z));
zs = r2strip(z,z,L);
zs = real(zs) + i*round(imag(zs));	% put them *exactly* on edges

if nargin < 7
  qdat = scqdata(beta,8);
elseif length(qdat)==1
  qdat = scqdata(beta,max(ceil(-log10(qdat)),8));
end
wp = zeros(size(zp));
zp = zp(:);
p = length(zp);

% Trap points which map to +/-Inf on the strip.
bad = abs(zp) < 2*eps;
zp(bad) = zp(bad) + 100*eps;
bad = abs(zp-i*Kp) < 2*eps;
zp(bad) = zp(bad) - i*100*eps*Kp;

% Map from rectangle to strip.
yp = r2strip(zp,z,L);

% Now map from strip to polygon.
i1 = 1:corners(3)-1;
i2 = corners(3):n;
ws = [NaN; w(i1); NaN; w(i2)];
bs = [0; beta(i1); 0; beta(i2)];
zs = [Inf; zs(i1); Inf; zs(i2)];
wp(:) = stmap(yp,ws,bs,zs,c,qdat);


