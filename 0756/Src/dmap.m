function wp = dmap(zp,w,beta,z,c,qdat)
%DMAP  Schwarz-Christoffel disk map.
%       DMAP(ZP,W,BETA,Z,C,QDAT) computes the values of the Schwarz-
%       Christoffel disk map at the points in vector ZP.  The arguments 
%       W, BETA, Z, C, and QDAT are as in DPARAM.  DMAP returns a vector
%       the same size as ZP.
%
%       DMAP(ZP,W,BETA,Z,C,TOL) uses quadrature data intended to give an
%       answer accurate to within roughly TOL.
%	
%	DMAP(ZP,W,BETA,Z,C) uses a tolerance of 1e-8.
%
%	See also DPARAM, DPLOT, DINVMAP.
%
%	Written by Toby Driscoll.  Last updated 5/31/95.

n = length(z);
w = w(:);
beta = beta(:);
z = z(:);
if nargin < 6
  qdat = scqdata(beta,8);
elseif length(qdat)==1
  qdat = scqdata(beta,max(ceil(-log10(qdat)),8));
end
wp = zeros(size(zp));
zp = zp(:);
p = length(zp);

% For each point in zp, find nearest prevertex.
[tmp,sing] = min(abs(zp(:,ones(n,1)).'-z(:,ones(1,p))));
sing = sing(:);				% indices of prevertices
atinf = find(isinf(w)); 		% infinite vertices 
atinf = atinf(:);
ninf = length(atinf);			% # of inf vertices
if ninf > 0
  % "Bad" points are closest to a prevertex of infinity.
  bad = sing(:,ones(ninf,1))' == atinf(:,ones(1,p));
  % Can be closest to any pre-infinity.
  if ninf > 1
    bad = any(bad);
  end
  % Exclude cases which are exactly those prevertices.
  bad = bad(:) & (abs(zp-z(sing)) > 10*eps);
  % Can't integrate starting at pre-infinity: find conformal center to use
  % as integration basis.
  if ~isinf(w(n-1))
    wc = w(n-1) + c*dquad(z(n-1),0,n-1,z,beta,qdat);
  else
    wc = w(n) + c*dquad(z(n),0,n,z,beta,qdat);
  end
else
  bad = zeros(p,1);		% all clear
  wc = [];				% don't need it
end

% zs = the starting singularities
% A MATLAB technicality could cause a mistake if sing is all ones and same
% length as z, hence a workaround.
zs = wp(:);
zs(1:p+1) = z([sing;2]);
zs = zs(1:p);
% ws = SCmap(zs)
ws = wp(:);
ws(1:p+1) = w([sing;2]);
ws = ws(1:p);

% Compute the map directly at "normal" points.
wp(~bad) = ws(~bad) + c*dquad(zs(~bad),zp(~bad),sing(~bad),z,beta,qdat);
% Compute map at "bad" points, using conformal center as basis, to avoid
% integration where right endpoint is too close to a singularity.
wp(bad) = wc - c*dquad(zp(bad),zeros(sum(bad),1),zeros(sum(bad),1),...
    z,beta,qdat);


