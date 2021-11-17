function wp = hpmap(zp,w,beta,x,c,qdat)
%HPMAP  Schwarz-Christoffel half-plane map.
%	HPMAP(ZP,W,BETA,X,C,QDAT) computes the values of the
%	Schwarz-Christoffel half-plane map at the points in vector ZP.  The
%	polygon's vertices should be given in W and the arguments X, C, and
%	QDAT should be computed by HPPARAM.  HPMAP returns a vector the same
%	size as ZP. 
%	
%       HPMAP(ZP,W,BETA,X,C,TOL) uses quadrature data intended to give an
%       answer accurate to within roughly TOL.
%	
%	HPMAP(ZP,W,BETA,X,C) uses a tolerance of 1e-8.
%
%	See also HPPARAM, HPPLOT, HPINVMAP.
%
%	Written by Toby Driscoll.  Last updated 5/31/95.

n = length(w);
w = w(:);
beta = beta(:);
x = x(:);
if any(isinf(x))
  x(n) = [];
  beta(n) = [];
end

if nargin < 6
  qdat = scqdata(beta,8);
elseif length(qdat)==1
  qdat = scqdata(beta,max(ceil(-log10(qdat)),8));
end
wp = zeros(size(zp));
zp = zp(:);
p = length(zp);

% For each point in zp, find nearest prevertex.
[tmp,sing] = min(abs(zp(:,ones(length(x),1)).'-x(:,ones(1,p))));
sing = sing(:);				% indices of prevertices
atinf = find(isinf(w(1:length(x)))); 	% infinite vertices
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
  bad = bad(:) & (abs(zp-x(sing)) > 10*eps);
  % Can't integrate starting at pre-infinity: which neighboring prevertex
  % to use?
  direcn = real(zp(bad)-x(sing(bad)));
  sing(bad) = sing(bad) + sign(direcn) + (direcn==0);
  % Midpoints of these integrations 
  mid = (x(sing(bad)) + zp(bad)) / 2;
else
  bad = zeros(p,1);
end
  
% xs = the starting singularities
% A MATLAB technicality could cause a mistake if sing is all ones and same
% length as x, hence a workaround.
xs = wp(:); xs(1:p+1) = x([sing;2]); xs = xs(1:p);
% ws = f(xs)
ws = wp(:); ws(1:p+1) = w([sing;2]); ws = ws(1:p);

% Compute the map directly at "normal" points.
if any(~bad)
  wp(~bad) = ws(~bad) + c*hpquad(xs(~bad),zp(~bad),sing(~bad),...
      x,beta,qdat);
end
% Compute map at "bad" points, stopping at midpoint to avoid integration
% where right endpoint is close to a singularity.
if any(bad)
  wp(bad) = ws(bad) + c*...
      (hpquad(xs(bad),mid,sing(bad),x,beta,qdat) -...
      hpquad(zp(bad),mid,zeros(sum(bad),1),x,beta,qdat));
end
  


