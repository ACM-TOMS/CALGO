function wp = stmap(zp,w,beta,z,c,qdat)
%STMAP  Schwarz-Christoffel strip map.
%       STMAP(ZP,W,BETA,Z,C,QDAT) computes the values of the Schwarz-
%       Christoffel strip map at the points in vector ZP.  The arguments 
%       W, BETA, Z, C, and QDAT are as in STPARAM.  STMAP returns a vector
%       the same size as ZP.
%
%       STMAP(ZP,W,BETA,Z,C,TOL) uses quadrature data intended to give an
%       answer accurate to within roughly TOL.
%	
%	STMAP(ZP,W,BETA,Z,C) uses a tolerance of 1e-8.
%
%	See also STPARAM, STPLOT, STINVMAP.
%
%	Written by Toby Driscoll.  Last updated 5/31/95.

N = length(w);
n = N-2;
w = w(:);
beta = beta(:);
z = z(:);
% Renumber vertices so that the ends of the strip map to w([1,k])
wend = [find(isinf(z)&(z<0)),find(isinf(z)&(z>0))];
renum = [wend(1):N,1:wend(1)-1];
w = w(renum);
beta = beta(renum);
z = z(renum);
k = find(renum==wend(2));
% nb = Number of prevertices on bottom edge of strip
nb = k-2;
z([1,k]) = [];
w([1,k]) = [];

if nargin < 6
  qdat = scqdata(beta([2:k-1,k+1:N]),8);
elseif length(qdat)==1
  qdat = scqdata(beta([2:k-1,k+1:N]),max(ceil(-log10(qdat)),8));
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
  % Can't integrate starting at pre-infinity: which prevertex
  % is next closest?
  zf = z(~isinf(w));
  [tmp,s2] = min(abs(zp(bad,ones(n-ninf,1)).'-zf(:,ones(1,sum(bad)))));
  shift = cumsum(isinf(w));
  shift(atinf) = [];
  sing(bad) = s2(:) + shift(s2(:));
  % Midpoints of these integrations 
  mid = (z(sing(bad)) + zp(bad)) / 2;
else
  bad = zeros(p,1);		% all clear
end

% zs = the starting singularities
% A MATLAB technicality could cause a mistake if sing is all ones and same
% length as z, hence a workaround.
zs = wp(:);
zs(1:p+1) = z([sing;2]);
zs = zs(1:p);
% ws = map(zs)
ws = wp(:);
ws(1:p+1) = w([sing;2]);
ws = ws(1:p);

% Compute the map directly at "normal" points.
if any(~bad)
  wp(~bad) = ws(~bad) + ...
      c*stquad(zs(~bad),zp(~bad),sing(~bad),z,beta,qdat);
end
% Compute map at "bad" points, stopping at midpoint to avoid integration
% where right endpoint is close to a singularity.
if any(bad)
  wp(bad) = ws(bad) + c*...
      (stquad(zs(bad),mid,sing(bad),z,beta,qdat) -...
      stquad(zp(bad),mid,zeros(sum(bad),1),z,beta,qdat));
end


