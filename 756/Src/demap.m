function wp = demap(zp,w,beta,z,c,qdat)
%DEMAP  Schwarz-Christoffel exterior map.
%       DEMAP(ZP,W,BETA,Z,C,QDAT) computes the values of the Schwarz-
%       Christoffel exterior map at the points in vector ZP.  The arguments 
%       W, BETA, Z, C, and QDAT are as in DEPARAM.  DEMAP returns a vector
%       the same size as ZP.
%
%       DEMAP(ZP,W,BETA,Z,C,TOL) uses quadrature data intended to give an
%       answer accurate to within roughly TOL.
%
%	See also DEPARAM, DEPLOT, DEINVMAP.
%
%	Written by Toby Driscoll.  Last updated 5/31/95.

if nargin < 6
  qdat = scqdata(beta,8);
elseif length(qdat)==1
  qdat = scqdata(beta,max(ceil(-log10(qdat)),8));
end

n = length(w);
beta = beta(:);
z = z(:);
p = length(zp);
wp = zeros(size(zp));
ws = wp;
zs = wp;

% For each point in zp, find nearest prevertex.
[mindist,sing] = min(abs(ones(n,1)*zp(:).'-z(:,ones(1,p))));

% zs = the starting singularities
% A MATLAB technicality could cause a mistake if sing is all ones and same
% length as z, hence a workaround.
zs(1:p+1) = z([sing,2]);
zs = zs(1:p);
% ws = SCmap(zs)
ws(1:p+1) = w([sing,2]);
ws = ws(1:p);

% Must be careful about the singularity at the origin, since the
% quadrature routine doesn't pay attention to the right endpoint.

abszp = abs(zp); 			% dist to sing at 0
zp2zs = abs(zp-zs);			% dist from zp to zs
bad = zp2zs < 10*eps;
unf = ones(size(zp2zs));		% unfinished?
dist = unf;
znew = unf;
% Take care of "bad" ones explicitly.
wp(bad) = ws(bad);
unf(bad) = zeros(size(unf(bad)));
% Integrate for the rest.
dist(unf) = min(1,2*abszp(unf)./zp2zs(unf));	% how far may we go?
znew(unf) = zs(unf) + dist(unf).*(zp(unf)-zs(unf));
wp(unf) = ws(unf) + c*dequad(zs(unf),znew(unf),sing(unf),z,beta,qdat);
unf = (dist<1); 			% unfinished positions
while any(unf)
  zold = znew;
  dist(unf) = min(1,2*abszp(unf)./abs(zp(unf)-zold(unf)));
  znew(unf) = zold(unf) + dist(unf).*(zp(unf)-zold(unf));
  wp(unf) = wp(unf) + c*dequad(zold(unf),znew(unf),[],z,beta,qdat);
  unf = (dist<1);
end

