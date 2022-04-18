function wp = eval(M,zp,tol)
%EVAL Evaluate Schwarz-Christoffel rectangle map at points.
%   EVAL(M,ZP) evaluates the Schwarz-Christoffel map M at the points ZP
%   in the source rectangle of M. The default tolerance of M is used.
%   
%   EVAL(M,ZP,TOL) attempts to give an answer accurate to TOL. If TOL is
%   less than the accuracy of M, this is unlikely to be met.
%   
%   See also RECTMAP, EVALINV.

%   Copyright 1998 by Toby Driscoll.
%   $Id: eval.m,v 2.2 2003/01/10 15:19:42 driscoll Exp $

p = polygon(M);
n = length(p);
z = M.prevertex;
zr = z(corners(M));

if nargin < 3
  qdata = M.qdata;
else 
  qdata = tol;
end

wp = NaN*zp;
idx = find(isinpoly(zp,polygon(zr)));
wp(idx) = ...
    rmap(zp(idx),vertex(p),angle(p)-1,z,M.constant,M.stripL,qdata);
