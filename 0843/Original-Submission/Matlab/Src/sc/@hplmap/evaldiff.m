function fp = evaldiff(M,zp)
%EVALDIFF Derivative of Schwarz-Christoffel half-plane map at points.
%   EVALDIFF(M,ZP) computes the derivative of the Schwarz-Christoffel
%   half-plane map M at the points ZP.
%   
%   See also HPLMAP, EVAL.

%   Copyright 1998 by Toby Driscoll.
%   $Id: evaldiff.m,v 2.1 1998/05/10 04:17:38 tad Exp $

z = M.prevertex;
c = M.constant;
beta = angle(polygon(M)) - 1;

fp = hpderiv(zp,z,beta,c);
