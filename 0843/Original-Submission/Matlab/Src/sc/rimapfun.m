function zdot = rimapfun(wp,yp,flag,scale,z,beta,c,zs,L);
%   Used by RINVMAP for solution of an ODE.

%   Copyright 1998 by Toby Driscoll.
%   $Id: rimapfun.m,v 2.1 1998/05/10 04:51:53 tad Exp $

lenyp = length(yp);
lenzp = lenyp/2;
zp = yp(1:lenzp) + i*yp(lenzp+1:lenyp);

f = scale./rderiv(zp,z,beta,c,L,zs);
zdot = [real(f);imag(f)];
