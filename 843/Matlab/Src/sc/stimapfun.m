function zdot = stimapfun(wp,yp,flag,scale,z,beta,c);

%   Used by STINVMAP for solution of an ODE.

%   Copyright 1998 by Toby Driscoll.
%   $Id: stimapfun.m,v 2.1 1998/05/10 04:54:17 tad Exp $

lenyp = length(yp);
lenzp = lenyp/2;
zp = yp(1:lenzp) + i*yp(lenzp+1:lenyp);

f = scale./stderiv(zp,z,beta,c);
zdot = [real(f);imag(f)];
